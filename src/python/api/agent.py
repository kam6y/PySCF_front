"""
AI Agent API endpoints.
Handles chat interactions using simple Gemini API for molecular analysis assistance.
"""

import logging
import json
from typing import Any, Iterator

from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse

from generated_models import AgentChatRequest
from services.chat_history_service import get_chat_history_service
from services.settings_service import SettingsService

# Set up logging
logger = logging.getLogger(__name__)

# --- Input size caps ---
# Maximum length of a single chat message in characters.
# 10 000 chars is generous for interactive chat while preventing abuse.
MAX_MESSAGE_LENGTH: int = 10_000
# Maximum number of history items in a single chat request.
MAX_HISTORY_ITEMS: int = 200
# Maximum total characters across all history text parts before sending to the
# LLM.  Prevents extremely long context windows from being forwarded.
MAX_TOTAL_HISTORY_CHARS: int = 500_000
# Maximum length of individual history part text (same cap as a message).
MAX_HISTORY_PART_TEXT_LENGTH: int = MAX_MESSAGE_LENGTH

router = APIRouter(prefix="/api/agent")


def _format_sse_event(event: dict[str, Any]) -> str:
    """
    Format a Server-Sent Event message.

    Args:
        event: Event data to serialize

    Returns:
        Formatted SSE message string
    """
    return f"data: {json.dumps(event)}\n\n"


def _extract_history_role(message: Any) -> str:
    """Extract role from dict or Pydantic history message."""
    if isinstance(message, dict):
        return message.get("role", "")

    role_attr = getattr(message, "role", "")
    return (
        str(role_attr).split(".")[-1] if hasattr(role_attr, "value") else str(role_attr)
    )


def _extract_text_parts(parts: list[Any]) -> list[str]:
    """Extract text values from dict or Pydantic history parts.

    Only ``str`` values are included.  ``None`` and non-str values
    (which can arise from ``Part.text: Optional[str]`` or untyped
    dicts) are silently skipped so the caller always receives a
    ``list[str]`` that is safe to hand to the Gemini SDK.
    """
    text_parts: list[str] = []
    for part in parts:
        if isinstance(part, dict) and "text" in part:
            val = part["text"]
        elif hasattr(part, "text"):
            val = part.text
        else:
            continue
        if isinstance(val, str):
            text_parts.append(val)
    return text_parts


def _convert_history_to_gemini_format(history: list[Any]) -> list[dict[str, Any]]:
    """
    Convert frontend message history format to Gemini API format.

    Frontend format:
        [
            {"role": "user", "parts": [{"text": "Hello"}]},
            {"role": "model", "parts": [{"text": "Hi there!"}]}
        ]

    Gemini format:
        [
            {"role": "user", "parts": ["Hello"]},
            {"role": "model", "parts": ["Hi there!"]}
        ]

    Args:
        history: List of message dictionaries in frontend format

    Returns:
        List of message dictionaries in Gemini format
    """
    converted = []

    for msg in history:
        role = _extract_history_role(msg)
        parts = (
            msg.get("parts", []) if isinstance(msg, dict) else getattr(msg, "parts", [])
        ) or []
        text_parts = _extract_text_parts(parts)

        if text_parts:
            converted.append({"role": role, "parts": text_parts})

    return converted


def _validate_chat_request(request: AgentChatRequest) -> None:
    if not request.message or not request.message.strip():
        raise HTTPException(status_code=400, detail="Message cannot be empty")

    if len(request.message) > MAX_MESSAGE_LENGTH:
        raise HTTPException(
            status_code=400,
            detail=f"Message is too long (maximum {MAX_MESSAGE_LENGTH} characters)",
        )

    history = request.history or []

    # Cap on number of history items
    if len(history) > MAX_HISTORY_ITEMS:
        raise HTTPException(
            status_code=400,
            detail=f"Too many history items (maximum {MAX_HISTORY_ITEMS})",
        )

    # Cap on total characters across all history text parts
    total_chars = 0
    for msg in history:
        parts = (
            msg.get("parts", []) if isinstance(msg, dict) else getattr(msg, "parts", [])
        ) or []
        for part in parts:
            text = (
                part.get("text", "")
                if isinstance(part, dict)
                else getattr(part, "text", "")
            )
            if text:
                part_len = len(text)
                if part_len > MAX_HISTORY_PART_TEXT_LENGTH:
                    raise HTTPException(
                        status_code=400,
                        detail=f"History item text too long (maximum {MAX_HISTORY_PART_TEXT_LENGTH} characters)",
                    )
                total_chars += part_len

    if total_chars > MAX_TOTAL_HISTORY_CHARS:
        raise HTTPException(
            status_code=400,
            detail=f"Total history text too long (maximum {MAX_TOTAL_HISTORY_CHARS} characters)",
        )


def stream_chat_response(request: AgentChatRequest) -> Iterator[dict[str, Any]]:
    """
    Create event data for simple Gemini API chat responses.

    Args:
        request: Chat request body

    Yields:
        Event dictionaries to be formatted as SSE messages by the route
    """
    # NOTE: Validation is performed by the caller (chat() route handler)
    # before this generator is invoked. Do NOT re-validate here — raising
    # HTTPException inside a generator after headers are sent is fragile.

    message = request.message
    history = request.history or []
    session_id = request.session_id

    # Save user message to database if session_id is provided
    if session_id:
        try:
            chat_service = get_chat_history_service()
            chat_service.add_message(session_id, "user", message)
            logger.debug(f"Saved user message to session: {session_id}")
        except Exception as e:
            logger.warning(f"Failed to save user message to session {session_id}: {e}")

    # Accumulate AI response for saving to database
    accumulated_response = []
    db_save_successful = False
    client_aborted = False
    stream_completed = False

    try:
        logger.debug(
            "Starting Gemini chat stream (message_length=%d, session=%s)",
            len(message),
            session_id,
        )

        # Get API key from settings
        settings_service = SettingsService()
        settings = settings_service.get_settings()
        api_key = settings.get("gemini_api_key")

        if not api_key:
            error_msg = "Gemini API key is not configured. Please set it in Settings."
            logger.error(error_msg)
            yield {"type": "error", "payload": {"message": error_msg}}
            return

        # Import and configure Gemini
        import google.generativeai as genai

        genai.configure(api_key=api_key)

        # Create model with system instruction
        model = genai.GenerativeModel(
            "gemini-2.5-flash",
            system_instruction="""You are a helpful AI assistant for quantum chemistry calculations.
You can help users understand molecular structures, explain calculation results,
and provide guidance on using the PySCF quantum chemistry application.
Be concise and helpful. When discussing chemistry concepts, be accurate and educational.""",
        )

        # Convert history to Gemini format
        gemini_history = _convert_history_to_gemini_format(history)
        logger.debug(f"Converted {len(history)} history messages to Gemini format")

        # Start chat with history
        chat = model.start_chat(history=gemini_history)

        # Send agent status
        yield {
            "type": "agent_status",
            "payload": {
                "status": "responding",
                "agent": "chat",
            },
        }

        # Stream response
        response = chat.send_message(message, stream=True)

        for chunk in response:
            if chunk.text:
                accumulated_response.append(chunk.text)
                yield {"type": "chunk", "payload": {"text": chunk.text}}

        stream_completed = True
        logger.debug("Stream completed successfully")

        # Save AI response to database BEFORE sending completion event
        if session_id and accumulated_response:
            try:
                complete_response = "".join(accumulated_response)
                chat_service = get_chat_history_service()
                chat_service.add_message(session_id, "model", complete_response)
                db_save_successful = True
                logger.info(
                    f"Saved AI response to session: {session_id} (length: {len(complete_response)} chars)"
                )
            except Exception as e:
                logger.error(
                    f"Failed to save AI response to session {session_id}: {e}",
                    exc_info=True,
                )

        # Send completion event
        yield {"type": "done"}

    except GeneratorExit:
        # Client disconnected - clean up gracefully
        client_aborted = True
        logger.info("Client disconnected from SSE stream (GeneratorExit)")
        raise

    except (BrokenPipeError, ConnectionResetError) as e:
        # Client disconnected while streaming; do not persist partial responses.
        client_aborted = True
        logger.info(f"Client disconnected from SSE stream ({type(e).__name__})")
        raise

    except Exception as e:
        logger.error(f"Error during Gemini chat streaming: {e}", exc_info=True)

        try:
            yield {
                "type": "error",
                "payload": {"message": "AI chat failed. Check settings and logs."},
            }
        except (BrokenPipeError, ConnectionResetError, GeneratorExit):
            logger.debug("Unable to send error message - connection closed")

    finally:
        # Fallback: Save AI response to database if not already saved.
        # Only persist when the stream completed fully -- a mid-stream
        # exception must NOT save a truncated response as if it were complete.
        if (
            session_id
            and accumulated_response
            and not db_save_successful
            and not client_aborted
            and stream_completed
        ):
            try:
                complete_response = "".join(accumulated_response)
                chat_service = get_chat_history_service()
                chat_service.add_message(session_id, "model", complete_response)
                logger.warning(
                    f"Fallback save: AI response saved to session {session_id} after error"
                )
            except Exception as e:
                logger.error(
                    f"Fallback save failed for session {session_id}: {e}", exc_info=True
                )


@router.post("/chat")
def chat(request: AgentChatRequest) -> StreamingResponse:
    """Chat with AI agent using simple Gemini API with Server-Sent Events."""
    _validate_chat_request(request)

    session_id = request.session_id
    logger.info(
        "Processing chat request - Message length: %s, History entries: %s, Session ID: %s",
        len(request.message),
        len(request.history or []),
        session_id,
    )

    def event_generator() -> Iterator[str]:
        for event in stream_chat_response(request):
            yield _format_sse_event(event)

    return StreamingResponse(event_generator(), media_type="text/event-stream")
