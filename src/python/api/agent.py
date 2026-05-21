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

# Constants
MAX_MESSAGE_LENGTH = 100000  # Maximum allowed message length in characters
router = APIRouter(prefix='/api/agent')


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
    return str(role_attr).split('.')[-1] if hasattr(role_attr, 'value') else str(role_attr)


def _extract_text_parts(parts: list) -> list:
    """Extract text values from dict or Pydantic history parts."""
    text_parts = []
    for part in parts:
        if isinstance(part, dict) and "text" in part:
            text_parts.append(part["text"])
        elif hasattr(part, "text"):
            text_parts.append(part.text)
    return text_parts


def _convert_history_to_gemini_format(history: list) -> list:
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
        parts = msg.get("parts", []) if isinstance(msg, dict) else getattr(msg, "parts", [])
        text_parts = _extract_text_parts(parts)

        if text_parts:
            converted.append({
                "role": role,
                "parts": text_parts
            })

    return converted


def _validate_chat_request(request: AgentChatRequest) -> None:
    if not request.message or not request.message.strip():
        raise HTTPException(status_code=400, detail="Message cannot be empty")

    if len(request.message) > MAX_MESSAGE_LENGTH:
        raise HTTPException(
            status_code=400,
            detail=f"Message is too long (maximum {MAX_MESSAGE_LENGTH} characters)",
        )


def stream_chat_response(request: AgentChatRequest) -> Iterator[dict[str, Any]]:
    """
    Create event data for simple Gemini API chat responses.

    Args:
        request: Chat request body

    Yields:
        Event dictionaries to be formatted as SSE messages by the route
    """
    _validate_chat_request(request)

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

    try:
        logger.debug(f"Starting Gemini chat stream for message: {message[:100]}{'...' if len(message) > 100 else ''}")

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
            'gemini-2.5-flash',
            system_instruction="""You are a helpful AI assistant for quantum chemistry calculations.
You can help users understand molecular structures, explain calculation results,
and provide guidance on using the PySCF quantum chemistry application.
Be concise and helpful. When discussing chemistry concepts, be accurate and educational."""
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

        logger.debug("Stream completed successfully")

        # Save AI response to database BEFORE sending completion event
        if session_id and accumulated_response:
            try:
                complete_response = ''.join(accumulated_response)
                chat_service = get_chat_history_service()
                chat_service.add_message(session_id, "model", complete_response)
                db_save_successful = True
                logger.info(f"Saved AI response to session: {session_id} (length: {len(complete_response)} chars)")
            except Exception as e:
                logger.error(f"Failed to save AI response to session {session_id}: {e}", exc_info=True)

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
        accumulated_response.append(f"\n\n[Error: {str(e)}]")

        try:
            yield {
                "type": "error",
                "payload": {"message": f"An error occurred during the stream: {str(e)}"},
            }
        except (BrokenPipeError, ConnectionResetError, GeneratorExit):
            logger.debug("Unable to send error message - connection closed")

    finally:
        # Fallback: Save AI response to database if not already saved
        if session_id and accumulated_response and not db_save_successful and not client_aborted:
            try:
                complete_response = ''.join(accumulated_response)
                chat_service = get_chat_history_service()
                chat_service.add_message(session_id, "model", complete_response)
                logger.warning(f"Fallback save: AI response saved to session {session_id} after error")
            except Exception as e:
                logger.error(f"Fallback save failed for session {session_id}: {e}", exc_info=True)


@router.post('/chat')
def chat(request: AgentChatRequest):
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

    return StreamingResponse(event_generator(), media_type='text/event-stream')
