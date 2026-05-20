"""
AI Agent API endpoints.
Handles chat interactions using simple Gemini API for molecular analysis assistance.
"""

import logging
import json
from typing import Iterator, Dict, Any
from flask import Blueprint, jsonify, Response, stream_with_context
from flask_pydantic import validate

from generated_models import AgentChatRequest
from services.chat_history_service import get_chat_history_service
from services.settings_service import SettingsService

# Set up logging
logger = logging.getLogger(__name__)

# Constants
MAX_MESSAGE_LENGTH = 100000  # Maximum allowed message length in characters


def _format_sse_event(event_type: str, payload: Dict[str, Any] = None) -> str:
    """
    Format a Server-Sent Event message.

    Args:
        event_type: Type of event ('chunk', 'done', 'error')
        payload: Optional payload data

    Returns:
        Formatted SSE message string
    """
    event_data = {"type": event_type}
    if payload:
        event_data["payload"] = payload
    return f"data: {json.dumps(event_data)}\n\n"


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


def _create_simple_chat_stream(message: str, history: list, session_id: str = None) -> Iterator[str]:
    """
    Create an SSE stream for simple Gemini API chat responses.

    Args:
        message: User's message
        history: Chat history in frontend format (dict list)
        session_id: Optional session ID for persisting conversation history

    Yields:
        SSE formatted strings
    """
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
            yield _format_sse_event("error", {"message": error_msg})
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
        yield _format_sse_event("agent_status", {
            "status": "responding",
            "agent": "chat"
        })

        # Stream response
        response = chat.send_message(message, stream=True)

        for chunk in response:
            if chunk.text:
                accumulated_response.append(chunk.text)
                yield _format_sse_event("chunk", {"text": chunk.text})

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
        yield _format_sse_event("done")

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
            yield _format_sse_event("error", {"message": f"An error occurred during the stream: {str(e)}"})
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


# Create agent blueprint
agent_bp = Blueprint('agent', __name__)


@agent_bp.route('/api/agent/chat', methods=['POST'])
@validate()
def chat_with_agent(body: AgentChatRequest):
    """Chat with AI agent using simple Gemini API with Server-Sent Events."""
    try:
        # Input validation
        if not body.message or not body.message.strip():
            raise ValueError("Message cannot be empty")

        if len(body.message) > MAX_MESSAGE_LENGTH:
            raise ValueError(f"Message is too long (maximum {MAX_MESSAGE_LENGTH} characters)")

        session_id = body.session_id if hasattr(body, 'session_id') else None
        logger.info(f"Processing chat request - Message length: {len(body.message)}, History entries: {len(body.history or [])}, Session ID: {session_id}")

        # Use simple Gemini chat
        return Response(
            stream_with_context(_create_simple_chat_stream(body.message, body.history or [], session_id)),
            content_type='text/event-stream'
        )

    except ValueError as e:
        logger.warning(f"Validation error in agent chat: {e}")
        return jsonify({
            'success': False,
            'error': f'Invalid input: {str(e)}'
        }), 400

    except Exception as e:
        logger.error(f"Unexpected error in agent chat endpoint: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred. Please check your API key settings.'
        }), 500
