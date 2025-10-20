"""
AI Agent API endpoints.
Handles chat interactions with the multi-agent supervisor system for molecular analysis and research.
"""

import logging
import json
from typing import Iterator, Dict, Any, List
from flask import Blueprint, jsonify, Response, stream_with_context
from flask_pydantic import validate
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage

from generated_models import AgentChatRequest, ExecuteConfirmedActionRequest, AgentActionType
from agent.quantum_calculator import tools
from agent.graph import get_compiled_graph
from services.chat_history_service import get_chat_history_service

# Set up logging
logger = logging.getLogger(__name__)

# Constants
MAX_MESSAGE_LENGTH = 100000  # Maximum allowed message length in characters

# Global cache for compiled graph (lazy-loaded on first request)
# This prevents re-compilation on every request and ensures initialization
# happens only once, avoiding potential file system side effects on subsequent requests
_compiled_graph_cache = None


def _get_or_create_compiled_graph():
    """
    Get the compiled graph, creating it if necessary.

    This function implements lazy loading with caching to ensure:
    1. Graph is compiled only once during the application lifecycle
    2. First-time initialization side effects (file creation, module loading) occur
       on the first request rather than at import time
    3. Subsequent requests use the cached graph for better performance

    Returns:
        Compiled LangGraph supervisor application
    """
    global _compiled_graph_cache

    if _compiled_graph_cache is None:
        logger.info("Compiling supervisor graph (first-time initialization)")
        _compiled_graph_cache = get_compiled_graph()
        logger.info("Supervisor graph compiled and cached successfully")

    return _compiled_graph_cache


def _is_internal_transfer_message(content: str) -> bool:
    """
    Check if a message is an internal transfer message that should be filtered out.
    
    Args:
        content: Message content to check
        
    Returns:
        True if the message is an internal transfer message, False otherwise
    """
    internal_patterns = [
        "Successfully transferred to",
        "Successfully transferred back to",
    ]
    return any(pattern in content for pattern in internal_patterns)


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


def _normalize_markdown_codeblocks(text: str) -> str:
    """
    Normalize markdown to ensure code blocks are properly formatted.

    This function ensures that markdown code blocks (```) have proper spacing
    before them, which is required for correct markdown rendering.

    Specifically:
    - Ensures there is a newline before ``` code blocks
    - Handles cases like "text,```code" -> "text,\n```code"
    - Handles cases where content might be a list instead of a string

    Args:
        text: Raw markdown text that may have formatting issues

    Returns:
        Normalized markdown text with proper code block spacing

    Example:
        >>> text = "HOMO energy is -7.89 eV.,```orbital-viewer\ncalc_id: foo"
        >>> _normalize_markdown_codeblocks(text)
        "HOMO energy is -7.89 eV.,\n```orbital-viewer\ncalc_id: foo"
    """
    import re

    # Handle case where text might not be a string
    if not isinstance(text, str):
        # If it's a list, join the text content
        if isinstance(text, list):
            text_parts = []
            for item in text:
                if isinstance(item, str):
                    text_parts.append(item)
                elif isinstance(item, dict) and 'text' in item:
                    text_parts.append(item['text'])
                elif hasattr(item, 'text'):
                    text_parts.append(item.text)
            text = ''.join(text_parts)
        else:
            # Convert to string as fallback
            text = str(text)

    # Pattern: any non-whitespace character followed immediately by ```
    # Replace with the same character, newline, then ```
    # This handles cases like: "text,```code" or "text.```code"
    text = re.sub(r'([^\s\n])(```)', r'\1\n\2', text)

    return text


def _convert_dict_to_langchain_format(history: List[Dict[str, Any]]) -> List[BaseMessage]:
    """
    Convert frontend message history format to LangChain BaseMessage format.

    Frontend format:
        [
            {"role": "user", "parts": [{"text": "Hello"}]},
            {"role": "model", "parts": [{"text": "Hi there!"}]}
        ]

    LangChain format:
        [
            HumanMessage(content="Hello"),
            AIMessage(content="Hi there!")
        ]

    Args:
        history: List of message dictionaries in frontend format

    Returns:
        List of LangChain BaseMessage objects
    """
    converted = []

    for msg in history:
        # Handle both dict and Pydantic model objects
        if isinstance(msg, dict):
            role = msg.get("role", "")
            parts = msg.get("parts", [])
        else:
            # Pydantic model (HistoryItem)
            role_attr = getattr(msg, "role", "")
            # Handle enum values (convert to string)
            role = str(role_attr).split('.')[-1] if hasattr(role_attr, 'value') else str(role_attr)
            parts = getattr(msg, "parts", [])

        # Extract text content from parts
        text_content = ""
        for part in parts:
            if isinstance(part, dict) and "text" in part:
                text_content += part["text"]
            elif hasattr(part, "text"):
                # Pydantic model
                text_content += part.text

        # Convert to appropriate message type
        if role == "user":
            converted.append(HumanMessage(content=text_content))
        elif role == "model":
            converted.append(AIMessage(content=text_content))
        else:
            logger.warning(f"Unknown role '{role}' in message history, skipping")

    return converted


def _create_supervisor_stream(message: str, history: list, session_id: str = None) -> Iterator[str]:
    """
    Create an SSE stream for Supervisor multi-agent responses with worker status tracking.

    This function integrates the LangGraph Supervisor to coordinate specialized worker agents
    and provides real-time status updates about which agent is currently executing.

    New SSE event types:
    - agent_status: Worker execution status updates
      {"type": "agent_status", "payload": {"status": "running"|"completed"|"responding", "agent": "agent_name"}}
    - chunk: Text chunks from AI responses
    - done: Stream completion
    - error: Error messages

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
    # Track whether DB save was successful (for fallback logic)
    db_save_successful = False

    try:
        logger.debug(f"Starting LangGraph Supervisor stream for message: {message[:100]}{'...' if len(message) > 100 else ''}")

        # Get compiled supervisor graph (cached after first call)
        graph_app = _get_or_create_compiled_graph()

        # Convert history to LangChain format
        langchain_history = _convert_dict_to_langchain_format(history)
        logger.debug(f"Converted {len(history)} history messages to LangChain format")

        # Add current user message
        all_messages = langchain_history + [HumanMessage(content=message)]

        # Prepare graph input
        graph_input = {"messages": all_messages}

        # Track current worker agent
        current_worker = None
        supervisor_response_started = False

        # Stream graph execution with updates mode to track node executions
        logger.debug("Starting LangGraph Supervisor stream with updates mode and subgraphs enabled")
        for chunk in graph_app.stream(
            graph_input,
            stream_mode="updates",
            subgraphs=True  # Enable subgraph (worker) updates
        ):
            logger.debug(f"Received chunk: {type(chunk)}, is_tuple: {isinstance(chunk, tuple)}")

            # Check if this is a subgraph update (worker agent execution)
            if isinstance(chunk, tuple):
                namespace, node_updates = chunk

                # Subgraph execution detected
                if len(namespace) > 0:
                    # Extract worker name from namespace
                    # Format: ('quantum_calculation_worker:uuid',) or ('research_expert:uuid',)
                    worker_id = namespace[-1]
                    worker_name = worker_id.split(":")[0] if ":" in worker_id else worker_id

                    # Worker execution started
                    if current_worker != worker_name:
                        logger.info(f"Worker agent '{worker_name}' started execution")
                        current_worker = worker_name
                        yield _format_sse_event("agent_status", {
                            "status": "running",
                            "agent": worker_name
                        })

                    # Note: We'll send "completed" when we detect the parent graph update after worker finishes
                else:
                    # Parent graph update (Supervisor)
                    logger.debug(f"Parent graph update: {list(node_updates.keys())}")

                    # If we were tracking a worker, mark it as completed
                    if current_worker:
                        logger.info(f"Worker agent '{current_worker}' completed execution")
                        yield _format_sse_event("agent_status", {
                            "status": "completed",
                            "agent": current_worker
                        })
                        current_worker = None

                    # Process supervisor's response
                    for node_name, data in node_updates.items():
                        if "messages" in data:
                            messages = data["messages"]
                            if messages:
                                # Get the last message
                                last_message = messages[-1]
                                if hasattr(last_message, "content") and last_message.content:
                                    # Filter out internal transfer messages
                                    if _is_internal_transfer_message(last_message.content):
                                        logger.debug(f"Filtering internal transfer message: {last_message.content[:100]}")
                                        continue
                                    
                                    # Notify that supervisor is responding
                                    if not supervisor_response_started:
                                        logger.info("Supervisor started generating response")
                                        yield _format_sse_event("agent_status", {
                                            "status": "responding",
                                            "agent": "supervisor"
                                        })
                                        supervisor_response_started = True

                                    # Stream the supervisor's response (with markdown normalization)
                                    normalized_content = _normalize_markdown_codeblocks(last_message.content)
                                    logger.info(f"Streaming supervisor response, length: {len(normalized_content)}, preview: {normalized_content[:200]}...")

                                    # Accumulate for saving
                                    accumulated_response.append(normalized_content)

                                    yield _format_sse_event("chunk", {"text": normalized_content})
            else:
                # Regular dict format (parent graph without namespace)
                logger.debug(f"Dict format update: {list(chunk.keys())}")

                # If we were tracking a worker, mark it as completed
                if current_worker:
                    logger.info(f"Worker agent '{current_worker}' completed execution")
                    yield _format_sse_event("agent_status", {
                        "status": "completed",
                        "agent": current_worker
                    })
                    current_worker = None

                # Process updates
                for node_name, data in chunk.items():
                    if "messages" in data:
                        messages = data["messages"]
                        if messages:
                            last_message = messages[-1]
                            if hasattr(last_message, "content") and last_message.content:
                                # Filter out internal transfer messages
                                if _is_internal_transfer_message(last_message.content):
                                    logger.debug(f"Filtering internal transfer message: {last_message.content[:100]}")
                                    continue
                                
                                # Notify that supervisor is responding
                                if not supervisor_response_started:
                                    logger.info("Supervisor started generating response")
                                    yield _format_sse_event("agent_status", {
                                        "status": "responding",
                                        "agent": "supervisor"
                                    })
                                    supervisor_response_started = True

                                # Stream the response (with markdown normalization)
                                normalized_content = _normalize_markdown_codeblocks(last_message.content)
                                logger.info(f"Streaming response from {node_name}, length: {len(normalized_content)}, preview: {normalized_content[:200]}...")

                                # Accumulate for saving
                                accumulated_response.append(normalized_content)

                                yield _format_sse_event("chunk", {"text": normalized_content})

        logger.debug("Stream completed successfully")

        # Save AI response to database BEFORE sending completion event
        # This ensures database write completes before frontend invalidates cache
        if session_id and accumulated_response:
            try:
                complete_response = ''.join(accumulated_response)
                chat_service = get_chat_history_service()
                chat_service.add_message(session_id, "model", complete_response)
                db_save_successful = True
                logger.info(f"Saved AI response to session: {session_id} (length: {len(complete_response)} chars)")
            except Exception as e:
                logger.error(f"Failed to save AI response to session {session_id}: {e}", exc_info=True)
                # Don't fail the stream due to DB error, continue with completion event

        # Send completion event AFTER database write completes
        # This ensures 'done' is sent when the stream completes normally and DB is updated
        logger.debug("Sending stream completion event")
        yield _format_sse_event("done")

    except GeneratorExit:
        # Client disconnected - clean up gracefully
        logger.info("Client disconnected from SSE stream (GeneratorExit)")
        # Don't yield anything here - the connection is already closed
        # Yielding after GeneratorExit causes RuntimeError
        raise  # Re-raise to properly close the generator

    except Exception as e:
        logger.error(f"Error during LangGraph Supervisor streaming: {e}", exc_info=True)
        # Add error message to accumulated response
        error_msg = f"\n\n[Error: {str(e)}]"
        accumulated_response.append(error_msg)

        try:
            yield _format_sse_event("error", {"message": f"An error occurred during the stream: {str(e)}"})
        except (BrokenPipeError, ConnectionResetError, GeneratorExit):
            # Connection already closed, can't send error message
            logger.debug("Unable to send error message - connection closed")

    finally:
        # Fallback: Save AI response to database if not already saved
        # This runs only if an error occurred before the normal save point
        if session_id and accumulated_response and not db_save_successful:
            try:
                complete_response = ''.join(accumulated_response)
                chat_service = get_chat_history_service()
                chat_service.add_message(session_id, "model", complete_response)
                logger.warning(f"Fallback save: AI response saved to session {session_id} after error (length: {len(complete_response)} chars)")
            except Exception as e:
                logger.error(f"Fallback save failed for session {session_id}: {e}", exc_info=True)


# Create agent blueprint
agent_bp = Blueprint('agent', __name__)


@agent_bp.route('/api/agent/chat', methods=['POST'])
@validate()
def chat_with_agent(body: AgentChatRequest):
    """Chat with AI agent using multi-agent supervisor architecture with Server-Sent Events."""
    try:
        # Input validation
        if not body.message or not body.message.strip():
            raise ValueError("Message cannot be empty")

        if len(body.message) > MAX_MESSAGE_LENGTH:
            raise ValueError(f"Message is too long (maximum {MAX_MESSAGE_LENGTH} characters)")

        session_id = body.session_id if hasattr(body, 'session_id') else None
        logger.info(f"Processing agent chat request - Message length: {len(body.message)}, History entries: {len(body.history or [])}, Session ID: {session_id}")
        logger.debug(f"Message preview: {body.message[:100]}{'...' if len(body.message) > 100 else ''}")

        # Use LangGraph Supervisor multi-agent system
        return Response(
            stream_with_context(_create_supervisor_stream(body.message, body.history or [], session_id)),
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
            'error': 'An internal server error occurred in the agent. Please contact support if this persists.'
        }), 500


@agent_bp.route('/api/agent/execute-confirmed-action', methods=['POST'])
@validate()
def execute_confirmed_action(body: ExecuteConfirmedActionRequest):
    """Execute a confirmed destructive action requested by AI agent."""
    try:
        logger.info(f"Processing confirmed action execution - Action type: {body.action_type}")

        # Handle different action types
        if body.action_type == AgentActionType.delete_calculation:
            # Validate that calculation_id is provided
            if not body.calculation_id:
                raise ValueError("calculation_id is required for delete_calculation action")

            logger.info(f"Executing confirmed deletion for calculation: {body.calculation_id}")

            # Execute the confirmed deletion using the internal function
            result = tools._execute_confirmed_deletion(body.calculation_id)

            if result['success']:
                logger.info(f"Deletion successful: {result['message']}")
                return jsonify({
                    'success': True,
                    'data': {
                        'message': result['message'],
                        'action_type': 'delete_calculation',
                        'calculation_id': result['calculation_id']
                    }
                }), 200
            else:
                logger.warning(f"Deletion failed: {result['message']}")
                return jsonify({
                    'success': False,
                    'error': result['message']
                }), 400

        else:
            # Unknown action type
            logger.warning(f"Unknown action type requested: {body.action_type}")
            raise ValueError(f"Unknown action type: {body.action_type}")

    except ValueError as e:
        logger.warning(f"Validation error in execute confirmed action: {e}")
        return jsonify({
            'success': False,
            'error': f'Invalid input: {str(e)}'
        }), 400

    except Exception as e:
        logger.error(f"Unexpected error in execute confirmed action endpoint: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred while executing the action. Please contact support if this persists.'
        }), 500
