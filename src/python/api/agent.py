"""
AI Agent API endpoints.
Handles chat interactions with the AI agent for molecular analysis and assistance.
"""

import logging
import json
from typing import Iterator, Dict, Any, List
from flask import Blueprint, jsonify, Response, stream_with_context
from flask_pydantic import validate
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage

from generated_models import AgentChatRequest, ExecuteConfirmedActionRequest, AgentActionType
from agent.molecular_agent import MolecularAgent
from agent import tools
from agent.graph import get_compiled_graph

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


# Create agent blueprint
agent_bp = Blueprint('agent', __name__)

# Global agent instance (will be initialized lazily)
_molecular_agent = None
_agent_lock = __import__('threading').Lock()


def _create_fallback_stream() -> Iterator[str]:
    """Create a fallback SSE stream when the molecular agent is not available."""
    fallback_message = "The AI agent service is currently not available. Please check the server configuration."
    yield _format_sse_event("chunk", {"text": fallback_message})
    yield _format_sse_event("done")


def _create_agent_stream(molecular_agent: MolecularAgent, message: str, history: list) -> Iterator[str]:
    """
    Create an SSE stream for agent responses.
    
    Args:
        molecular_agent: The MolecularAgent instance to use for generating responses
        message: User's message
        history: Chat history
        
    Yields:
        SSE formatted strings
    """
    chunks_count = 0
    try:
        logger.debug(f"Starting stream for message: {message[:100]}{'...' if len(message) > 100 else ''}")
        
        # Call the agent's chat method
        chunks_iterator = molecular_agent.chat(message, history)
        for chunk in chunks_iterator:
            chunks_count += 1
            logger.debug(f"Streaming chunk #{chunks_count}, length: {len(chunk)}")
            yield _format_sse_event("chunk", {"text": chunk})
        
        logger.debug(f"Stream completed successfully with {chunks_count} chunks")
        
    except Exception as e:
        logger.error(f"Error during streaming after {chunks_count} chunks: {e}", exc_info=True)
        yield _format_sse_event("error", {"message": "An error occurred during the stream."})
    finally:
        logger.debug("Sending stream completion event")
        yield _format_sse_event("done")


def _create_graph_stream(message: str, history: list) -> Iterator[str]:
    """
    Create an SSE stream for LangGraph multi-agent responses.

    This function integrates the LangGraph dispatcher to route queries
    to the appropriate specialist agent (MolecularAgent or ResearchAgent).
    Uses real-time token streaming with stream_mode="messages".

    Args:
        message: User's message
        history: Chat history in frontend format (dict list)

    Yields:
        SSE formatted strings
    """
    chunks_count = 0
    try:
        logger.debug(f"Starting LangGraph stream for message: {message[:100]}{'...' if len(message) > 100 else ''}")

        # Get compiled graph
        graph_app = get_compiled_graph()

        # Convert history to LangChain format
        langchain_history = _convert_dict_to_langchain_format(history)
        logger.debug(f"Converted {len(history)} history messages to LangChain format")

        # Add current user message
        all_messages = langchain_history + [HumanMessage(content=message)]

        # Prepare graph input
        graph_input = {"messages": all_messages}

        # Stream graph execution with custom mode to capture tokens from nodes
        logger.debug("Starting LangGraph stream with custom mode")
        for chunk in graph_app.stream(
            graph_input,
            stream_mode="custom"
        ):
            # Stream custom tokens emitted by nodes via get_stream_writer()
            if isinstance(chunk, dict) and "token" in chunk:
                token = chunk["token"]
                chunks_count += 1
                logger.debug(f"Streaming token #{chunks_count}, length: {len(token)}")
                yield _format_sse_event("chunk", {"text": token})

        if chunks_count == 0:
            logger.warning("No tokens streamed from graph, returning fallback")
            yield _format_sse_event("chunk", {"text": "応答を生成できませんでした。もう一度お試しください。"})

        logger.debug(f"Stream completed successfully with {chunks_count} tokens")

    except Exception as e:
        logger.error(f"Error during LangGraph streaming after {chunks_count} chunks: {e}", exc_info=True)
        yield _format_sse_event("error", {"message": f"An error occurred during the stream: {str(e)}"})
    finally:
        logger.debug("Sending stream completion event")
        yield _format_sse_event("done")


def get_molecular_agent():
    """
    Get or create the MolecularAgent instance (lazy initialization).
    This approach ensures that the agent is initialized after Gunicorn fork,
    avoiding issues with Google Genai SDK and multiprocessing.
    
    Returns:
        MolecularAgent or None: The agent instance, or None if initialization fails
    """
    global _molecular_agent
    
    # Fast path: if already initialized, return immediately
    if _molecular_agent is not None:
        return _molecular_agent
    
    # Slow path: acquire lock and initialize
    with _agent_lock:
        # Double-check pattern: another thread might have initialized while we waited
        if _molecular_agent is not None:
            return _molecular_agent
        
        try:
            logger.info("Initializing MolecularAgent (lazy initialization after fork)")
            _molecular_agent = MolecularAgent()
            logger.info("MolecularAgent initialized successfully")
            return _molecular_agent
        except Exception as e:
            logger.error(f"Failed to initialize MolecularAgent: {e}", exc_info=True)
            return None


@agent_bp.route('/api/agent/chat', methods=['POST'])
@validate()
def chat_with_agent(body: AgentChatRequest):
    """Chat with AI agent for molecular analysis and assistance using Server-Sent Events."""
    try:
        # Input validation
        if not body.message or not body.message.strip():
            raise ValueError("Message cannot be empty")
        
        if len(body.message) > MAX_MESSAGE_LENGTH:
            raise ValueError(f"Message is too long (maximum {MAX_MESSAGE_LENGTH} characters)")
        
        logger.info(f"Processing agent chat request - Message length: {len(body.message)}, History entries: {len(body.history or [])}")
        logger.debug(f"Message preview: {body.message[:100]}{'...' if len(body.message) > 100 else ''}")
        
        # Use LangGraph multi-agent dispatcher for routing
        return Response(
            stream_with_context(_create_graph_stream(body.message, body.history or [])),
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
        
        # Note: This endpoint doesn't require molecular agent, it directly calls tool functions
        
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