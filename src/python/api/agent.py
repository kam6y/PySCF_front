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


def _create_supervisor_stream(message: str, history: list) -> Iterator[str]:
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

    Yields:
        SSE formatted strings
    """
    try:
        logger.debug(f"Starting LangGraph Supervisor stream for message: {message[:100]}{'...' if len(message) > 100 else ''}")

        # Get compiled supervisor graph
        graph_app = get_compiled_graph()

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
                                    # Notify that supervisor is responding
                                    if not supervisor_response_started:
                                        logger.info("Supervisor started generating response")
                                        yield _format_sse_event("agent_status", {
                                            "status": "responding",
                                            "agent": "supervisor"
                                        })
                                        supervisor_response_started = True

                                    # Stream the supervisor's response
                                    logger.debug(f"Streaming supervisor response, length: {len(last_message.content)}")
                                    yield _format_sse_event("chunk", {"text": last_message.content})
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
                                # Notify that supervisor is responding
                                if not supervisor_response_started:
                                    logger.info("Supervisor started generating response")
                                    yield _format_sse_event("agent_status", {
                                        "status": "responding",
                                        "agent": "supervisor"
                                    })
                                    supervisor_response_started = True

                                # Stream the response
                                logger.debug(f"Streaming response from {node_name}, length: {len(last_message.content)}")
                                yield _format_sse_event("chunk", {"text": last_message.content})

        logger.debug("Stream completed successfully")

    except Exception as e:
        logger.error(f"Error during LangGraph Supervisor streaming: {e}", exc_info=True)
        yield _format_sse_event("error", {"message": f"An error occurred during the stream: {str(e)}"})
    finally:
        logger.debug("Sending stream completion event")
        yield _format_sse_event("done")


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

        logger.info(f"Processing agent chat request - Message length: {len(body.message)}, History entries: {len(body.history or [])}")
        logger.debug(f"Message preview: {body.message[:100]}{'...' if len(body.message) > 100 else ''}")

        # Use LangGraph Supervisor multi-agent system
        return Response(
            stream_with_context(_create_supervisor_stream(body.message, body.history or [])),
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
