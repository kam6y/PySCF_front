"""
AI Agent API endpoints.
Handles chat interactions with the AI agent for molecular analysis and assistance.
"""

import logging
import json
from typing import Iterator, Dict, Any
from flask import Blueprint, jsonify, Response, stream_with_context
from flask_pydantic import validate

from generated_models import AgentChatRequest, ExecuteConfirmedActionRequest, AgentActionType
from agent.molecular_agent import MolecularAgent
from agent import tools

# Set up logging
logger = logging.getLogger(__name__)

# Constants
MAX_MESSAGE_LENGTH = 10000  # Maximum allowed message length in characters


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


# Create agent blueprint
agent_bp = Blueprint('agent', __name__)

def _create_fallback_stream() -> Iterator[str]:
    """Create a fallback SSE stream when the molecular agent is not available."""
    fallback_message = "The AI agent service is currently not available. Please check the server configuration."
    yield _format_sse_event("chunk", {"text": fallback_message})
    yield _format_sse_event("done")


def _create_agent_stream(message: str, history: list) -> Iterator[str]:
    """
    Create an SSE stream for agent responses.
    
    Args:
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


# Initialize the molecular agent (singleton instance)
try:
    molecular_agent = MolecularAgent()
    logger.info("MolecularAgent initialized successfully")
except Exception as e:
    logger.error(f"Failed to initialize MolecularAgent: {e}")
    molecular_agent = None


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
        
        # Check if molecular agent is available
        if not molecular_agent:
            logger.warning("MolecularAgent not available, returning fallback response")
            return Response(
                stream_with_context(_create_fallback_stream()),
                content_type='text/event-stream'
            )

        # Return streaming response
        return Response(
            stream_with_context(_create_agent_stream(body.message, body.history or [])),
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