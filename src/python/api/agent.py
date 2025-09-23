"""
AI Agent API endpoints.
Handles chat interactions with the AI agent for molecular analysis and assistance.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from generated_models import AgentChatRequest, AgentChatResponse
from agent.molecular_agent import MolecularAgent

# Set up logging
logger = logging.getLogger(__name__)

# Create agent blueprint
agent_bp = Blueprint('agent', __name__)

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
    """Chat with AI agent for molecular analysis and assistance."""
    try:
        # Input validation
        if not body.message or not body.message.strip():
            raise ValueError("Message cannot be empty")
        
        if len(body.message) > 10000:  # Set reasonable message length limit
            raise ValueError("Message is too long (maximum 10,000 characters)")
        
        logger.info(f"Processing agent chat request - Message length: {len(body.message)}, History entries: {len(body.history or [])}")
        logger.debug(f"Message preview: {body.message[:100]}{'...' if len(body.message) > 100 else ''}")
        
        # Check if molecular agent is available
        if not molecular_agent:
            logger.warning("MolecularAgent not available, returning fallback response")
            reply = "The AI agent service is currently not available. Please check the server configuration."
            response = AgentChatResponse(
                success=True,
                data={"reply": reply}
            )
            return jsonify(response.model_dump()), 200
        
        # Use the molecular agent to generate response
        reply = molecular_agent.chat(body.message, body.history or [])
        
        logger.info(f"Agent chat completed successfully - Response length: {len(reply)}")
        logger.debug(f"Response preview: {reply[:100]}{'...' if len(reply) > 100 else ''}")
        
        response = AgentChatResponse(
            success=True,
            data={"reply": reply}
        )
        return jsonify(response.model_dump()), 200

    except ValueError as e:
        logger.warning(f"Validation error in agent chat: {e}")
        return jsonify({
            'success': False, 
            'error': f'Invalid input: {str(e)}'
        }), 400
    
    except ConnectionError as e:
        logger.error(f"Connection error in agent chat: {e}")
        return jsonify({
            'success': False, 
            'error': 'Unable to connect to AI service. Please try again later.'
        }), 503
    
    except Exception as e:
        logger.error(f"Unexpected error in agent chat endpoint: {e}", exc_info=True)
        return jsonify({
            'success': False, 
            'error': 'An internal server error occurred in the agent. Please contact support if this persists.'
        }), 500