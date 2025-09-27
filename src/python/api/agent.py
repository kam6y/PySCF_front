"""
AI Agent API endpoints.
Handles chat interactions with the AI agent for molecular analysis and assistance.
"""

import logging
import json
from flask import Blueprint, jsonify, Response, stream_with_context
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
    """Chat with AI agent for molecular analysis and assistance using Server-Sent Events."""
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
            
            def fallback_stream():
                fallback_message = "The AI agent service is currently not available. Please check the server configuration."
                sse_data = json.dumps({"type": "chunk", "payload": {"text": fallback_message}})
                yield f"data: {sse_data}\n\n"
                done_data = json.dumps({"type": "done"})
                yield f"data: {done_data}\n\n"
            
            return Response(stream_with_context(fallback_stream()), content_type='text/event-stream')

        def stream():
            """エージェントからの応答をSSE形式でストリーミングするジェネレータ"""
            try:
                # ジェネレータ版のchatメソッドを呼び出す
                chunks_iterator = molecular_agent.chat(body.message, body.history or [])
                for chunk in chunks_iterator:
                    # SSE形式でデータをフォーマットしてyield
                    sse_data = json.dumps({"type": "chunk", "payload": {"text": chunk}})
                    yield f"data: {sse_data}\n\n"
            except Exception as e:
                logger.error(f"Error during streaming: {e}", exc_info=True)
                error_data = json.dumps({"type": "error", "payload": {"message": "An error occurred during the stream."}})
                yield f"data: {error_data}\n\n"
            finally:
                # ストリームの終了を通知するイベント
                done_data = json.dumps({"type": "done"})
                yield f"data: {done_data}\n\n"

        # Responseオブジェクトをストリームとして返す
        return Response(stream_with_context(stream()), content_type='text/event-stream')

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