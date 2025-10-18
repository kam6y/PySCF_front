"""
Chat History API endpoints.

Provides REST API for managing chat history sessions and messages.

Note:
    - Input validation is handled by flask_pydantic @validate decorator
    - Business logic errors are handled by service layer (ServiceError)
    - All endpoints return consistent JSON response format
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from generated_models import (
    CreateChatSessionRequest,
    UpdateChatSessionRequest
)
from services.chat_history_service import get_chat_history_service
from services.exceptions import ServiceError

# Set up logging
logger = logging.getLogger(__name__)

# Create chat history blueprint
chat_history_bp = Blueprint('chat_history', __name__)


@chat_history_bp.route('/api/chat-history/sessions', methods=['GET'])
def get_chat_sessions():
    """Get all chat sessions, ordered by most recent update."""
    try:
        logger.info("Processing request to get all chat sessions")
        service = get_chat_history_service()
        data = service.list_sessions()

        return jsonify({
            'success': True,
            'data': data
        }), 200

    except ServiceError as e:
        logger.error(f"Service error in get_chat_sessions: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

    except Exception as e:
        logger.error(f"Unexpected error in get_chat_sessions: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred. Please contact support if this persists.'
        }), 500


@chat_history_bp.route('/api/chat-history/sessions', methods=['POST'])
@validate()
def create_chat_session(body: CreateChatSessionRequest):
    """
    Create a new chat session.

    Request validation is handled by @validate decorator.
    Default value for 'name' is provided by Pydantic model.
    """
    try:
        # body.name already has default value from Pydantic model ("新しいチャット")
        logger.info(f"Processing request to create chat session: '{body.name}'")
        service = get_chat_history_service()
        session = service.create_session(name=body.name)

        return jsonify({
            'success': True,
            'data': {
                'session': session
            }
        }), 201

    except ServiceError as e:
        logger.error(f"Service error in create_chat_session: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

    except Exception as e:
        logger.error(f"Unexpected error in create_chat_session: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred. Please contact support if this persists.'
        }), 500


@chat_history_bp.route('/api/chat-history/sessions/<session_id>', methods=['GET'])
def get_chat_session(session_id: str):
    """Get a specific chat session with all messages."""
    try:
        logger.info(f"Processing request to get chat session: {session_id}")
        service = get_chat_history_service()
        session_data = service.get_session_with_messages(session_id)

        if session_data is None:
            logger.warning(f"Chat session not found: {session_id}")
            return jsonify({
                'success': False,
                'error': f'Chat session not found: {session_id}'
            }), 404

        return jsonify({
            'success': True,
            'data': session_data
        }), 200

    except ServiceError as e:
        logger.error(f"Service error in get_chat_session: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

    except Exception as e:
        logger.error(f"Unexpected error in get_chat_session: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred. Please contact support if this persists.'
        }), 500


@chat_history_bp.route('/api/chat-history/sessions/<session_id>', methods=['PATCH'])
@validate()
def update_chat_session(session_id: str, body: UpdateChatSessionRequest):
    """
    Update a chat session's metadata (name).

    Request validation is handled by @validate decorator.
    """
    try:
        logger.info(f"Processing request to update chat session: {session_id}")
        service = get_chat_history_service()
        session = service.update_session(session_id, body.name)

        if session is None:
            logger.warning(f"Chat session not found: {session_id}")
            return jsonify({
                'success': False,
                'error': f'Chat session not found: {session_id}'
            }), 404

        return jsonify({
            'success': True,
            'data': {
                'session': session
            }
        }), 200

    except ServiceError as e:
        logger.error(f"Service error in update_chat_session: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

    except Exception as e:
        logger.error(f"Unexpected error in update_chat_session: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred. Please contact support if this persists.'
        }), 500


@chat_history_bp.route('/api/chat-history/sessions/<session_id>', methods=['DELETE'])
def delete_chat_session(session_id: str):
    """Delete a chat session and all its messages."""
    try:
        logger.info(f"Processing request to delete chat session: {session_id}")
        service = get_chat_history_service()
        deleted = service.delete_session(session_id)

        if not deleted:
            logger.warning(f"Chat session not found: {session_id}")
            return jsonify({
                'success': False,
                'error': f'Chat session not found: {session_id}'
            }), 404

        return jsonify({
            'success': True,
            'data': {
                'message': 'Chat session deleted successfully',
                'deleted_id': session_id
            }
        }), 200

    except ServiceError as e:
        logger.error(f"Service error in delete_chat_session: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

    except Exception as e:
        logger.error(f"Unexpected error in delete_chat_session: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred. Please contact support if this persists.'
        }), 500
