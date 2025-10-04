"""
Settings management API endpoints.
Handles application settings retrieval and updates.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from services import get_settings_service, ServiceError
from generated_models import SettingsUpdateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create settings blueprint
settings_bp = Blueprint('settings', __name__)


@settings_bp.route('/api/settings', methods=['GET'])
def get_settings():
    """Get current application settings."""
    try:
        settings_service = get_settings_service()
        
        # Call service layer
        settings = settings_service.get_settings()
        
        return jsonify({
            'success': True,
            'data': {
                'settings': settings
            }
        })
    except ServiceError as e:
        logger.error(f"Service error retrieving settings: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Failed to retrieve settings: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500


@settings_bp.route('/api/settings', methods=['PUT'])
@validate()
def update_settings(body: SettingsUpdateRequest):
    """Update application settings."""
    try:
        settings_service = get_settings_service()
        
        # Extract settings from root model
        new_settings = body.root if hasattr(body, 'root') else body
        
        # Call service layer
        updated_settings = settings_service.update_settings(new_settings.model_dump())
        
        return jsonify({
            'success': True,
            'data': {
                'settings': updated_settings
            }
        })
    except ServiceError as e:
        logger.error(f"Service error updating settings: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Failed to update settings: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500