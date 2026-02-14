"""
Settings management API endpoints.
Handles application settings retrieval and updates.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from services import get_settings_service
from generated_models import SettingsUpdateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create settings blueprint
settings_bp = Blueprint('settings', __name__)


@settings_bp.route('/api/settings', methods=['GET'])
def get_settings():
    """Get current application settings."""
    settings_service = get_settings_service()

    # Call service layer
    settings = settings_service.get_settings()

    return jsonify({
        'success': True,
        'data': {
            'settings': settings
        }
    })


@settings_bp.route('/api/settings', methods=['PUT'])
@validate()
def update_settings(body: SettingsUpdateRequest):
    """Update application settings."""
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
