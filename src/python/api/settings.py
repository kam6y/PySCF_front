"""
Settings management API endpoints.
Handles application settings retrieval and updates.
"""


import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from quantum_calc import get_process_manager, get_current_settings, update_app_settings
from quantum_calc.resource_manager import get_resource_manager
from generated_models import SettingsUpdateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create settings blueprint
settings_bp = Blueprint('settings', __name__)


@settings_bp.route('/api/settings', methods=['GET'])
def get_settings():
    """Get current application settings."""
    try:
        logger.info("Getting application settings")
        
        # Get current settings
        settings = get_current_settings()
        
        logger.info(f"Successfully retrieved settings: {settings}")
        return jsonify({
            'success': True,
            'data': {
                'settings': settings.model_dump()
            }
        })
    except Exception as e:
        logger.error(f"Failed to retrieve settings: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve settings: {str(e)}'
        }), 500


@settings_bp.route('/api/settings', methods=['PUT'])
@validate()
def update_settings(body: SettingsUpdateRequest):
    """Update application settings."""
    try:
        logger.info(f"Updating application settings: {body}")
        
        # Extract settings from root model
        new_settings = body.root if hasattr(body, 'root') else body
        
        # Update settings
        updated_settings = update_app_settings(new_settings.model_dump())
        
        # Update process manager with new parallel instance limit
        process_manager = get_process_manager()
        process_manager.set_max_parallel_instances(updated_settings.max_parallel_instances)
        
        # Update resource manager with new resource constraints
        resource_manager = get_resource_manager()
        resource_manager.update_resource_constraints(
            max_cpu_utilization_percent=updated_settings.max_cpu_utilization_percent,
            max_memory_utilization_percent=updated_settings.max_memory_utilization_percent
        )
        
        logger.info(f"Successfully updated settings: {updated_settings}")
        return jsonify({
            'success': True,
            'data': {
                'settings': updated_settings.model_dump()
            }
        })
    except Exception as e:
        logger.error(f"Failed to update settings: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to update settings: {str(e)}'
        }), 400 if isinstance(e, ValueError) else 500