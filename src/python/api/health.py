"""
Health check API endpoints.
Provides basic health monitoring for the Flask application.
"""

from flask import Blueprint, jsonify

# Create health blueprint
health_bp = Blueprint('health', __name__)


@health_bp.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    # Import here to avoid circular imports when using SERVER_CONFIG
    import os
    import sys
    import json
    
    # Load server configuration for version info
    try:
        config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'config', 'server-config.json')
        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as f:
                config = json.load(f)
            version = config.get('app_info', {}).get('version', 'unknown')
        else:
            version = 'unknown'
    except Exception:
        version = 'unknown'
    
    return jsonify({
        'status': 'ok',
        'service': 'pyscf-front-api',
        'version': version
    })