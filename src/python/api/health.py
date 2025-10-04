"""
Health check API endpoints.
Provides basic health monitoring for the Flask application.
"""

from flask import Blueprint, jsonify, current_app

# Create health blueprint
health_bp = Blueprint('health', __name__)


@health_bp.route('/health', methods=['GET'])
def health_check():
    """
    Health check endpoint.

    Returns application status and version information from Flask app.config,
    which serves as the single source of truth for configuration.
    """
    try:
        # Get version from Flask app.config (single source of truth)
        version = current_app.config.get('APP_VERSION', 'unknown')
    except RuntimeError:
        # Outside Flask context (should not happen in normal operation)
        version = 'unknown'

    return jsonify({
        'status': 'ok',
        'service': 'pyscf-front-api',
        'version': version
    })