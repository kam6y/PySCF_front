import pytest
import os
from unittest import mock
from app import create_app

@pytest.fixture
def production_app():
    """
    Create a Flask app instance mimicking production environment.
    """
    # Patch environment variables to simulate production
    # We explicitly do NOT set PYSCF_AUTH_TOKEN here to test the failure case
    with mock.patch.dict(os.environ, {'PYSCF_ENV': 'production'}, clear=True):
        # Create app with testing config
        test_config = {
            'TESTING': True,
            'WEBSOCKET_WATCHER_ENABLED': False,
            'SOCKETIO': {'cors_allowed_origins': '*'}
        }
        app = create_app(server_port=5002, test_config=test_config)
        yield app

@pytest.fixture
def production_client(production_app):
    """Create a test client for the production_app."""
    return production_app.test_client()

def test_missing_token_in_production(production_client):
    """Test that requests without a token are rejected in production."""
    # In production, missing token should result in 401
    response = production_client.get('/health')
    assert response.status_code == 401
    assert response.json['success'] is False
    assert 'Missing authentication token' in response.json['error']

def test_missing_token_in_development():
    """Test that requests without a token are allowed (with warning) in development."""
    with mock.patch.dict(os.environ, {'PYSCF_ENV': 'development'}, clear=True):
        test_config = {
            'TESTING': True,
            'WEBSOCKET_WATCHER_ENABLED': False
        }
        # We need to recreate app to pick up the env var change if it's checked at creation time
        # But verify_auth_token checks env var at request time, so it should be fine.
        # However, create_app reads config at startup.
        
        app = create_app(server_port=5003, test_config=test_config)
        client = app.test_client()
        
        # In development, missing token should be allowed (200 OK for health check)
        # Note: The health endpoint might not exist or return 200, checking 401 is the key.
        # If /health doesn't exist it returns 404, which is not 401.
        response = client.get('/health')
        assert response.status_code != 401
