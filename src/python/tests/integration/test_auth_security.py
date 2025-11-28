import pytest
import os
from unittest import mock
from app import create_app

# Define a specific token for testing
TEST_TOKEN = "pytest-test-token-12345"

@pytest.fixture
def auth_app():
    """
    Create a Flask app instance specifically for auth testing.
    This is separate from the global app fixture to ensure we can
    inject the environment variable before app creation.
    """
    # Patch environment variables before creating the app
    with mock.patch.dict(os.environ, {'PYSCF_AUTH_TOKEN': TEST_TOKEN}):
        # Create app with testing config
        test_config = {
            'TESTING': True,
            # Disable other unrelated features to speed up test
            'WEBSOCKET_WATCHER_ENABLED': False, 
            'SOCKETIO': {'cors_allowed_origins': '*'}
        }
        app = create_app(server_port=5001, test_config=test_config)
        
        yield app

@pytest.fixture
def auth_client(auth_app):
    """Create a test client for the auth_app."""
    return auth_app.test_client()

def test_request_without_token(auth_client):
    """Test that requests without a token are rejected."""
    response = auth_client.get('/health')
    assert response.status_code == 401
    assert response.json['success'] is False
    assert response.json['error'] == 'Unauthorized'

def test_request_with_wrong_token(auth_client):
    """Test that requests with an incorrect token are rejected."""
    response = auth_client.get('/health', headers={'X-Auth-Token': 'wrong-token'})
    assert response.status_code == 401
    assert response.json['success'] is False

def test_request_with_correct_token(auth_client):
    """Test that requests with the correct token are accepted."""
    response = auth_client.get('/health', headers={'X-Auth-Token': TEST_TOKEN})
    assert response.status_code == 200
    # Health endpoint might return different structure, but 200 is key

def test_options_request_cors(auth_client):
    """Test that OPTIONS requests are allowed without token (for CORS)."""
    response = auth_client.options('/health')
    # Flask/CORS might return 200 or 204
    assert response.status_code in [200, 204]
