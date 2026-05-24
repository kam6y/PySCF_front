from fastapi.testclient import TestClient

from app import create_fastapi_app


def make_auth_client(monkeypatch, token='test-token', env='development'):
    if token is None:
        monkeypatch.delenv('PYSCF_AUTH_TOKEN', raising=False)
    else:
        monkeypatch.setenv('PYSCF_AUTH_TOKEN', token)
    monkeypatch.setenv('PYSCF_ENV', env)
    app = create_fastapi_app(server_port=5000, test_config={'TESTING': env != 'production'})
    return TestClient(app)


def test_missing_token_in_production(monkeypatch):
    """Test that requests without a token are rejected in production."""
    with make_auth_client(monkeypatch, token=None, env='production') as client:
        response = client.get('/health')

    assert response.status_code == 401
    assert response.json() == {
        'success': False,
        'error': 'Unauthorized: Missing authentication token',
    }


def test_missing_token_in_development(monkeypatch):
    """Test that requests without a token are allowed in development."""
    with make_auth_client(monkeypatch, token=None, env='development') as client:
        response = client.get('/health')

    assert response.status_code != 401
