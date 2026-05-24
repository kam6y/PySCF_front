from fastapi import Query
from fastapi.testclient import TestClient
from pydantic import BaseModel


class SampleBody(BaseModel):
    required_name: str


def test_validation_error_returns_400_envelope(app):
    @app.post('/_test/validation-body')
    def validation_body(body: SampleBody):
        return {'ok': True}

    with TestClient(app) as client:
        response = client.post('/_test/validation-body', json={})

    assert response.status_code == 400
    assert response.json()['success'] is False
    assert response.json()['error'].startswith('Validation failed:')


def test_malformed_json_returns_400_envelope(app):
    @app.post('/_test/validation-json')
    def validation_json(body: SampleBody):
        return {'ok': True}

    with TestClient(app) as client:
        response = client.post(
            '/_test/validation-json',
            content=b'{"required_name":',
            headers={'Content-Type': 'application/json'},
        )

    assert response.status_code == 400
    assert response.json()['success'] is False
    assert response.json()['error'].startswith('Validation failed:')


def test_invalid_query_parameter_returns_400_envelope(app):
    @app.get('/_test/validation-query')
    def validation_query(limit: int = Query(10)):
        return {'limit': limit}

    with TestClient(app) as client:
        response = client.get('/_test/validation-query?limit=not-an-int')

    assert response.status_code == 400
    assert response.json()['success'] is False
    assert response.json()['error'].startswith('Validation failed:')


def test_options_preflight_succeeds_without_auth(monkeypatch, app):
    monkeypatch.setenv('PYSCF_AUTH_TOKEN', 'secret-token')

    with TestClient(app) as client:
        response = client.options(
            '/health',
            headers={
                'Origin': 'http://127.0.0.1:3000',
                'Access-Control-Request-Method': 'PATCH',
                'Access-Control-Request-Headers': 'X-Auth-Token, Content-Type',
            },
        )

    assert response.status_code in {200, 204}
    assert 'X-Auth-Token' in response.headers['access-control-allow-headers']


def test_auth_error_includes_cors_headers(monkeypatch, app):
    monkeypatch.setenv('PYSCF_AUTH_TOKEN', 'secret-token')

    with TestClient(app) as client:
        response = client.get(
            '/health',
            headers={
                'Origin': 'http://127.0.0.1:3000',
                'X-Auth-Token': 'wrong-token',
            },
        )

    assert response.status_code == 401
    assert response.json()['success'] is False
    assert response.headers['access-control-allow-origin'] == 'http://127.0.0.1:3000'


def test_unhandled_exception_returns_500_envelope(app):
    @app.get('/_test/unhandled-error')
    def unhandled_error():
        raise RuntimeError('boom')

    with TestClient(app, raise_server_exceptions=False) as client:
        response = client.get('/_test/unhandled-error')

    assert response.status_code == 500
    assert response.json() == {
        'success': False,
        'error': 'An internal server error occurred.',
    }


def test_quantum_calculate_requires_calculation_method_before_model_validation(client):
    response = client.post('/api/quantum/calculate', json={'molecule': {'atoms': []}})

    assert response.status_code == 400
    body = response.json()
    assert body['success'] is False
    assert 'calculation_method' in body['error']


def test_quantum_calculate_rejects_inapplicable_method_parameter(client):
    response = client.post(
        '/api/quantum/calculate',
        json={
            'calculation_method': 'HF',
            'xyz': 'H 0 0 0',
            'basis_function': 'sto-3g',
            'exchange_correlation': 'b3lyp',
        },
    )

    assert response.status_code == 400
    body = response.json()
    assert body['success'] is False
    assert 'exchange_correlation' in body['error']


def test_gpu4pyscf_install_rejects_non_loopback_client(client, monkeypatch):
    import api.system as system_api

    monkeypatch.setattr(system_api, '_get_client_host', lambda request: '203.0.113.10')

    response = client.post('/api/system/gpu4pyscf-install')

    assert response.status_code == 403
    assert response.json() == {
        'success': False,
        'error': 'GPU4PySCF installation is only available from the local machine.',
    }
