import importlib
import os

from fastapi import FastAPI


def test_thread_control_vars_are_set_before_app_import(monkeypatch):
    thread_vars = [
        'OMP_NUM_THREADS',
        'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS',
        'BLIS_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS',
        'NUMEXPR_NUM_THREADS',
    ]
    for name in thread_vars:
        monkeypatch.delenv(name, raising=False)

    module = importlib.reload(importlib.import_module('app'))

    for name in thread_vars:
        assert os.environ[name] == '1'
    assert hasattr(module, 'fastapi_app')
    assert hasattr(module, 'sio')
    assert hasattr(module, 'app')


def test_create_fastapi_app_disables_default_docs():
    from app import create_fastapi_app

    app = create_fastapi_app(server_port=5000, test_config={'TESTING': True})

    assert isinstance(app, FastAPI)
    routes = {route.path for route in app.routes}
    assert '/docs' not in routes
    assert '/redoc' not in routes
    assert '/openapi.json' not in routes
    assert '/api-docs' in routes
    assert '/api-docs/spec.json' in routes


def test_config_is_stored_on_app_state():
    from app import create_fastapi_app

    app = create_fastapi_app(
        server_port=5123,
        test_config={
            'TESTING': True,
            'CALCULATIONS_DIR': '/tmp/pyscf-fastapi-test',
            'WEBSOCKET_WATCHER_ENABLED': False,
        },
    )

    assert app.state.TESTING is True
    assert app.state.SERVER_PORT == 5123
    assert app.state.CALCULATIONS_DIR == '/tmp/pyscf-fastapi-test'
    assert app.state.WEBSOCKET_WATCHER_ENABLED is False
