import importlib
import os
from unittest import mock

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
    assert hasattr(module, 'app')
    assert isinstance(module.app, FastAPI)
    assert not hasattr(module, 'sio')


def test_app_import_does_not_initialize_process_manager(monkeypatch):
    import app as app_module
    import quantum_calc

    initialize_mock = mock.Mock()
    monkeypatch.setattr(
        quantum_calc,
        'initialize_process_manager_with_callback',
        initialize_mock,
    )

    importlib.reload(app_module)

    initialize_mock.assert_not_called()


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
