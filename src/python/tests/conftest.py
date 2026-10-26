"""
Central test configuration and fixtures for PySCF Front backend tests.

This module contains reusable pytest fixtures that are automatically available
to all test modules.
"""

import os
import logging
import shutil
import socket
import tempfile
import threading
import time
from collections.abc import Generator
from concurrent.futures import Executor, Future
from unittest import mock

import pytest
import uvicorn
from fastapi.testclient import TestClient

from app import create_fastapi_app

logger = logging.getLogger(__name__)


def _get_free_port() -> int:
    """Return an available localhost TCP port for a short-lived test server."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(('127.0.0.1', 0))
        return int(sock.getsockname()[1])


def _wait_for_server_start(
    server: uvicorn.Server,
    thread: threading.Thread,
    timeout: float = 5.0,
) -> None:
    """Wait until Uvicorn reports startup or fail with a clear error."""
    deadline = time.monotonic() + timeout
    while not server.started:
        if not thread.is_alive():
            raise RuntimeError('Uvicorn server thread exited before startup.')
        if time.monotonic() >= deadline:
            raise TimeoutError('Timed out waiting for Uvicorn server startup.')
        time.sleep(0.05)


def _clear_calculation_update_stream_hub() -> None:
    """Reset the SSE hub singleton between tests."""
    from services.calculation_update_stream import get_calculation_update_stream_hub

    hub = get_calculation_update_stream_hub()
    hub._global_subscribers.clear()
    hub._calculation_subscribers.clear()


# ============================================================================
# Core Application Fixtures
# ============================================================================


class DummyExecutor(Executor):
    """
    A synchronous executor that mimics the ProcessPoolExecutor interface.

    This executor runs tasks immediately and synchronously in the same thread,
    which is essential for testing asynchronous workflows without the complexity
    of actual multiprocessing. It allows tests to verify the complete workflow
    from task submission to completion in a predictable, deterministic manner.

    Usage:
        Use mocker.patch to replace ProcessPoolExecutor with DummyExecutor:
        ```
        mocker.patch(
            'module.ProcessPoolExecutor',
            new=DummyExecutor
        )
        ```
    """

    def __init__(self, *args, **kwargs):
        """Initialize the dummy executor."""
        super().__init__()
        self._shutdown = False

    def submit(self, fn, *args, **kwargs):
        """
        Execute the function immediately and return a Future with the result.

        Args:
            fn: The function to execute.
            *args: Positional arguments for the function.
            **kwargs: Keyword arguments for the function.

        Returns:
            Future: A Future object with the result or exception.
        """
        if self._shutdown:
            raise RuntimeError('Executor has been shutdown.')

        future = Future()
        try:
            result = fn(*args, **kwargs)
            future.set_result(result)
        except Exception as e:
            future.set_exception(e)

        return future

    def shutdown(self, wait=True, *, cancel_futures=False):
        """Mark the executor as shut down."""
        self._shutdown = True


@pytest.fixture(autouse=True)
def reset_process_manager_between_tests():
    """Ensure process manager state does not leak between tests."""
    from quantum_calc.process_manager import shutdown_process_manager

    shutdown_process_manager()
    yield
    shutdown_process_manager()


@pytest.fixture(scope='function')
def app():
    """
    Create and configure a FastAPI application instance for testing.

    Yields:
        FastAPI: A configured FastAPI application instance in TESTING mode.
    """
    temp_dir = tempfile.mkdtemp(prefix='pyscf_test_')
    test_config = {
        'TESTING': True,
        'CALCULATIONS_DIR': temp_dir,
    }

    import services as services_module
    import quantum_calc.settings_manager as settings_manager_module
    from quantum_calc.settings_manager import SettingsManager
    from quantum_calc.process_manager import shutdown_process_manager

    test_settings_manager = SettingsManager(
        settings_file=os.path.join(temp_dir, "app-settings.json")
    )
    test_settings = test_settings_manager.get_default_settings().model_copy(
        update={'calculations_directory': temp_dir}
    )
    test_settings_manager.save_settings(test_settings)

    with (
        mock.patch.dict(os.environ, {'PYSCF_ENV': 'development', 'PYSCF_ENABLE_DEBUG_ENDPOINTS': 'true'}),
        mock.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor),
        mock.patch.object(settings_manager_module, "_settings_manager", test_settings_manager),
        mock.patch.multiple(
            services_module,
            _quantum_service=None,
            _pubchem_service=None,
            _smiles_service=None,
            _settings_service=None,
            _system_service=None,
        ),
    ):
        shutdown_process_manager()
        try:
            _app = create_fastapi_app(server_port=5000, test_config=test_config)
            yield _app
        finally:
            shutdown_process_manager()
            _clear_calculation_update_stream_hub()

    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture(scope='function')
def asgi_server() -> Generator[str, None, None]:
    """Run the real FastAPI ASGI app on a local Uvicorn server."""
    port = _get_free_port()
    temp_dir = tempfile.mkdtemp(prefix='pyscf_asgi_test_')
    test_config = {
        'TESTING': True,
        'CALCULATIONS_DIR': temp_dir,
    }

    import services as services_module
    import quantum_calc.settings_manager as settings_manager_module
    from app import create_app
    from quantum_calc.process_manager import shutdown_process_manager
    from quantum_calc.settings_manager import SettingsManager

    test_settings_manager = SettingsManager(
        settings_file=os.path.join(temp_dir, 'app-settings.json')
    )
    test_settings = test_settings_manager.get_default_settings().model_copy(
        update={'calculations_directory': temp_dir}
    )
    test_settings_manager.save_settings(test_settings)

    try:
        with (
            mock.patch.dict(
                os.environ,
                {'PYSCF_AUTH_TOKEN': 'test-token', 'PYSCF_ENV': 'development', 'PYSCF_ENABLE_DEBUG_ENDPOINTS': 'true'},
            ),
            mock.patch(
                'quantum_calc.process_manager.ProcessPoolExecutor',
                new=DummyExecutor,
            ),
            mock.patch.object(
                settings_manager_module,
                '_settings_manager',
                test_settings_manager,
            ),
            mock.patch.multiple(
                services_module,
                _quantum_service=None,
                _pubchem_service=None,
                _smiles_service=None,
                _settings_service=None,
                _system_service=None,
            ),
        ):
            shutdown_process_manager()
            uvicorn_app = create_app(server_port=port, test_config=test_config)
            config = uvicorn.Config(
                uvicorn_app,
                host='127.0.0.1',
                port=port,
                log_level='warning',
                access_log=False,
                lifespan='on',
                timeout_graceful_shutdown=1,
            )
            server = uvicorn.Server(config)
            thread = threading.Thread(target=server.run, daemon=True)
            thread.start()

            try:
                _wait_for_server_start(server, thread)
                yield f'http://127.0.0.1:{port}'
            finally:
                server.should_exit = True
                thread.join(timeout=5)
                if thread.is_alive():
                    server.force_exit = True
                    thread.join(timeout=1)
                shutdown_process_manager()
                _clear_calculation_update_stream_hub()
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture(scope='function')
def client(app):
    """
    Create a test client for making HTTP requests to the application.

    This fixture provides a FastAPI TestClient that can be used to simulate
    HTTP requests without running a real server. It's function-scoped to
    ensure each test gets a fresh client.

    Args:
        app: The FastAPI application instance (from app fixture).

    Returns:
        TestClient: A test client for the application.
    """
    with TestClient(app) as test_client:
        yield test_client


# ============================================================================
# Helper Classes and Utilities
# ============================================================================

# ============================================================================
# Test Data Fixtures
# ============================================================================

@pytest.fixture
def sample_h2_xyz():
    """
    Provide sample H2 molecule coordinates in valid XYZ format.

    Returns:
        str: XYZ coordinates for a hydrogen molecule with proper format:
             Line 1: number of atoms
             Line 2: comment line
             Lines 3+: atom symbol and coordinates
    """
    return "2\nHydrogen molecule\nH 0 0 0\nH 0 0 0.74"


@pytest.fixture
def valid_dft_params(sample_h2_xyz):
    """
    Provide valid DFT calculation parameters for testing.

    Args:
        sample_h2_xyz: Sample H2 molecule coordinates.

    Returns:
        dict: Valid parameters for a DFT calculation.
    """
    return {
        "name": "Test H2 DFT",
        "xyz": sample_h2_xyz,
        "calculation_method": "DFT",
        "basis_function": "sto-3g",
        "exchange_correlation": "b3lyp",
        "charges": 0,
        "spin": 0,
        "optimize_geometry": False
    }


@pytest.fixture
def valid_hf_params(sample_h2_xyz):
    """
    Provide valid Hartree-Fock calculation parameters for testing.

    Args:
        sample_h2_xyz: Sample H2 molecule coordinates.

    Returns:
        dict: Valid parameters for a Hartree-Fock calculation.
    """
    return {
        "name": "Test H2 HF",
        "xyz": sample_h2_xyz,
        "calculation_method": "HF",
        "basis_function": "sto-3g",
        "charges": 0,
        "spin": 0,
        "optimize_geometry": False
    }
