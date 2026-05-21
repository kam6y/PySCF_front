"""
Central test configuration and fixtures for PySCF Front backend tests.

This module contains reusable pytest fixtures that are automatically available
to all test modules.
"""

import os
import logging
import tempfile
from concurrent.futures import Executor, Future
from unittest import mock

import pytest
from fastapi.testclient import TestClient

from app import create_fastapi_app

logger = logging.getLogger(__name__)


# ============================================================================
# Core Application Fixtures
# ============================================================================


@pytest.fixture(autouse=True)
def _configure_application():
    """Disable pytest-flask app.config mutation for FastAPI tests."""


@pytest.fixture(autouse=True)
def _monkeypatch_response_class():
    """Disable pytest-flask response_class patching for FastAPI tests."""


@pytest.fixture(autouse=True)
def _push_request_context():
    """Disable pytest-flask request-context setup for FastAPI tests."""


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
        'WEBSOCKET_WATCHER_ENABLED': False,
        'SOCKETIO': {
            'cors_allowed_origins': ['http://127.0.0.1:*', 'http://localhost:*', 'file://'],
            'logger': False,
            'engineio_logger': False,
        },
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
        mock.patch.dict(os.environ, {'PYSCF_ENV': 'development'}),
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
        _app = create_fastapi_app(server_port=5000, test_config=test_config)
        yield _app
        shutdown_process_manager()

    import shutil
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture(scope='function')
def client(app):
    """
    Create a test client for making HTTP requests to the application.

    This fixture provides a Flask test client that can be used to simulate
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

@pytest.fixture
def dummy_executor():
    """
    Provide a DummyExecutor instance for testing async workflows.

    Returns:
        DummyExecutor: A synchronous executor for testing.
    """
    return DummyExecutor()


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
def sample_water_xyz():
    """
    Provide sample water molecule coordinates in valid XYZ format.

    Returns:
        str: XYZ coordinates for a water molecule with proper format:
             Line 1: number of atoms
             Line 2: comment line
             Lines 3+: atom symbol and coordinates
    """
    return "3\nWater molecule\nO 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692"


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
