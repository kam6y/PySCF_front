"""
Central test configuration and fixtures for PySCF Front backend tests.

This module contains reusable pytest fixtures that are automatically available
to all test modules. Fixtures defined here follow best practices for Flask
application testing using the Application Factory pattern.
"""

import pytest
import tempfile
import os
from pathlib import Path
from concurrent.futures import Executor, Future

# Import application factory and socketio instance
from app import create_app, socketio


# ============================================================================
# Core Application Fixtures
# ============================================================================

@pytest.fixture(scope='session')
def app():
    """
    Create and configure a Flask application instance for testing.

    This is a session-scoped fixture, meaning it's created once per test session.
    The application is configured with TESTING=True and uses a temporary directory
    for calculations to ensure test isolation.

    Yields:
        Flask: A configured Flask application instance in TESTING mode.
    """
    # Create a temporary directory for test calculations
    temp_dir = tempfile.mkdtemp(prefix='pyscf_test_')

    # Test configuration
    test_config = {
        'TESTING': True,
        'CALCULATIONS_DIR': temp_dir,
        'WTF_CSRF_ENABLED': False,  # Disable CSRF for testing
        'SERVER_NAME': 'localhost:5000',  # Required for url_for() in tests
        # Disable WebSocket file watcher in tests
        'WEBSOCKET_WATCHER_ENABLED': False,
        # Use simple SocketIO configuration for testing
        'SOCKETIO': {
            'cors_allowed_origins': '*',
            'async_mode': 'threading',
            'logger': False,  # Reduce noise in test output
            'engineio_logger': False,
        }
    }

    # Create app using Application Factory with test configuration
    _app = create_app(server_port=5000, test_config=test_config)

    # Establish application context for the test session
    with _app.app_context():
        yield _app

    # Cleanup: Remove temporary directory
    import shutil
    try:
        shutil.rmtree(temp_dir)
    except Exception as e:
        print(f"Warning: Failed to cleanup temp directory {temp_dir}: {e}")


@pytest.fixture(scope='function')
def client(app):
    """
    Create a test client for making HTTP requests to the application.

    This fixture provides a Flask test client that can be used to simulate
    HTTP requests without running a real server. It's function-scoped to
    ensure each test gets a fresh client.

    Args:
        app: The Flask application instance (from app fixture).

    Returns:
        FlaskClient: A test client for the application.
    """
    return app.test_client()


@pytest.fixture(scope='function')
def socketio_client(app):
    """
    Create a SocketIO test client for testing WebSocket functionality.

    This fixture provides a SocketIO test client that can simulate WebSocket
    connections, emit events, and receive messages without a running server.

    Args:
        app: The Flask application instance (from app fixture).

    Returns:
        SocketIOTestClient: A test client for WebSocket communication.
    """
    return socketio.test_client(app, namespace=None)


@pytest.fixture
def runner(app):
    """
    Create a test runner for Flask CLI commands.

    This fixture is useful for testing custom Flask CLI commands.

    Args:
        app: The Flask application instance (from app fixture).

    Returns:
        FlaskCliRunner: A test runner for CLI commands.
    """
    return app.test_cli_runner()


# ============================================================================
# Helper Classes and Utilities
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

    def shutdown(self, wait=True):
        """Mark the executor as shut down."""
        self._shutdown = True


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
    Provide sample H2 molecule coordinates.

    Returns:
        str: XYZ coordinates for a hydrogen molecule.
    """
    return "H 0 0 0\nH 0 0 0.74"


@pytest.fixture
def sample_water_xyz():
    """
    Provide sample water molecule coordinates.

    Returns:
        str: XYZ coordinates for a water molecule.
    """
    return "O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692"


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
        "spin": 0
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
        "spin": 0
    }
