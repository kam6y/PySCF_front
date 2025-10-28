"""
Sample test to verify that the test infrastructure is working correctly.

This test module validates that all fixtures defined in conftest.py
are properly configured and functional.
"""

import pytest


def test_app_fixture(app):
    """
    GIVEN the app fixture
    WHEN the application is accessed
    THEN it should be in TESTING mode with correct configuration
    """
    assert app.config['TESTING'] is True
    assert 'CALCULATIONS_DIR' in app.config
    assert app.config['WTF_CSRF_ENABLED'] is False


def test_client_fixture(client):
    """
    GIVEN the client fixture
    WHEN a simple GET request is made
    THEN it should return a response
    """
    response = client.get('/health')
    assert response is not None
    # Health endpoint should return 200
    assert response.status_code == 200


def test_socketio_client_fixture(socketio_client):
    """
    GIVEN the socketio_client fixture
    WHEN the client is accessed
    THEN it should be properly initialized and connected
    """
    assert socketio_client is not None
    assert socketio_client.is_connected()


def test_sample_data_fixtures(sample_h2_xyz, sample_water_xyz):
    """
    GIVEN sample molecule fixtures
    WHEN the fixtures are accessed
    THEN they should contain valid XYZ coordinate data
    """
    # H2 should have 4 lines (atom count + comment + 2 atoms)
    assert 'H' in sample_h2_xyz
    h2_lines = [line for line in sample_h2_xyz.split('\n') if line.strip()]
    assert len(h2_lines) == 4

    # Water should have 5 lines (atom count + comment + 3 atoms)
    assert 'O' in sample_water_xyz
    assert 'H' in sample_water_xyz
    water_lines = [line for line in sample_water_xyz.split('\n') if line.strip()]
    assert len(water_lines) == 5


def test_valid_params_fixtures(valid_dft_params, valid_hf_params):
    """
    GIVEN valid parameter fixtures
    WHEN the fixtures are accessed
    THEN they should contain all required calculation parameters
    """
    # Check DFT parameters
    assert 'name' in valid_dft_params
    assert 'xyz' in valid_dft_params
    assert valid_dft_params['calculation_method'] == 'DFT'
    assert 'basis_function' in valid_dft_params
    assert 'exchange_correlation' in valid_dft_params

    # Check HF parameters
    assert 'name' in valid_hf_params
    assert 'xyz' in valid_hf_params
    assert valid_hf_params['calculation_method'] == 'HF'
    assert 'basis_function' in valid_hf_params


def test_dummy_executor_fixture(dummy_executor):
    """
    GIVEN the dummy_executor fixture
    WHEN a task is submitted
    THEN it should execute synchronously and return the result
    """
    def sample_task(x, y):
        return x + y

    future = dummy_executor.submit(sample_task, 2, 3)
    assert future.done()
    assert future.result() == 5


def test_dummy_executor_exception_handling(dummy_executor):
    """
    GIVEN the dummy_executor fixture
    WHEN a task raises an exception
    THEN the exception should be captured in the Future
    """
    def failing_task():
        raise ValueError("Test error")

    future = dummy_executor.submit(failing_task)
    assert future.done()

    with pytest.raises(ValueError, match="Test error"):
        future.result()
