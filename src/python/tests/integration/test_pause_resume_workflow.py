"""
Integration tests for pause/resume workflow.

Tests the complete pause→resume→complete workflow including:
- Starting a calculation
- Pausing it mid-execution
- Resuming from checkpoint
- Completing the calculation
- Error cases for invalid operations
"""

import pytest
import os
import json
import time
from pathlib import Path


# ============================================================================
# Helper Functions
# ============================================================================

def wait_for_status(client, calc_id, expected_status, timeout=300, poll_interval=0.5):
    """
    Poll the calculation status until it matches expected_status or timeout.

    Args:
        client: Flask test client
        calc_id: Calculation ID
        expected_status: Expected status string (or list of strings)
        timeout: Maximum wait time in seconds
        poll_interval: Time between polls in seconds

    Returns:
        dict: The calculation data when status matches

    Raises:
        TimeoutError: If status doesn't match within timeout
    """
    if isinstance(expected_status, str):
        expected_status = [expected_status]

    start_time = time.time()
    last_status = None

    while time.time() - start_time < timeout:
        response = client.get(f'/api/quantum/calculations/{calc_id}')
        assert response.status_code == 200, f"Failed to get calculation: {response.get_json()}"

        data = response.get_json()
        calc = data['data']['calculation']
        current_status = calc['status']

        if current_status != last_status:
            print(f"Status changed: {last_status} → {current_status}")
            last_status = current_status

        if current_status in expected_status:
            return calc

        # If calculation ended in error, fail immediately
        if current_status == 'error':
            error_msg = calc.get('error_message', 'Unknown error')
            raise AssertionError(f"Calculation failed with error: {error_msg}")

        time.sleep(poll_interval)

    raise TimeoutError(
        f"Timeout waiting for status {expected_status}. "
        f"Last status: {last_status}, elapsed: {time.time() - start_time:.1f}s"
    )


def get_calculation_dir(app, calc_id):
    """Get the directory path for a calculation."""
    # Import settings_manager to get actual calculation directory
    from quantum_calc.settings_manager import get_settings_manager
    settings_manager = get_settings_manager()
    settings = settings_manager.load_settings()
    return Path(settings.calculations_directory) / calc_id


def check_file_exists(calc_dir, filename):
    """Check if a file exists in the calculation directory."""
    file_path = calc_dir / filename
    exists = file_path.exists()
    print(f"File check: {filename} → {'EXISTS' if exists else 'NOT FOUND'}")
    return exists


def read_json_file(calc_dir, filename):
    """Read and return the contents of a JSON file."""
    file_path = calc_dir / filename
    if not file_path.exists():
        return None
    with open(file_path, 'r') as f:
        return json.load(f)


# ============================================================================
# Test Class
# ============================================================================

class TestPauseResumeWorkflow:
    """
    Integration tests for pause/resume functionality.

    These tests use real ProcessPoolExecutor to test actual async behavior.
    They verify the complete workflow including file system state changes.
    """

    @pytest.fixture
    def quick_DFT_params(self):
        """
        Provide calculation parameters optimized for pause/resume testing.

        Uses ethanol molecule with cc-pVDZ basis and DFT (B3LYP) method with
        geometry optimization enabled to ensure the calculation runs long enough
        to test pause/resume functionality. This heavier calculation provides
        sufficient time to pause mid-execution.
        """
        # Ethanol molecule (C2H5OH) - larger and heavier than water
        ethanol_xyz = """9
Ethanol molecule
C   -0.8883  0.1670 -0.0273
C    0.4658 -0.5116 -0.0368
O    1.4311  0.3229  0.5867
H   -0.8487  1.1175 -0.5695
H   -1.6471 -0.4704 -0.4896
H   -1.1964  0.3978  0.9977
H    0.7920 -0.7224 -1.0597
H    0.4246 -1.4559  0.5138
H    1.4671  1.1550  0.0848"""

        return {
            "name": "Test Ethanol DFT for Pause/Resume",
            "xyz": ethanol_xyz,
            "calculation_method": "DFT",
            "functional": "b3lyp",  # Hybrid functional
            "basis_function": "cc-pvdz",  # Larger basis set for longer calculation
            "charges": 0,
            "spin": 0,
            "optimize_geometry": True,  # Enable to make calculation longer
            "cpu_cores": 1,
            "memory_mb": 2048
        }

    def test_pause_resume_full_workflow(self, client, app, quick_DFT_params):
        """
        GIVEN a running quantum calculation
        WHEN pause is requested, then resume is requested
        THEN the calculation pauses, resumes, and completes successfully

        This test verifies:
        1. Calculation starts and enters 'running' state
        2. Pause request triggers 'pausing' → 'paused' transition
        3. Pause state files are created (pause_state.json, calculation.chk)
        4. Resume request restarts calculation
        5. Calculation completes successfully
        6. Pause state files are cleaned up

        NOTE: This test uses real PySCF calculation (Ethanol with cc-pVDZ and DFT)
        which takes several seconds to complete, providing enough time to pause mid-execution.
        """
        # ARRANGE
        # Use real PySCF calculation for authentic pause behavior
        # Ethanol with cc-pVDZ and DFT is heavy enough to provide time to pause

        # ACT & ASSERT
        # Step 1: Submit calculation
        print("\n=== Step 1: Submitting calculation ===")
        response = client.post('/api/quantum/calculate', json=quick_DFT_params)
        assert response.status_code == 202

        data = response.get_json()
        assert data['success'] is True
        calc_id = data['data']['calculation']['id']
        print(f"Calculation ID: {calc_id}")

        calc_dir = get_calculation_dir(app, calc_id)

        # Step 2: Wait for calculation to start running
        print("\n=== Step 2: Waiting for 'running' status ===")
        calc = wait_for_status(client, calc_id, 'running', timeout=10)
        assert calc['status'] == 'running'

        # Step 3: Request pause while running
        print("\n=== Step 3: Requesting pause ===")
        # Wait briefly for geometry optimization to start, then pause
        # This gives the calculation time to enter a pausable state
        time.sleep(0.5)

        # Check current status before attempting to pause
        current_calc = client.get(f'/api/quantum/calculations/{calc_id}').get_json()['data']['calculation']
        if current_calc['status'] in ['completed', 'error']:
            pytest.skip(
                f"Calculation finished too quickly (status: {current_calc['status']}) to test pause/resume workflow. "
                "This is not a test failure, just means the calculation was too fast."
            )

        pause_response = client.post(f'/api/quantum/calculations/{calc_id}/pause')
        assert pause_response.status_code == 202

        pause_data = pause_response.get_json()
        assert pause_data['success'] is True
        assert 'pause' in pause_data['data']['message'].lower()

        # Verify .pause_requested flag file is created
        # (It may already be gone if pause happened very quickly)
        print("\n=== Checking for pause flag file ===")
        time.sleep(0.1)

        # Step 4: Wait for paused status
        print("\n=== Step 4: Waiting for 'paused' status ===")
        calc = wait_for_status(client, calc_id, 'paused', timeout=20)
        assert calc['status'] == 'paused'

        # Verify pause state files exist
        print("\n=== Step 5: Verifying pause state files ===")

        # Check if pause_state.json exists
        pause_state_exists = check_file_exists(calc_dir, 'pause_state.json')

        if not pause_state_exists:
            # If pause_state.json doesn't exist, the calculation might have completed too quickly
            # Check if calculation is still paused or already completed
            current_calc = client.get(f'/api/quantum/calculations/{calc_id}').get_json()['data']['calculation']
            print(f"WARNING: pause_state.json not found. Current status: {current_calc['status']}")

            # For this test to be meaningful, we need the calculation to actually pause
            # If it completed too quickly, we should skip the rest of the test
            if current_calc['status'] == 'completed':
                pytest.skip("Calculation completed too quickly to test pause/resume workflow. "
                           "This is not a test failure, just means the calculation was too fast.")

        assert pause_state_exists, "pause_state.json should exist after pausing"
        assert check_file_exists(calc_dir, 'calculation.chk'), "calculation.chk should exist"
        assert not check_file_exists(calc_dir, '.pause_requested'), ".pause_requested should be cleaned up"

        # Verify pause_state.json contents
        pause_state = read_json_file(calc_dir, 'pause_state.json')
        assert pause_state is not None
        assert 'paused_at' in pause_state
        assert 'checkpoint_exists' in pause_state
        print(f"Pause state: {pause_state}")

        # Step 6: Request resume
        print("\n=== Step 6: Requesting resume ===")
        resume_response = client.post(f'/api/quantum/calculations/{calc_id}/resume')
        assert resume_response.status_code == 202

        resume_data = resume_response.get_json()
        assert resume_data['success'] is True
        assert 'resume' in resume_data['data']['message'].lower()

        # Step 7: Wait for completion
        print("\n=== Step 7: Waiting for 'completed' status ===")
        calc = wait_for_status(client, calc_id, 'completed', timeout=600)
        assert calc['status'] == 'completed'

        # Step 8: Verify completion and cleanup
        print("\n=== Step 8: Verifying completion and cleanup ===")
        assert check_file_exists(calc_dir, 'results.json'), "results.json should exist"
        assert not check_file_exists(calc_dir, 'pause_state.json'), "pause_state.json should be deleted"

        # Verify results
        results = read_json_file(calc_dir, 'results.json')
        assert results is not None
        # For geometry optimization, energy is stored as 'scf_energy'
        assert 'scf_energy' in results or 'energy' in results
        energy = results.get('scf_energy', results.get('energy'))
        print(f"Final energy: {energy}")

        print("\n=== Test completed successfully ===")

    def test_pause_non_running_calculation(self, client, app, quick_DFT_params, mocker):
        """
        GIVEN a calculation that is not in 'running' state
        WHEN pause is requested
        THEN a 400 error is returned
        """
        # ARRANGE - Create a completed calculation
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.06
        mock_scf.mo_energy = [-0.5, 0.3]
        mock_scf.mo_occ = [2.0, 0.0]

        mocker.patch('quantum_calc.hf_calculator.gto.M', return_value=mock_mol)
        mocker.patch('quantum_calc.hf_calculator.scf.RHF', return_value=mock_scf)

        # Submit and wait for completion
        response = client.post('/api/quantum/calculate', json=quick_DFT_params)
        calc_id = response.get_json()['data']['calculation']['id']

        # Wait for completion
        calc = wait_for_status(client, calc_id, ['completed', 'error'], timeout=100)
        assert calc['status'] == 'completed'

        # ACT - Try to pause completed calculation
        pause_response = client.post(f'/api/quantum/calculations/{calc_id}/pause')

        # ASSERT
        assert pause_response.status_code == 400
        error_data = pause_response.get_json()
        assert error_data['success'] is False
        assert 'not running' in error_data['error'].lower() or 'cannot pause' in error_data['error'].lower()

    def test_resume_non_paused_calculation(self, client, app, quick_DFT_params, mocker):
        """
        GIVEN a calculation that is not in 'paused' state
        WHEN resume is requested
        THEN a 400 error is returned
        """
        # ARRANGE - Create a completed calculation
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.06
        mock_scf.mo_energy = [-0.5, 0.3]
        mock_scf.mo_occ = [2.0, 0.0]

        mocker.patch('quantum_calc.hf_calculator.gto.M', return_value=mock_mol)
        mocker.patch('quantum_calc.hf_calculator.scf.RHF', return_value=mock_scf)

        # Submit and wait for completion
        response = client.post('/api/quantum/calculate', json=quick_DFT_params)
        calc_id = response.get_json()['data']['calculation']['id']

        calc = wait_for_status(client, calc_id, ['completed', 'error'], timeout=100)
        assert calc['status'] == 'completed'

        # ACT - Try to resume completed calculation
        resume_response = client.post(f'/api/quantum/calculations/{calc_id}/resume')

        # ASSERT
        assert resume_response.status_code == 400
        error_data = resume_response.get_json()
        assert error_data['success'] is False
        assert 'not paused' in error_data['error'].lower() or 'cannot resume' in error_data['error'].lower()

    def test_pause_nonexistent_calculation(self, client):
        """
        GIVEN a non-existent calculation ID
        WHEN pause is requested
        THEN a 400 error is returned

        NOTE: The implementation returns 400 (ValidationError) instead of 404
        because the error occurs during validation before checking existence.
        """
        # ACT
        pause_response = client.post('/api/quantum/calculations/nonexistent_id_12345/pause')

        # ASSERT
        assert pause_response.status_code == 400
        error_data = pause_response.get_json()
        assert error_data['success'] is False
        assert 'not found' in error_data['error'].lower()

    def test_resume_nonexistent_calculation(self, client):
        """
        GIVEN a non-existent calculation ID
        WHEN resume is requested
        THEN a 400 error is returned

        NOTE: The implementation returns 400 (ValidationError) instead of 404
        because the error occurs during validation before checking existence.
        """
        # ACT
        resume_response = client.post('/api/quantum/calculations/nonexistent_id_12345/resume')

        # ASSERT
        assert resume_response.status_code == 400
        error_data = resume_response.get_json()
        assert error_data['success'] is False
        assert 'not found' in error_data['error'].lower()
