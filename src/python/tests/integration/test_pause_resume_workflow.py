"""
Integration tests for pause/resume workflow.

Tests the complete pause→resume→complete workflow including:
- Starting a calculation
- Pausing it mid-execution
- Resuming from checkpoint
- Completing the calculation
"""

import json
import time
from pathlib import Path

import pytest

from quantum_calc import get_current_settings
from quantum_calc._calculation_repository import CalculationRepository


# ============================================================================
# Helper Functions
# ============================================================================

def wait_for_status(client, calc_id, expected_status, timeout=300, poll_interval=0.5):
    """
    Poll the calculation status until it matches expected_status or timeout.

    Args:
        client: FastAPI test client
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
        assert response.status_code == 200, f"Failed to get calculation: {response.json()}"

        data = response.json()
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


def get_calculation_dir(calc_id):
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


class ControllablePauseResumeProcessManager:
    """Small fake process manager for deterministic pause/resume workflow tests."""

    def __init__(self):
        self.complete_submissions = False
        self.running_ids = set()
        self.queued_ids = set()

    def get_active_calculations(self):
        return list(self.running_ids)

    def get_queued_calculations(self):
        return list(self.queued_ids)

    def submit_calculation(self, calculation_id, parameters):
        repository = self._repository()
        calc_dir = self._calc_dir(repository, calculation_id)

        if parameters.get('resume_from_pause') or self.complete_submissions:
            self.running_ids.discard(calculation_id)
            repository.save_calculation_results(calc_dir, {'energy': -76.0})
            repository.save_calculation_status(calc_dir, 'completed')
            repository.delete_pause_state(calc_dir)
            return True, 'running', None

        self.running_ids.add(calculation_id)
        Path(calc_dir, 'calculation.chk').touch()
        repository.save_calculation_status(calc_dir, 'running')
        return True, 'running', None

    def pause_calculation(self, calculation_id):
        repository = self._repository()
        calc_dir = self._calc_dir(repository, calculation_id)
        if not Path(calc_dir).exists():
            raise ValueError(f"Calculation not found: {calculation_id}")

        status, _ = repository.read_calculation_status_details(calc_dir)
        if status != 'running':
            raise ValueError(f"Calculation is not running (status: {status})")

        flag_file = Path(calc_dir, '.pause_requested')
        flag_file.touch()
        Path(calc_dir, 'calculation.chk').touch()
        repository.save_pause_state(
            calc_dir,
            {
                'calculation_phase': 'scf_calculation',
                'checkpoint_exists': True,
            },
        )
        flag_file.unlink(missing_ok=True)
        repository.save_calculation_status(calc_dir, 'paused')
        self.running_ids.discard(calculation_id)
        return True

    def resume_calculation(self, calculation_id):
        repository = self._repository()
        calc_dir = self._calc_dir(repository, calculation_id)
        if not Path(calc_dir).exists():
            raise ValueError(f"Calculation not found: {calculation_id}")

        status, _ = repository.read_calculation_status_details(calc_dir)
        if status != 'paused':
            raise ValueError(f"Calculation is not paused (status: {status})")

        params = repository.read_calculation_parameters(calc_dir)
        if not params:
            raise ValueError("No calculation parameters found")

        params = {
            **params,
            'resume_from_pause': True,
            'pause_state': repository.load_pause_state(calc_dir),
        }
        self.submit_calculation(calculation_id, params)
        return {'calculation_id': calculation_id, 'status': 'completed'}

    def _repository(self):
        settings = get_current_settings()
        return CalculationRepository(base_dir=settings.calculations_directory)

    def _calc_dir(self, repository, calculation_id):
        return str(repository.resolve_calculation_path(calculation_id))


# ============================================================================
# Test Class
# ============================================================================

class TestPauseResumeWorkflow:
    """
    Integration tests for pause/resume functionality.

    These tests use a controllable fake process manager to avoid depending on
    multiprocessing or real PySCF runtime. They verify the API workflow and
    persisted file system state changes.
    """

    @pytest.fixture
    def process_manager(self, mocker):
        manager = ControllablePauseResumeProcessManager()
        mocker.patch(
            'services.calculation_service_context.get_process_manager',
            return_value=manager,
        )
        return manager

    @pytest.fixture
    def quick_DFT_params(self):
        """
        Provide calculation parameters optimized for pause/resume testing.

        Uses ethanol-like DFT parameters. The process manager is faked in these
        tests, so the molecule size no longer controls timing.
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
            "exchange_correlation": "b3lyp",  # Hybrid functional (correct field name)
            "basis_function": "cc-pvdz",  # Larger basis set for longer calculation
            "charges": 0,
            "spin": 0,
            "optimize_geometry": True,  # Enable to make calculation longer
            "cpu_cores": 1,
            "memory_mb": 2048
        }

    @pytest.fixture(autouse=True)
    def mock_system_resources(self, mocker):
        """
        Mock system resources to ensure tests don't fail due to CI load.
        
        This prevents 'System CPU usage exceeds limit' errors when running
        tests on busy CI runners.
        """
        # Mock psutil to return low CPU usage
        mocker.patch('quantum_calc.resource_manager.psutil.cpu_percent', return_value=10.0)
        
        # Mock virtual memory to ensure sufficient memory is available
        mock_memory = mocker.Mock()
        mock_memory.total = 16 * 1024 * 1024 * 1024  # 16 GB
        mock_memory.available = 8 * 1024 * 1024 * 1024  # 8 GB
        mock_memory.percent = 50.0
        mocker.patch('quantum_calc.resource_manager.psutil.virtual_memory', return_value=mock_memory)

    def test_pause_resume_full_workflow(self, client, quick_DFT_params, process_manager):
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

        NOTE: This test uses a controllable fake process manager so it does not
        depend on real PySCF runtime or sandbox multiprocessing support.
        """
        # ARRANGE
        process_manager.complete_submissions = False

        # ACT & ASSERT
        # Step 1: Submit calculation
        print("\n=== Step 1: Submitting calculation ===")
        response = client.post('/api/quantum/calculate', json=quick_DFT_params)
        assert response.status_code == 202

        data = response.json()
        assert data['success'] is True
        calc_id = data['data']['calculation']['id']
        print(f"Calculation ID: {calc_id}")

        calc_dir = get_calculation_dir(calc_id)

        # Step 2: Wait for calculation to start running
        print("\n=== Step 2: Waiting for 'running' status ===")
        calc = wait_for_status(client, calc_id, 'running', timeout=10)
        assert calc['status'] == 'running'

        # Step 3: Request pause while running
        print("\n=== Step 3: Requesting pause ===")
        pause_response = client.post(f'/api/quantum/calculations/{calc_id}/pause')
        assert pause_response.status_code == 202

        pause_data = pause_response.json()
        assert pause_data['success'] is True
        assert 'pause' in pause_data['data']['message'].lower()

        # Step 4: Wait for paused status
        print("\n=== Step 4: Waiting for 'paused' status ===")
        calc = wait_for_status(client, calc_id, 'paused', timeout=20)
        assert calc['status'] == 'paused'

        # Verify pause state files exist
        print("\n=== Step 5: Verifying pause state files ===")

        # Check if pause_state.json exists
        pause_state_exists = check_file_exists(calc_dir, 'pause_state.json')
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

        resume_data = resume_response.json()
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
