"""
Unit tests for Quantum service.

Tests the QuantumService class parameter validation logic with mocked dependencies.
Focuses on validate_calculation_parameters() which is core business logic.
"""

import pytest

from quantum_calc import CalculationError
from quantum_calc._calculation_repository import CalculationRepository
from services.quantum_service import QuantumService
from services.exceptions import ServiceError, ValidationError


# ============================================================================
# validate_calculation_parameters() - HF Method Tests
# ============================================================================

def test_validate_hf_method_valid_params():
    """
    GIVEN valid HF calculation parameters
    WHEN validate_calculation_parameters is called
    THEN it should return None (validation passes)
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': 'sto-3g',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None  # No validation errors


def test_validate_hf_with_exchange_correlation_warning(caplog):
    """
    GIVEN HF parameters with exchange_correlation
    WHEN validate_calculation_parameters is called
    THEN it should reject the request with an error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0,
        'exchange_correlation': 'PBE0'  # Not applicable to HF
    }

    # ACT
    result = service.validate_calculation_parameters(params)

    # ASSERT
    assert result is not None  # Should fail validation
    assert 'exchange_correlation' in result
    assert 'not applicable' in result.lower()
    assert 'DFT' in result or 'TDDFT' in result


def test_validate_hf_with_tddft_params_warning(caplog):
    """
    GIVEN HF parameters with TDDFT-specific parameters
    WHEN validate_calculation_parameters is called
    THEN it should reject the request with an error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': 'sto-3g',
        'charges': 0,
        'spin': 0,
        'tddft_nstates': 20  # Not applicable to HF
    }

    # ACT
    result = service.validate_calculation_parameters(params)

    # ASSERT
    assert result is not None  # Should fail validation
    assert 'tddft_nstates' in result
    assert 'not applicable' in result.lower()
    assert 'TDDFT' in result


# ============================================================================
# validate_calculation_parameters() - DFT Method Tests
# ============================================================================

def test_validate_dft_method_valid_params():
    """
    GIVEN valid DFT calculation parameters with XC functional
    WHEN validate_calculation_parameters is called
    THEN it should return None (validation passes)
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'DFT',
        'basis_function': '6-31G(d)',
        'exchange_correlation': 'B3LYP',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None


def test_validate_dft_missing_exchange_correlation():
    """
    GIVEN DFT parameters without exchange_correlation functional
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'DFT',
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
        # Missing exchange_correlation
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'exchange-correlation functional' in result.lower()


@pytest.mark.parametrize("xc_functional", [
    'B3LYP',
    'PBE0',
    'M06-2X',
    'CAM-B3LYP',
    'wB97X-D'
])
def test_validate_dft_various_functionals(xc_functional):
    """
    GIVEN DFT parameters with various valid XC functionals
    WHEN validate_calculation_parameters is called
    THEN it should return None for all
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'DFT',
        'exchange_correlation': xc_functional,
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None


# ============================================================================
# validate_calculation_parameters() - TDDFT Method Tests
# ============================================================================

def test_validate_tddft_valid_params():
    """
    GIVEN valid TDDFT calculation parameters
    WHEN validate_calculation_parameters is called
    THEN it should return None
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'TDDFT',
        'basis_function': '6-31G(d)',
        'exchange_correlation': 'B3LYP',
        'tddft_nstates': 10,
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None


def test_validate_tddft_missing_exchange_correlation():
    """
    GIVEN TDDFT parameters without exchange_correlation
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'TDDFT',
        'tddft_nstates': 5,
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'exchange-correlation functional' in result.lower()


def test_validate_tddft_missing_nstates():
    """
    GIVEN TDDFT parameters without tddft_nstates
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'TDDFT',
        'exchange_correlation': 'B3LYP',
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'tddft_nstates' in result.lower()


def test_validate_tddft_nstates_zero():
    """
    GIVEN TDDFT parameters with tddft_nstates = 0
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'TDDFT',
        'exchange_correlation': 'B3LYP',
        'tddft_nstates': 0,  # Invalid: must be > 0
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'greater than 0' in result.lower()


# ============================================================================
# validate_calculation_parameters() - CASCI Method Tests
# ============================================================================

def test_validate_casci_valid_params():
    """
    GIVEN valid CASCI calculation parameters
    WHEN validate_calculation_parameters is called
    THEN it should return None
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASCI',
        'basis_function': '6-31G',
        'ncas': 4,
        'nelecas': 4,
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None


def test_validate_casci_missing_ncas():
    """
    GIVEN CASCI parameters without ncas
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASCI',
        'nelecas': 4,
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'ncas' in result.lower()


def test_validate_casci_missing_nelecas():
    """
    GIVEN CASCI parameters without nelecas
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASCI',
        'ncas': 6,
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'nelecas' in result.lower()


def test_validate_casci_nelecas_exceeds_limit():
    """
    GIVEN CASCI parameters where nelecas > 2 * ncas
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASCI',
        'ncas': 4,
        'nelecas': 10,  # Exceeds 2 * ncas (8)
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'cannot exceed 2 * ncas' in result.lower()


def test_validate_casci_ncas_zero():
    """
    GIVEN CASCI parameters with ncas = 0
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASCI',
        'ncas': 0,  # Invalid
        'nelecas': 0,
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'greater than 0' in result.lower()


# ============================================================================
# validate_calculation_parameters() - CASSCF Method Tests
# ============================================================================

def test_validate_casscf_valid_params():
    """
    GIVEN valid CASSCF calculation parameters
    WHEN validate_calculation_parameters is called
    THEN it should return None
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASSCF',
        'basis_function': '6-31G',
        'ncas': 6,
        'nelecas': 6,
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None


def test_validate_casscf_same_constraints_as_casci():
    """
    GIVEN CASSCF parameters with invalid ncas/nelecas
    WHEN validate_calculation_parameters is called
    THEN it should return error message (same constraints as CASCI)
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASSCF',
        'ncas': 3,
        'nelecas': 10,  # Exceeds 2 * ncas
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'cannot exceed 2 * ncas' in result.lower()


# ============================================================================
# validate_calculation_parameters() - General Parameter Tests
# ============================================================================

def test_validate_negative_spin():
    """
    GIVEN parameters with negative spin
    WHEN validate_calculation_parameters is called
    THEN it should return OpenAPI-aligned bounds error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': 'sto-3g',
        'charges': 0,
        'spin': -1  # Invalid: cannot be negative
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is not None
    assert 'spin' in result
    assert 'below minimum' in result.lower()


def test_validate_high_charge_rejected():
    """
    GIVEN parameters with charge above the OpenAPI maximum
    WHEN validate_calculation_parameters is called
    THEN it should return error message
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': 'sto-3g',
        'charges': 11,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)

    # ASSERT
    assert result is not None
    assert 'charges' in result
    assert 'exceeds maximum' in result.lower()


@pytest.mark.parametrize("method", ['HF', 'DFT', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF'])
def test_validate_all_methods_accept_basic_params(method):
    """
    GIVEN basic parameters for any calculation method
    WHEN validate_calculation_parameters is called
    THEN it should not crash (may return errors for method-specific params)
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': method,
        'basis_function': 'sto-3g',
        'charges': 0,
        'spin': 0,
        # Add method-specific params to make them valid
        'exchange_correlation': 'B3LYP' if method in ['DFT', 'TDDFT'] else None,
        'tddft_nstates': 5 if method == 'TDDFT' else None,
        'ncas': 4 if method in ['CASCI', 'CASSCF'] else None,
        'nelecas': 4 if method in ['CASCI', 'CASSCF'] else None,
    }
    
    # ACT - Should not raise exception
    result = service.validate_calculation_parameters(params)
    
    # ASSERT - Just ensure it doesn't crash
    # Result can be None (valid) or error message (invalid)
    assert result is None or isinstance(result, str)


# ============================================================================
# get_supported_parameters() Tests
# ============================================================================

def test_get_supported_parameters_success(mocker):
    """
    GIVEN get_all_supported_parameters returns data
    WHEN get_supported_parameters is called
    THEN it should return the parameters dict
    """
    # ARRANGE
    mock_params = {
        'calculation_methods': ['HF', 'DFT', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF'],
        'basis_sets': ['STO-3G', '6-31G', '6-31G(d)', 'cc-pVDZ'],
        'exchange_correlations': ['B3LYP', 'PBE0', 'M06-2X']
    }
    
    mocker.patch('services.quantum_service.get_all_supported_parameters', return_value=mock_params)
    
    service = QuantumService()
    
    # ACT
    result = service.get_supported_parameters()
    
    # ASSERT
    assert result == mock_params
    assert 'calculation_methods' in result
    assert 'basis_sets' in result


def test_get_supported_parameters_error(mocker):
    """
    GIVEN get_all_supported_parameters raises an exception
    WHEN get_supported_parameters is called
    THEN it should raise ServiceError
    """
    # ARRANGE
    mocker.patch(
        'services.quantum_service.get_all_supported_parameters',
        side_effect=RuntimeError('Module error')
    )
    
    service = QuantumService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError, match="Failed to retrieve supported parameters"):
        service.get_supported_parameters()


def test_start_calculation_preserves_terminal_status_written_during_submit(tmp_path, mocker):
    """
    GIVEN a process manager completes a calculation during submit
    WHEN start_calculation returns
    THEN the terminal status is not overwritten by the initial running status
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    params = {
        "name": "Fast Calc",
        "created_at": "2026-05-20T00:00:00",
        "calculation_method": "HF",
        "basis_function": "sto-3g",
        "charges": 0,
        "spin": 0,
    }

    def complete_during_submit(calculation_id, _params):
        calc_dir = tmp_path / calculation_id
        service.repository.save_calculation_results(
            str(calc_dir),
            {"energy": -1.0, "success": True},
        )
        service.repository.save_calculation_status(str(calc_dir), "completed")
        return True, "running", None

    process_manager = mocker.Mock()
    process_manager.get_active_calculations.return_value = []
    process_manager.get_queued_calculations.return_value = []
    process_manager.submit_calculation.side_effect = complete_during_submit
    mocker.patch("services.quantum_service.get_process_manager", return_value=process_manager)

    response = service.start_calculation(params)
    calc_dir = tmp_path / response["id"]

    assert response["status"] == "completed"
    assert service.repository.read_calculation_status_details(str(calc_dir)) == (
        "completed",
        None,
    )
    assert service.repository.read_calculation_results(str(calc_dir)) == {
        "energy": -1.0,
        "success": True,
    }


def test_start_calculation_does_not_recover_new_pending_before_submit(
    tmp_path,
    mocker,
):
    """
    GIVEN stale recovery sees no active or queued calculations
    WHEN start_calculation creates a new calculation and submit succeeds
    THEN the new pending calculation is not marked as restart-interrupted
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    params = {
        "name": "New Calc",
        "created_at": "2026-05-20T00:00:00",
        "calculation_method": "HF",
        "basis_function": "sto-3g",
        "charges": 0,
        "spin": 0,
    }

    process_manager = mocker.Mock()
    process_manager.get_active_calculations.return_value = []
    process_manager.get_queued_calculations.return_value = []
    process_manager.submit_calculation.return_value = (True, "running", None)
    mocker.patch("services.quantum_service.get_process_manager", return_value=process_manager)

    response = service.start_calculation(params)
    calc_dir = tmp_path / response["id"]

    assert response["status"] == "running"
    assert service.repository.read_calculation_status_details(str(calc_dir)) == (
        "running",
        None,
    )
    assert service.repository.read_calculation_results(str(calc_dir)) != {
        "error": QuantumService.RESTART_INTERRUPTED_MESSAGE,
    }


def test_recover_stale_non_terminal_calculations_marks_only_stale_as_error(
    tmp_path,
    mocker,
):
    """
    GIVEN persisted non-terminal statuses without active or queued work
    WHEN stale calculations are recovered
    THEN only those non-terminal calculations are marked as error
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    stale_statuses = ["pending", "running", "waiting", "pausing"]
    terminal_statuses = ["completed", "paused"]

    for status in stale_statuses + terminal_statuses:
        calc_id = f"calc-{status}"
        calc_dir = tmp_path / calc_id
        calc_dir.mkdir()
        service.repository.save_calculation_parameters(
            str(calc_dir),
            {"name": calc_id, "created_at": "2026-05-20T00:00:00"},
        )
        if status != "pending":
            service.repository.save_calculation_status(str(calc_dir), status)

    process_manager = mocker.Mock()
    process_manager.get_active_calculations.return_value = []
    process_manager.get_queued_calculations.return_value = []

    service._recover_stale_non_terminal_calculations(process_manager)

    for status in stale_statuses:
        calc_dir = tmp_path / f"calc-{status}"
        assert service.repository.read_calculation_status_details(str(calc_dir)) == (
            "error",
            None,
        )
        assert service.repository.read_calculation_results(str(calc_dir)) == {
            "error": "Calculation interrupted because the backend was restarted.",
        }

    assert service.repository.read_calculation_status_details(
        str(tmp_path / "calc-completed")
    ) == ("completed", None)
    assert service.repository.read_calculation_status_details(
        str(tmp_path / "calc-paused")
    ) == ("paused", None)


def test_recover_stale_non_terminal_calculations_keeps_active_and_queued(
    tmp_path,
    mocker,
):
    """
    GIVEN non-terminal calculations still known by the process manager
    WHEN stale recovery runs
    THEN active and queued calculations keep their current statuses
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    for calc_id, status in {
        "active-calc": "running",
        "queued-calc": "waiting",
        "stale-calc": "pausing",
    }.items():
        calc_dir = tmp_path / calc_id
        calc_dir.mkdir()
        service.repository.save_calculation_parameters(
            str(calc_dir),
            {"name": calc_id, "created_at": "2026-05-20T00:00:00"},
        )
        service.repository.save_calculation_status(str(calc_dir), status)

    process_manager = mocker.Mock()
    process_manager.get_active_calculations.return_value = ["active-calc"]
    process_manager.get_queued_calculations.return_value = ["queued-calc"]

    service._recover_stale_non_terminal_calculations(process_manager)

    assert service.repository.read_calculation_status_details(
        str(tmp_path / "active-calc")
    ) == ("running", None)
    assert service.repository.read_calculation_status_details(
        str(tmp_path / "queued-calc")
    ) == ("waiting", None)
    assert service.repository.read_calculation_status_details(
        str(tmp_path / "stale-calc")
    ) == ("error", None)


def test_recover_stale_non_terminal_calculations_preserves_existing_results(
    tmp_path,
    mocker,
):
    """
    GIVEN a stale non-terminal calculation already has valid results
    WHEN stale recovery marks the calculation as error
    THEN the existing results file should not be overwritten
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_dir = tmp_path / "stale-with-results"
    calc_dir.mkdir()
    expected_results = {
        "energy": -75.0,
        "success": True,
        "mulliken_charges": [0.1, -0.1],
    }
    service.repository.save_calculation_parameters(
        str(calc_dir),
        {"name": "stale-with-results", "created_at": "2026-05-20T00:00:00"},
    )
    service.repository.save_calculation_status(str(calc_dir), "running")
    service.repository.save_calculation_results(str(calc_dir), expected_results)

    process_manager = mocker.Mock()
    process_manager.get_active_calculations.return_value = []
    process_manager.get_queued_calculations.return_value = []

    service._recover_stale_non_terminal_calculations(process_manager)

    assert service.repository.read_calculation_status_details(str(calc_dir)) == (
        "error",
        None,
    )
    assert service.repository.read_calculation_results(str(calc_dir)) == expected_results


@pytest.mark.parametrize("status", ["pending", "running", "waiting", "pausing"])
def test_delete_calculation_rejects_non_terminal_status(tmp_path, mocker, status):
    """
    GIVEN a calculation directory has a non-terminal status
    WHEN delete_calculation is called
    THEN deletion is rejected and the directory remains
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_id = f"calc-{status}"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    service.repository.save_calculation_parameters(
        str(calc_dir),
        {"name": "Queued Calc", "created_at": "2026-05-20T00:00:00"},
    )
    if status != "pending":
        service.repository.save_calculation_status(str(calc_dir), status)

    process_manager = mocker.Mock()
    process_manager.is_running.return_value = False
    mocker.patch("services.quantum_service.get_process_manager", return_value=process_manager)

    with pytest.raises(ValidationError, match="Cannot delete calculation"):
        service.delete_calculation(calc_id)

    assert calc_dir.is_dir()


@pytest.mark.parametrize(
    ("grid_size", "isovalue_pos", "isovalue_neg", "match"),
    [
        (39, None, None, "grid_size"),
        (121, None, None, "grid_size"),
        (80, 0.0009, None, "isovalue_pos"),
        (80, 0.1001, None, "isovalue_pos"),
        (80, None, -0.1001, "isovalue_neg"),
        (80, None, -0.0009, "isovalue_neg"),
    ],
)
def test_generate_orbital_cube_rejects_out_of_range_parameters(
    tmp_path,
    mocker,
    grid_size,
    isovalue_pos,
    isovalue_neg,
    match,
):
    """
    GIVEN orbital CUBE parameters outside the OpenAPI contract
    WHEN generate_orbital_cube is called
    THEN service-layer validation rejects them before generation
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_id = "completed-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    service.repository.save_calculation_status(str(calc_dir), "completed")
    orbital_generator = mocker.patch("services.quantum_service.MolecularOrbitalGenerator")

    with pytest.raises(ValidationError, match=match):
        service.generate_orbital_cube(
            calc_id,
            5,
            grid_size=grid_size,
            isovalue_pos=isovalue_pos,
            isovalue_neg=isovalue_neg,
        )

    orbital_generator.assert_not_called()


def test_generate_orbital_cube_accepts_contract_boundary_parameters(tmp_path, mocker):
    """
    GIVEN orbital CUBE parameters on the OpenAPI contract boundaries
    WHEN generate_orbital_cube is called
    THEN those values are passed through to the orbital generator
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_id = "completed-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    service.repository.save_calculation_status(str(calc_dir), "completed")

    generator = mocker.Mock()
    generator.validate_calculation.return_value = True
    generator.generate_cube_file.return_value = {
        "generation_params": {"file_size_kb": 1.0},
        "cached": False,
    }
    orbital_generator = mocker.patch(
        "services.quantum_service.MolecularOrbitalGenerator",
        return_value=generator,
    )

    result = service.generate_orbital_cube(
        calc_id,
        5,
        grid_size=40,
        isovalue_pos=0.001,
        isovalue_neg=-0.001,
    )

    orbital_generator.assert_called_once_with(str(calc_dir))
    generator.generate_cube_file.assert_called_once_with(
        orbital_index=5,
        grid_size=40,
        isovalue_pos=0.001,
        isovalue_neg=-0.001,
        return_content=True,
        save_to_disk=True,
    )
    assert result["generation_params"]["file_size_kb"] == 1.0


def test_generate_orbital_cube_invalid_orbital_index_returns_validation_error(
    tmp_path, mocker
):
    """
    GIVEN the orbital generator rejects an unavailable orbital index
    WHEN generate_orbital_cube is called
    THEN the service maps it to a validation error instead of a 500 error
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_id = "completed-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    service.repository.save_calculation_status(str(calc_dir), "completed")

    generator = mocker.Mock()
    generator.validate_calculation.return_value = True
    generator.generate_cube_file.side_effect = CalculationError(
        "Invalid orbital index: 12. Available range: 0-5"
    )
    mocker.patch(
        "services.quantum_service.MolecularOrbitalGenerator",
        return_value=generator,
    )

    with pytest.raises(ValidationError, match="Invalid orbital index"):
        service.generate_orbital_cube(calc_id, 12)


def test_generate_orbital_cube_generation_error_stays_service_error(tmp_path, mocker):
    """
    GIVEN the orbital generator fails for a non-validation reason
    WHEN generate_orbital_cube is called
    THEN the service keeps the error mapped to a service error
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_id = "completed-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    service.repository.save_calculation_status(str(calc_dir), "completed")

    generator = mocker.Mock()
    generator.validate_calculation.return_value = True
    generator.generate_cube_file.side_effect = CalculationError(
        "Failed to generate CUBE file"
    )
    mocker.patch(
        "services.quantum_service.MolecularOrbitalGenerator",
        return_value=generator,
    )

    with pytest.raises(ServiceError, match="Failed to generate CUBE file") as exc_info:
        service.generate_orbital_cube(calc_id, 2)

    assert exc_info.value.status_code == 500


def test_resume_calculation_returns_updated_waiting_status(tmp_path, mocker):
    """
    GIVEN a paused calculation is resumed and queued
    WHEN resume_calculation returns
    THEN the response calculation reflects waiting status and reason
    """
    service = QuantumService()
    service.repository = CalculationRepository(base_dir=str(tmp_path))

    calc_id = "paused-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    params = {
        "name": "Paused Calc",
        "created_at": "2026-05-20T00:00:00",
        "calculation_method": "HF",
    }
    service.repository.save_calculation_parameters(str(calc_dir), params)
    service.repository.save_calculation_status(str(calc_dir), "paused")

    def resume_and_queue(calculation_id):
        service.repository.save_calculation_status(
            str(calc_dir), "waiting", "All slots are busy"
        )
        return {
            "calculation_id": calculation_id,
            "status": "waiting",
            "waiting_reason": "All slots are busy",
        }

    process_manager = mocker.Mock()
    process_manager.resume_calculation.side_effect = resume_and_queue
    mocker.patch("services.quantum_service.get_process_manager", return_value=process_manager)

    response = service.resume_calculation(calc_id)

    assert response["calculation"]["status"] == "waiting"
    assert response["calculation"]["waitingReason"] == "All slots are busy"


# ============================================================================
# Edge Cases and Boundary Tests
# ============================================================================

@pytest.mark.parametrize("ncas,nelecas,should_pass", [
    (4, 4, True),   # Valid: exactly half filled
    (4, 8, True),   # Valid: fully filled
    (4, 2, True),   # Valid: partially filled
    (4, 9, False),  # Invalid: exceeds 2*ncas
    (5, 10, True),  # Valid: exactly at limit
    (5, 11, False), # Invalid: one over limit
])
def test_validate_casci_electron_orbital_relationships(ncas, nelecas, should_pass):
    """
    GIVEN various ncas/nelecas combinations
    WHEN validate_calculation_parameters is called for CASCI
    THEN it should validate electron-orbital relationship correctly
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'CASCI',
        'ncas': ncas,
        'nelecas': nelecas,
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    if should_pass:
        assert result is None, f"Expected validation to pass for ncas={ncas}, nelecas={nelecas}"
    else:
        assert result is not None, f"Expected validation to fail for ncas={ncas}, nelecas={nelecas}"
        assert 'cannot exceed' in result.lower()


def test_validate_params_empty_dict():
    """
    GIVEN an empty parameters dictionary
    WHEN validate_calculation_parameters is called
    THEN it should not crash (returns None as no method specified)
    """
    # ARRANGE
    service = QuantumService()
    params = {}
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    # Should return None as there's no calculation_method to validate
    assert result is None
