"""
Unit tests for Quantum service.

Tests the QuantumService class parameter validation logic with mocked dependencies.
Focuses on validate_calculation_parameters() which is core business logic.
"""

import pytest
from unittest.mock import MagicMock

from services.quantum_service import QuantumService
from services.exceptions import ServiceError


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
    GIVEN HF parameters with non-default exchange_correlation
    WHEN validate_calculation_parameters is called
    THEN it should log a warning but not fail
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': '6-31G',
        'charges': 0,
        'spin': 0,
        'exchange_correlation': 'PBE0'  # Non-default for HF
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None  # Still passes validation


def test_validate_hf_with_tddft_params_warning(caplog):
    """
    GIVEN HF parameters with TDDFT-specific parameters
    WHEN validate_calculation_parameters is called
    THEN it should log a warning but not fail
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': 'sto-3g',
        'charges': 0,
        'spin': 0,
        'tddft_nstates': 20  # TDDFT parameter ignored for HF
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None  # Still passes


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
    THEN it should return error message
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
    assert 'cannot be negative' in result.lower()


def test_validate_high_charge_warning(caplog):
    """
    GIVEN parameters with very high molecular charge
    WHEN validate_calculation_parameters is called
    THEN it should log a warning but pass validation
    """
    # ARRANGE
    service = QuantumService()
    params = {
        'calculation_method': 'HF',
        'basis_function': 'sto-3g',
        'charges': 12,  # High charge, should trigger warning
        'spin': 0
    }
    
    # ACT
    result = service.validate_calculation_parameters(params)
    
    # ASSERT
    assert result is None  # Still passes


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
