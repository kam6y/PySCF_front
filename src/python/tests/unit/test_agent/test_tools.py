"""
Unit tests for AI Agent tools.

Tests the agent tool wrapper functions with mocked services to verify
parameter validation, error handling, and JSON response formatting.
"""

import pytest
import json
from unittest.mock import MagicMock, patch

from agent.computational_chemist import tools
from services.exceptions import ServiceError, NotFoundError, ValidationError


# ============================================================================
# Helper Functions Tests
# ============================================================================

def test_validation_error_helper():
    """
    GIVEN an error message
    WHEN _validation_error is called
    THEN it should return formatted JSON error
    """
    # ACT
    result = tools._validation_error("Test validation error")
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert data['error'] == "Test validation error"


def test_handle_service_error_helper():
    """
    GIVEN a ServiceError with message
    WHEN _handle_service_error is called
    THEN it should return formatted JSON error
    """
    # ARRANGE
    error = ServiceError("Test service error")
    
    # ACT
    result = tools._handle_service_error(error, "test_function")
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "test_function" in data['error']
    assert "Test service error" in data['error']


# ============================================================================
# list_all_calculations() Tests
# ============================================================================

@patch('agent.tools.get_quantum_service')
def test_list_all_calculations_success(mock_get_service):
    """
    GIVEN quantum service returns calculations list
    WHEN list_all_calculations is called
    THEN it should return formatted JSON with calculations
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.list_calculations.return_value = {
        'base_directory': '/data/calculations',
        'calculations': [
            {'id': 'calc1', 'name': 'Test 1'},
            {'id': 'calc2', 'name': 'Test 2'}
        ],
        'count': 2
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.list_all_calculations()
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert 'data' in data
    assert data['data']['count'] == 2


@patch('agent.tools.get_quantum_service')
def test_list_all_calculations_service_error(mock_get_service):
    """
    GIVEN quantum service raises ServiceError
    WHEN list_all_calculations is called
    THEN it should return error JSON
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.list_calculations.side_effect = ServiceError("Database error")
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.list_all_calculations()
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "error" in data


# ============================================================================
# get_calculation_details() Tests
# ============================================================================

@patch('agent.tools.get_quantum_service')
def test_get_calculation_details_success(mock_get_service):
    """
    GIVEN valid calculation ID and service returns details
    WHEN get_calculation_details is called
    THEN it should return formatted JSON with details
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.get_calculation_details.return_value = {
        'calculation': {
            'id': 'calc123',
            'name': 'Water DFT',
            'status': 'completed'
        }
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.get_calculation_details('calc123')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert data['data']['calculation']['id'] == 'calc123'


@patch('agent.tools.get_quantum_service')
def test_get_calculation_details_invalid_id(mock_get_service):
    """
    GIVEN invalid calculation ID (empty string)
    WHEN get_calculation_details is called
    THEN it should return validation error
    """
    # ACT
    result = tools.get_calculation_details('')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "required" in data['error'].lower()


@patch('agent.tools.get_quantum_service')
def test_get_calculation_details_not_found(mock_get_service):
    """
    GIVEN calculation ID that doesn't exist
    WHEN get_calculation_details is called
    THEN it should return error JSON
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.get_calculation_details.side_effect = NotFoundError("Calculation not found")
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.get_calculation_details('nonexistent')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False


# ============================================================================
# search_pubchem_by_name() Tests
# ============================================================================

@patch('agent.tools.get_pubchem_service')
def test_search_pubchem_by_name_success(mock_get_service):
    """
    GIVEN valid compound name and service returns data
    WHEN search_pubchem_by_name is called
    THEN it should return formatted JSON with compound info
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.search_compound.return_value = {
        'xyz': 'H 0 0 0',
        'compound_info': {'cid': 962, 'iupac_name': 'water'}
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.search_pubchem_by_name('water')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert data['data']['compound_info']['cid'] == 962


@patch('agent.tools.get_pubchem_service')
def test_search_pubchem_by_name_empty_name(mock_get_service):
    """
    GIVEN empty compound name
    WHEN search_pubchem_by_name is called
    THEN it should return validation error
    """
    # ACT
    result = tools.search_pubchem_by_name('')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "required" in data['error'].lower()


@patch('agent.tools.get_pubchem_service')
def test_search_pubchem_by_name_invalid_search_type(mock_get_service):
    """
    GIVEN invalid search_type parameter
    WHEN search_pubchem_by_name is called
    THEN it should return validation error
    """
    # ACT
    result = tools.search_pubchem_by_name('water', search_type='invalid')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "search_type" in data['error'].lower()


# ============================================================================
# convert_smiles_to_xyz() Tests
# ============================================================================

@patch('agent.tools.get_smiles_service')
def test_convert_smiles_to_xyz_success(mock_get_service):
    """
    GIVEN valid SMILES string and service returns XYZ
    WHEN convert_smiles_to_xyz is called
    THEN it should return formatted JSON with XYZ data
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.convert_smiles.return_value = {
        'xyz': '3\nEthanol\nC 0 0 0...'
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.convert_smiles_to_xyz('CCO')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert 'xyz' in data['data']


@patch('agent.tools.get_smiles_service')
def test_convert_smiles_to_xyz_empty_smiles(mock_get_service):
    """
    GIVEN empty SMILES string
    WHEN convert_smiles_to_xyz is called
    THEN it should return validation error
    """
    # ACT
    result = tools.convert_smiles_to_xyz('')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "required" in data['error'].lower()


@patch('agent.tools.get_smiles_service')
def test_convert_smiles_to_xyz_too_long(mock_get_service):
    """
    GIVEN SMILES string exceeding max length
    WHEN convert_smiles_to_xyz is called
    THEN it should return validation error
    """
    # ACT
    long_smiles = 'C' * 501
    result = tools.convert_smiles_to_xyz(long_smiles)
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "too long" in data['error'].lower()


# ============================================================================
# start_quantum_calculation() Tests
# ============================================================================

@patch('agent.tools.get_quantum_service')
def test_start_quantum_calculation_valid_params(mock_get_service):
    """
    GIVEN valid calculation parameters
    WHEN start_quantum_calculation is called
    THEN it should return formatted JSON with calculation instance
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.start_calculation.return_value = {
        'id': 'calc123',
        'name': 'Test Calculation',
        'status': 'pending'
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.start_quantum_calculation(
        xyz='H 0 0 0\nH 0 0 0.74',
        calculation_method='HF',
        basis_function='sto-3g',
        charges=0,
        spin=0
    )
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert data['data']['calculation']['id'] == 'calc123'


@patch('agent.tools.get_quantum_service')
def test_start_quantum_calculation_empty_xyz(mock_get_service):
    """
    GIVEN empty XYZ string
    WHEN start_quantum_calculation is called
    THEN it should return validation error
    """
    # ACT
    result = tools.start_quantum_calculation(
        xyz='',
        calculation_method='HF',
        basis_function='sto-3g',
        charges=0,
        spin=0
    )
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "xyz" in data['error'].lower()


@patch('agent.tools.get_quantum_service')
def test_start_quantum_calculation_invalid_method(mock_get_service):
    """
    GIVEN invalid calculation_method
    WHEN start_quantum_calculation is called
    THEN it should return validation error
    """
    # ACT
    result = tools.start_quantum_calculation(
        xyz='H 0 0 0',
        calculation_method='INVALID',
        basis_function='sto-3g',
        charges=0,
        spin=0
    )
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "calculation_method" in data['error'].lower()


@patch('agent.tools.get_quantum_service')
def test_start_quantum_calculation_invalid_charges(mock_get_service):
    """
    GIVEN charges outside valid range
    WHEN start_quantum_calculation is called
    THEN it should return validation error
    """
    # ACT
    result = tools.start_quantum_calculation(
        xyz='H 0 0 0',
        calculation_method='HF',
        basis_function='sto-3g',
        charges=100,  # Outside -10 to 10 range
        spin=0
    )
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "charges" in data['error'].lower()


@patch('agent.tools.get_quantum_service')
def test_start_quantum_calculation_dft_requires_xc(mock_get_service):
    """
    GIVEN DFT method with exchange_correlation
    WHEN start_quantum_calculation is called
    THEN it should include XC functional in parameters
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.start_calculation.return_value = {
        'id': 'calc123',
        'name': 'DFT Test',
        'status': 'pending'
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.start_quantum_calculation(
        xyz='H 0 0 0',
        calculation_method='DFT',
        basis_function='6-31G',
        exchange_correlation='B3LYP',
        charges=0,
        spin=0
    )
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    
    # Verify the service was called with XC functional
    call_args = mock_service.start_calculation.call_args[0][0]
    assert call_args['exchange_correlation'] == 'B3LYP'


# ============================================================================
# get_supported_parameters() Tests
# ============================================================================

@patch('agent.tools.get_quantum_service')
def test_get_supported_parameters_success(mock_get_service):
    """
    GIVEN quantum service returns supported parameters
    WHEN get_supported_parameters is called
    THEN it should return formatted JSON with parameters
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.get_supported_parameters.return_value = {
        'calculation_methods': ['HF', 'DFT'],
        'basis_sets': ['sto-3g', '6-31G']
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.get_supported_parameters()
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert 'calculation_methods' in data['data']


# ============================================================================
# delete_calculation() Tests
# ============================================================================

@patch('agent.tools.get_quantum_service')
def test_delete_calculation_returns_confirmation_request(mock_get_service):
    """
    GIVEN valid calculation ID
    WHEN delete_calculation is called
    THEN it should return confirmation request (not actually delete)
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.get_calculation_details.return_value = {
        'calculation': {
            'id': 'calc123',
            'name': 'Test Calculation'
        }
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.delete_calculation('calc123')
    
    # ASSERT
    data = json.loads(result)
    assert 'requires_confirmation' in data
    assert data['requires_confirmation'] is True
    assert data['action'] == 'delete_calculation'
    assert data['calculation_id'] == 'calc123'


@patch('agent.tools.get_quantum_service')
def test_delete_calculation_empty_id(mock_get_service):
    """
    GIVEN empty calculation ID
    WHEN delete_calculation is called
    THEN it should return validation error
    """
    # ACT
    result = tools.delete_calculation('')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "required" in data['error'].lower()


# ============================================================================
# validate_xyz_format() Tests
# ============================================================================

@patch('agent.tools.get_pubchem_service')
def test_validate_xyz_format_valid(mock_get_service):
    """
    GIVEN valid XYZ format string
    WHEN validate_xyz_format is called
    THEN it should return validation result with valid=true
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.validate_xyz.return_value = {
        'valid': True,
        'atom_count': 2
    }
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.validate_xyz_format('H 0 0 0\nH 0 0 0.74')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is True
    assert data['data']['valid'] is True


@patch('agent.tools.get_pubchem_service')
def test_validate_xyz_format_empty(mock_get_service):
    """
    GIVEN empty XYZ string
    WHEN validate_xyz_format is called
    THEN it should return validation error
    """
    # ACT
    result = tools.validate_xyz_format('')
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert "required" in data['error'].lower()


# ============================================================================
# JSON Response Format Tests
# ============================================================================

@pytest.mark.parametrize("tool_func,args", [
    (tools.list_all_calculations, []),
    (tools.get_supported_parameters, []),
])
@patch('agent.tools.get_quantum_service')
def test_tools_return_valid_json(mock_get_service, tool_func, args):
    """
    GIVEN any tool function
    WHEN the tool is called
    THEN it should return valid JSON string
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.list_calculations.return_value = {'count': 0}
    mock_service.get_supported_parameters.return_value = {}
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tool_func(*args)
    
    # ASSERT
    # Should be valid JSON
    data = json.loads(result)
    assert isinstance(data, dict)
    assert 'success' in data


@patch('agent.tools.get_quantum_service')
def test_error_responses_have_consistent_format(mock_get_service):
    """
    GIVEN tool encounters an error
    WHEN error response is returned
    THEN it should have consistent JSON format
    """
    # ARRANGE
    mock_service = MagicMock()
    mock_service.list_calculations.side_effect = ServiceError("Test error")
    mock_get_service.return_value = mock_service
    
    # ACT
    result = tools.list_all_calculations()
    
    # ASSERT
    data = json.loads(result)
    assert data['success'] is False
    assert 'error' in data
    assert isinstance(data['error'], str)
