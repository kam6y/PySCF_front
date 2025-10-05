"""
Unit tests for PubChem service.

Tests the PubChemService class with mocked external dependencies to verify
business logic, error handling, and data transformation without actual API calls.
"""

import pytest
from unittest.mock import MagicMock

from services.pubchem_service import PubChemService
from services.exceptions import ServiceError, NotFoundError, ValidationError
from pubchem.client import PubChemError, PubChemNotFoundError


# ============================================================================
# search_compound() Tests
# ============================================================================

def test_search_compound_by_name_success(mocker):
    """
    GIVEN PubChemClient returns valid compound data
    WHEN search_compound is called with a compound name
    THEN it should return formatted XYZ data with compound information
    """
    # ARRANGE
    # Create mock compound data
    mock_compound_data = MagicMock()
    mock_compound_data.cid = 962
    mock_compound_data.iupac_name = 'water'
    mock_compound_data.molecular_formula = 'H2O'
    mock_compound_data.molecular_weight = 18.015
    mock_compound_data.synonyms = ['water', 'oxidane', 'dihydrogen oxide']
    mock_compound_data.atoms = [
        ['O', 0.0000, 0.0000, 0.1173],
        ['H', 0.0000, 0.7572, -0.4692],
        ['H', 0.0000, -0.7572, -0.4692]
    ]
    
    # Mock the PubChemClient.search_compound method
    mocker.patch('services.pubchem_service.PubChemClient.search_compound', return_value=mock_compound_data)
    
    # Mock the parser functions
    mocker.patch('services.pubchem_service.xyz_parser.format_compound_title', return_value='Water (CID: 962)')
    mocker.patch('services.pubchem_service.xyz_parser.atoms_to_xyz', return_value='3\nWater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692')
    
    service = PubChemService()
    
    # ACT
    result = service.search_compound('water', 'name')
    
    # ASSERT
    assert result is not None
    assert 'xyz' in result
    assert 'compound_info' in result
    assert 'atom_count' in result
    
    assert result['compound_info']['cid'] == 962
    assert result['compound_info']['iupac_name'] == 'water'
    assert result['compound_info']['molecular_formula'] == 'H2O'
    assert result['atom_count'] == 3


def test_search_compound_by_cid_success(mocker):
    """
    GIVEN PubChemClient returns valid compound data for a CID
    WHEN search_compound is called with search_type='cid'
    THEN it should return formatted XYZ data
    """
    # ARRANGE
    mock_compound_data = MagicMock()
    mock_compound_data.cid = 241
    mock_compound_data.iupac_name = 'benzene'
    mock_compound_data.molecular_formula = 'C6H6'
    mock_compound_data.molecular_weight = 78.114
    mock_compound_data.synonyms = ['benzene']
    mock_compound_data.atoms = [['C', 0, 0, 0]] * 6 + [['H', 0, 0, 0]] * 6
    
    mocker.patch('services.pubchem_service.PubChemClient.search_compound', return_value=mock_compound_data)
    mocker.patch('services.pubchem_service.xyz_parser.format_compound_title', return_value='Benzene')
    mocker.patch('services.pubchem_service.xyz_parser.atoms_to_xyz', return_value='12\nBenzene\n...')
    
    service = PubChemService()
    
    # ACT
    result = service.search_compound('241', 'cid')
    
    # ASSERT
    assert result['compound_info']['cid'] == 241
    assert result['atom_count'] == 12


def test_search_compound_invalid_search_type(mocker):
    """
    GIVEN an invalid search_type parameter
    WHEN search_compound is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="Invalid search type"):
        service.search_compound('water', 'invalid_type')


def test_search_compound_not_found(mocker):
    """
    GIVEN PubChemClient raises PubChemNotFoundError
    WHEN search_compound is called
    THEN it should raise NotFoundError
    """
    # ARRANGE
    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        side_effect=PubChemNotFoundError('Compound not found')
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(NotFoundError, match="Compound not found"):
        service.search_compound('nonexistent_compound_xyz123', 'name')


def test_search_compound_no_3d_structure(mocker):
    """
    GIVEN PubChemClient returns compound data without atoms (no 3D structure)
    WHEN search_compound is called
    THEN it should raise NotFoundError
    """
    # ARRANGE
    mock_compound_data = MagicMock()
    mock_compound_data.cid = 123
    mock_compound_data.atoms = []  # No 3D structure
    
    mocker.patch('services.pubchem_service.PubChemClient.search_compound', return_value=mock_compound_data)
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(NotFoundError, match="No compound with a 3D structure found"):
        service.search_compound('some_compound', 'name')


def test_search_compound_api_error_404(mocker):
    """
    GIVEN PubChemClient raises PubChemError with 404 status
    WHEN search_compound is called
    THEN it should raise NotFoundError
    """
    # ARRANGE
    error = PubChemError('Not found')
    error.status_code = 404
    
    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        side_effect=error
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(NotFoundError):
        service.search_compound('unknown', 'name')


def test_search_compound_api_error_400(mocker):
    """
    GIVEN PubChemClient raises PubChemError with 400 status
    WHEN search_compound is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    error = PubChemError('Bad request')
    error.status_code = 400
    
    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        side_effect=error
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError):
        service.search_compound('invalid', 'name')


def test_search_compound_api_error_generic(mocker):
    """
    GIVEN PubChemClient raises PubChemError without specific status code
    WHEN search_compound is called
    THEN it should raise ServiceError
    """
    # ARRANGE
    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        side_effect=PubChemError('Connection failed')
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError):
        service.search_compound('water', 'name')


def test_search_compound_unexpected_error(mocker):
    """
    GIVEN PubChemClient raises an unexpected exception
    WHEN search_compound is called
    THEN it should raise ServiceError with generic message
    """
    # ARRANGE
    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        side_effect=RuntimeError('Unexpected error')
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError, match="An internal server error occurred"):
        service.search_compound('water', 'name')


# ============================================================================
# validate_xyz() Tests
# ============================================================================

def test_validate_xyz_success(mocker):
    """
    GIVEN a valid XYZ format string
    WHEN validate_xyz is called
    THEN it should return validation results with valid=True
    """
    # ARRANGE
    xyz_string = "3\nWater\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"
    
    mock_validation_result = {
        'valid': True,
        'atom_count': 3,
        'atoms': [
            {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.1173},
            {'element': 'H', 'x': 0.0, 'y': 0.7572, 'z': -0.4692},
            {'element': 'H', 'x': 0.0, 'y': -0.7572, 'z': -0.4692}
        ]
    }
    
    mocker.patch('services.pubchem_service.xyz_parser.validate_xyz', return_value=mock_validation_result)
    
    service = PubChemService()
    
    # ACT
    result = service.validate_xyz(xyz_string)
    
    # ASSERT
    assert result is not None
    assert result['valid'] is True
    assert result['atom_count'] == 3
    assert len(result['atoms']) == 3


def test_validate_xyz_invalid_format(mocker):
    """
    GIVEN an invalid XYZ format string
    WHEN validate_xyz is called
    THEN it should return validation results with valid=False
    """
    # ARRANGE
    xyz_string = "invalid xyz format"
    
    mock_validation_result = {
        'valid': False,
        'error': 'Invalid XYZ format'
    }
    
    mocker.patch('services.pubchem_service.xyz_parser.validate_xyz', return_value=mock_validation_result)
    
    service = PubChemService()
    
    # ACT
    result = service.validate_xyz(xyz_string)
    
    # ASSERT
    assert result is not None
    assert result['valid'] is False
    assert 'error' in result


def test_validate_xyz_empty_string(mocker):
    """
    GIVEN an empty XYZ string
    WHEN validate_xyz is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="XYZ string cannot be empty"):
        service.validate_xyz("")


def test_validate_xyz_whitespace_only(mocker):
    """
    GIVEN a whitespace-only XYZ string
    WHEN validate_xyz is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="XYZ string cannot be empty"):
        service.validate_xyz("   \n  \t  ")


def test_validate_xyz_none_input(mocker):
    """
    GIVEN None as input
    WHEN validate_xyz is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="XYZ string cannot be empty"):
        service.validate_xyz(None)


def test_validate_xyz_memory_error(mocker):
    """
    GIVEN xyz_parser.validate_xyz raises MemoryError
    WHEN validate_xyz is called
    THEN it should raise ServiceError with status 413
    """
    # ARRANGE
    mocker.patch(
        'services.pubchem_service.xyz_parser.validate_xyz',
        side_effect=MemoryError('Too large')
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError, match="XYZ string too large"):
        service.validate_xyz("very large xyz data")


def test_validate_xyz_unexpected_error(mocker):
    """
    GIVEN xyz_parser.validate_xyz raises an unexpected exception
    WHEN validate_xyz is called
    THEN it should raise ServiceError
    """
    # ARRANGE
    mocker.patch(
        'services.pubchem_service.xyz_parser.validate_xyz',
        side_effect=RuntimeError('Unexpected error')
    )
    
    service = PubChemService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError, match="An internal server error occurred"):
        service.validate_xyz("H 0 0 0")


# ============================================================================
# Custom Timeout Configuration Tests
# ============================================================================

def test_pubchem_service_custom_timeout(mocker):
    """
    GIVEN a custom timeout parameter
    WHEN PubChemService is initialized
    THEN it should pass the timeout to PubChemClient
    """
    # ARRANGE
    mock_client_class = mocker.patch('services.pubchem_service.PubChemClient')
    
    # ACT
    service = PubChemService(timeout=60)
    
    # ASSERT
    mock_client_class.assert_called_once_with(timeout=60)


def test_pubchem_service_default_timeout(mocker):
    """
    GIVEN no timeout parameter
    WHEN PubChemService is initialized
    THEN it should use the default timeout of 30 seconds
    """
    # ARRANGE
    mock_client_class = mocker.patch('services.pubchem_service.PubChemClient')
    
    # ACT
    service = PubChemService()
    
    # ASSERT
    mock_client_class.assert_called_once_with(timeout=30)
