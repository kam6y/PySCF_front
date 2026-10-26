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
    PubChemService(timeout=60)
    
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
    PubChemService()
    
    # ASSERT
    mock_client_class.assert_called_once_with(timeout=30)


# ============================================================================
# Security Regression Tests — Input Redaction
# ============================================================================

def test_search_compound_not_found_error_does_not_reflect_query(mocker):
    """
    GIVEN PubChemClient returns compound data without atoms (no 3D structure)
    WHEN search_compound is called with a proprietary query sentinel
    THEN the raised NotFoundError must NOT contain the sentinel value
    """
    # ARRANGE
    sentinel = 'SECRET_COMPOUND_NAME_%%%_PRIVATE'
    mock_compound_data = MagicMock()
    mock_compound_data.cid = 99999
    mock_compound_data.atoms = []  # No 3D structure

    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        return_value=mock_compound_data
    )

    service = PubChemService()

    # ACT & ASSERT
    with pytest.raises(NotFoundError) as exc_info:
        service.search_compound(sentinel, 'name')

    error_message = str(exc_info.value)
    assert sentinel not in error_message, (
        f"Raw query was reflected in not-found error: {error_message}"
    )


def test_search_compound_not_found_error_is_generic(mocker):
    """
    GIVEN PubChemClient returns compound data without atoms
    WHEN search_compound is called
    THEN the error message should be the generic 'No compound with a 3D structure found'
    """
    # ARRANGE
    mock_compound_data = MagicMock()
    mock_compound_data.cid = 99999
    mock_compound_data.atoms = []

    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        return_value=mock_compound_data
    )

    service = PubChemService()

    # ACT & ASSERT
    with pytest.raises(NotFoundError, match="No compound with a 3D structure found"):
        service.search_compound('anything', 'name')


def test_search_compound_pubchem_not_found_does_not_reflect_query(mocker):
    """
    GIVEN PubChemClient raises PubChemNotFoundError
    WHEN search_compound is called with a proprietary query sentinel
    THEN the raised NotFoundError must NOT contain the sentinel value
    """
    # ARRANGE
    sentinel = 'TOP_SECRET_FORMULA_%%%'
    mocker.patch(
        'services.pubchem_service.PubChemClient.search_compound',
        side_effect=PubChemNotFoundError(
            'A 3D structure is not available for the compound with CID 12345.'
        )
    )

    service = PubChemService()

    # ACT & ASSERT
    with pytest.raises(NotFoundError) as exc_info:
        service.search_compound(sentinel, 'name')

    error_message = str(exc_info.value)
    assert sentinel not in error_message, (
        f"Raw query was reflected in not-found error: {error_message}"
    )


# ============================================================================
# Security Regression Tests — Exception Chain Redaction
# ============================================================================

def test_network_error_chain_does_not_leak_query_via_exc_info(mocker, caplog):
    """
    GIVEN a name-search where a network error embeds the raw query in its message
    WHEN the error propagates through client and service, logged with exc_info=True
    THEN the sentinel must NOT appear anywhere in caplog (including chained tracebacks)

    This exercises the real client._make_request -> _find_cid -> search_compound
    -> pubchem_service exc_info=True chain. The from None on client line 54
    suppresses the chained ConnectionError that carries the URL with the query.
    """
    import requests.exceptions
    import logging

    sentinel = 'CHAIN_LEAK_SENTINEL_NET_%%%'
    fake_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{sentinel}/cids/JSON'
    simulated_err = requests.exceptions.ConnectionError(
        f'HTTPConnectionPool: Max retries exceeded with url: {fake_url}'
    )

    service = PubChemService()
    # Mock at session level so real _make_request except-block runs
    mocker.patch.object(service.client.session, 'get', side_effect=simulated_err)

    with caplog.at_level(logging.DEBUG, logger='services.pubchem_service'):
        with pytest.raises(ServiceError):
            service.search_compound(sentinel, 'name')

    # caplog.text includes exc_info formatted tracebacks
    assert sentinel not in caplog.text, (
        f'Sentinel leaked via exception chain in exc_info log'
    )


def test_processing_error_chain_does_not_leak_via_exc_info(mocker, caplog):
    """
    GIVEN a downstream ValueError containing a sentinel occurs during compound processing
    WHEN it is caught by search_compound's (ValueError, KeyError, TypeError) handler
    AND the resulting PubChemError is logged with exc_info=True at the service layer
    THEN the sentinel must NOT appear in caplog

    This exercises client.py lines 75-77 (the from None on the re-raise)
    and line 76 (the sanitized log that no longer includes {e}).
    """
    import logging

    sentinel = 'CHAIN_LEAK_SENTINEL_PROC_%%%'

    service = PubChemService()
    # Mock _find_cid to return a valid CID
    mocker.patch.object(service.client, '_find_cid', return_value=12345)
    # Mock _get_compound_properties to raise ValueError with sentinel
    mocker.patch.object(
        service.client, '_get_compound_properties',
        side_effect=ValueError(f'unexpected value: {sentinel}')
    )

    with caplog.at_level(logging.DEBUG, logger='services.pubchem_service'):
        with pytest.raises(ServiceError):
            service.search_compound(sentinel, 'name')

    assert sentinel not in caplog.text, (
        f'Sentinel leaked via exception chain in exc_info log'
    )
