"""
Unit tests for SMILES service.

Tests the SMILESService class with mocked external dependencies to verify
SMILES to XYZ conversion logic, validation, and error handling.
"""

import pytest

from services.smiles_service import SMILESService
from services.exceptions import ServiceError, ValidationError
from SMILES.smiles_converter import SMILESError


# ============================================================================
# convert_smiles() Success Tests
# ============================================================================

def test_convert_smiles_success(mocker):
    """
    GIVEN a valid SMILES string
    WHEN convert_smiles is called
    THEN it should return XYZ data
    """
    # ARRANGE
    smiles = 'CCO'  # Ethanol
    expected_xyz = "9\nMolecule from SMILES: CCO\nC 0.0 0.0 0.0\nC 1.5 0.0 0.0\nO 2.0 1.4 0.0\n..."
    
    mocker.patch('services.smiles_service.smiles_to_xyz', return_value=expected_xyz)
    
    service = SMILESService()
    
    # ACT
    result = service.convert_smiles(smiles)
    
    # ASSERT
    assert result is not None
    assert 'xyz' in result
    assert result['xyz'] == expected_xyz


def test_convert_smiles_default_title(mocker):
    """
    GIVEN a valid SMILES string without title
    WHEN convert_smiles is called
    THEN it should use default title format
    """
    # ARRANGE
    smiles = 'CC(=O)O'  # Acetic acid
    
    mock_converter = mocker.patch('services.smiles_service.smiles_to_xyz', return_value='xyz_data')
    
    service = SMILESService()
    
    # ACT
    service.convert_smiles(smiles)
    
    # ASSERT
    # Verify the converter was called with default title
    mock_converter.assert_called_once_with(smiles, title=f"Molecule from SMILES: {smiles}")

# ============================================================================
# convert_smiles() Validation Error Tests
# ============================================================================

def test_convert_smiles_empty_string(mocker):
    """
    GIVEN an empty SMILES string
    WHEN convert_smiles is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="SMILES string is required"):
        service.convert_smiles("")


def test_convert_smiles_whitespace_only(mocker):
    """
    GIVEN a whitespace-only SMILES string
    WHEN convert_smiles is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="cannot be empty or contain only whitespace"):
        service.convert_smiles("   \t\n   ")


def test_convert_smiles_none_input(mocker):
    """
    GIVEN None as SMILES input
    WHEN convert_smiles is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="SMILES string is required"):
        service.convert_smiles(None)


def test_convert_smiles_non_string_input(mocker):
    """
    GIVEN a non-string SMILES input
    WHEN convert_smiles is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="must be a non-empty string"):
        service.convert_smiles(12345)


def test_convert_smiles_too_long(mocker):
    """
    GIVEN a SMILES string exceeding maximum length
    WHEN convert_smiles is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    service = SMILESService()
    long_smiles = 'C' * 501  # Exceeds 500 character limit
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="too long.*Maximum length is 500"):
        service.convert_smiles(long_smiles)


# ============================================================================
# convert_smiles() SMILES Error Tests
# ============================================================================

def test_convert_smiles_invalid_structure(mocker):
    """
    GIVEN smiles_to_xyz raises SMILESError for invalid SMILES
    WHEN convert_smiles is called
    THEN it should raise ValidationError
    """
    # ARRANGE
    mocker.patch(
        'services.smiles_service.smiles_to_xyz',
        side_effect=SMILESError('Invalid SMILES string')
    )
    
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ValidationError, match="Invalid SMILES string"):
        service.convert_smiles('INVALID_SMILES_XYZ')


# ============================================================================
# convert_smiles() Unexpected Error Tests
# ============================================================================

def test_convert_smiles_unexpected_error(mocker):
    """
    GIVEN smiles_to_xyz raises an unexpected exception
    WHEN convert_smiles is called
    THEN it should raise ServiceError with generic message
    """
    # ARRANGE
    mocker.patch(
        'services.smiles_service.smiles_to_xyz',
        side_effect=RuntimeError('Unexpected error')
    )
    
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError, match="An internal server error occurred"):
        service.convert_smiles('CCO')


def test_convert_smiles_memory_error(mocker):
    """
    GIVEN smiles_to_xyz raises MemoryError
    WHEN convert_smiles is called
    THEN it should raise ServiceError
    """
    # ARRANGE
    mocker.patch(
        'services.smiles_service.smiles_to_xyz',
        side_effect=MemoryError('Out of memory')
    )
    
    service = SMILESService()
    
    # ACT & ASSERT
    with pytest.raises(ServiceError, match="An internal server error occurred"):
        service.convert_smiles('C' * 100)


# ============================================================================
# Edge Case Tests
# ============================================================================

def test_convert_smiles_strips_whitespace(mocker):
    """
    GIVEN a SMILES string with leading/trailing whitespace
    WHEN convert_smiles is called
    THEN it should strip whitespace before conversion
    """
    # ARRANGE
    smiles_with_whitespace = '  CCO  \t'
    expected_smiles = 'CCO'
    
    mock_converter = mocker.patch('services.smiles_service.smiles_to_xyz', return_value='xyz_data')
    
    service = SMILESService()
    
    # ACT
    service.convert_smiles(smiles_with_whitespace)
    
    # ASSERT
    # Verify that the stripped SMILES was passed to the converter
    # The first argument should be the stripped version
    call_args = mock_converter.call_args[0]
    assert call_args[0] == expected_smiles


def test_convert_smiles_preserves_internal_structure(mocker):
    """
    GIVEN a complex SMILES string with various features
    WHEN convert_smiles is called
    THEN it should preserve the structure
    """
    # ARRANGE
    complex_smiles = 'CC(C)C(=O)O'  # Isobutyric acid
    
    mock_converter = mocker.patch('services.smiles_service.smiles_to_xyz', return_value='xyz_data')
    
    service = SMILESService()
    
    # ACT
    service.convert_smiles(complex_smiles)
    
    # ASSERT
    mock_converter.assert_called_once()
    call_args = mock_converter.call_args[0]
    assert call_args[0] == complex_smiles


# ============================================================================
# Security Regression Tests — Input Redaction
# ============================================================================

def test_convert_smiles_invalid_error_does_not_reflect_input():
    """
    GIVEN an invalid SMILES string containing a proprietary sentinel
    WHEN convert_smiles is called (real converter, no mocks)
    THEN the raised ValidationError must NOT contain the sentinel value
    """
    # ARRANGE
    sentinel = 'PROPRIETARY_MOLECULE_%%%_SECRET'
    service = SMILESService()

    # ACT & ASSERT
    with pytest.raises(ValidationError) as exc_info:
        service.convert_smiles(sentinel)

    error_message = str(exc_info.value)
    assert sentinel not in error_message, (
        f"Raw SMILES input was reflected in service error: {error_message}"
    )
