"""
Regression tests for SMILES converter input redaction.

Ensures that invalid SMILES input values are NOT reflected in error messages,
preventing potential leakage of proprietary molecular structures via logs,
error strings, or crash reports.
"""

import pytest

from SMILES.smiles_converter import smiles_to_xyz, SMILESError


# A sentinel value that should never appear in any error message.
PROPRIETARY_SMILES_SENTINEL = "PROPRIETARY_MOLECULE_%%%_SECRET"


def test_invalid_smiles_error_does_not_reflect_input():
    """
    GIVEN an invalid SMILES string containing a proprietary sentinel
    WHEN smiles_to_xyz is called (real RDKit path, no mocks)
    THEN the raised SMILESError must NOT contain the sentinel value
    """
    with pytest.raises(SMILESError) as exc_info:
        smiles_to_xyz(PROPRIETARY_SMILES_SENTINEL)

    error_message = str(exc_info.value)
    assert (
        PROPRIETARY_SMILES_SENTINEL not in error_message
    ), f"Raw SMILES input was reflected in error message: {error_message}"


def test_invalid_smiles_error_is_generic():
    """
    GIVEN an invalid SMILES string
    WHEN smiles_to_xyz is called
    THEN the error message should be the generic 'Invalid SMILES string'
    """
    with pytest.raises(SMILESError, match="Invalid SMILES string"):
        smiles_to_xyz("NOT_A_VALID_SMILES_@@@")
