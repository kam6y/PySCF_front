"""Tests for PubChem client functionality, rewritten to match the current implementation."""

import pytest
from unittest.mock import patch, Mock
import requests

# We assume the package is installed in editable mode, so no sys.path modification is needed.

from pubchem.client import PubChemClient, PubChemError, CompoundData
from pubchem.parser import atoms_to_xyz

# --- Mock Data ---
MOCK_WATER_CID = 962
MOCK_WATER_PROPERTIES_JSON = {
    "PropertyTable": {
        "Properties": [
            {
                "CID": MOCK_WATER_CID,
                "MolecularFormula": "H2O",
                "MolecularWeight": "18.015",
                "IUPACName": "water"
            }
        ]
    }
}
MOCK_WATER_SYNONYMS_JSON = {
    "InformationList": {
        "Information": [
            {"CID": MOCK_WATER_CID, "Synonym": ["Water", "H2O", "Dihydrogen oxide"]}
        ]
    }
}
MOCK_WATER_SDF = f"""
Molecule from PubChem

  -OEChem-11111111113D

  3  2  0     0  0  0  0  0  0999 V2000
    -0.0000    0.0000    0.1193 O   0  0  0  0  0  0  0  0  0  0  0  0
     0.0000    0.7632   -0.4770 H   0  0  0  0  0  0  0  0  0  0  0  0
     0.0000   -0.7632   -0.4770 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""
MOCK_NOT_FOUND_JSON = {"Fault": {"Code": "PUG-Client-Error.NotFound"}}
MOCK_CIDS_JSON = {"IdentifierList": {"CID": [MOCK_WATER_CID]}}


# --- Helper for Mocking Requests ---
def mock_requests_get(*args, **kwargs):
    """A factory for creating mock responses for requests.get."""
    url = args[0]
    mock_response = Mock(spec=requests.Response)
    mock_response.status_code = 200
    mock_response.raise_for_status.return_value = None

    if f"/compound/name/water/cids/JSON" in url:
        mock_response.json.return_value = MOCK_CIDS_JSON
    elif f"/cid/{MOCK_WATER_CID}/property/" in url:
        mock_response.json.return_value = MOCK_WATER_PROPERTIES_JSON
    elif f"/cid/{MOCK_WATER_CID}/synonyms/" in url:
        mock_response.json.return_value = MOCK_WATER_SYNONYMS_JSON
    elif f"/cid/{MOCK_WATER_CID}/SDF" in url:
        mock_response.text = MOCK_WATER_SDF
    elif "/name/notfound/cids/JSON" in url:
        mock_response.status_code = 404
        mock_response.json.return_value = MOCK_NOT_FOUND_JSON
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError("Not Found")
    else:
        # Default 404 for any other unhandled URL
        mock_response.status_code = 404
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError("Mock Not Found")
        
    return mock_response


@pytest.fixture
def client():
    """Fixture for PubChemClient."""
    return PubChemClient(timeout=10)


@patch('requests.Session.get', side_effect=mock_requests_get)
def test_search_compound_by_name_success(mock_get, client):
    """Test successful compound search by name."""
    result = client.search_compound("water", "name")
    
    assert isinstance(result, CompoundData)
    assert result.cid == MOCK_WATER_CID
    assert result.iupac_name == "water"
    assert result.molecular_formula == "H2O"
    assert len(result.atoms) == 3
    assert result.atoms[0]['element'] == 'O'
    assert "Dihydrogen oxide" in result.synonyms
    assert mock_get.call_count == 4 # cids, properties, synonyms, sdf

def test_search_compound_by_cid(client):
    """Test search by a valid CID string."""
    # No mock needed for this part of the logic
    cid = client._find_cid("962", "cid")
    assert cid == 962

def test_search_compound_invalid_cid(client):
    """Test compound search with an invalid CID format."""
    with pytest.raises(PubChemError, match="Invalid CID format"):
        client.search_compound("invalid", "cid")

@patch('requests.Session.get', side_effect=mock_requests_get)
def test_search_compound_not_found(mock_get, client):
    """Test compound search where the compound is not found on PubChem."""
    result = client.search_compound("notfound", "name")
    assert result is None

@patch('requests.Session.get')
def test_get_3d_structure_from_sdf_unavailable(mock_get, client):
    """Test case where properties are found but 3D SDF is not."""
    # Simulate 404 for the SDF request
    def custom_mock(*args, **kwargs):
        url = args[0]
        if "SDF" in url:
            raise PubChemError("404 Not Found") # Simulate failure from _make_request
        return mock_requests_get(*args, **kwargs)

    mock_get.side_effect = custom_mock
    
    result = client.search_compound("water", "name")
    # The search should return None because the 3D structure is missing
    assert result is None


@patch('requests.Session.get')
def test_api_request_failure(mock_get, client):
    """Test a general request failure (e.g., network error)."""
    mock_get.side_effect = requests.exceptions.Timeout("Connection timed out")
    
    with pytest.raises(PubChemError, match="API request failed"):
        client.search_compound("water", "name")

def test_parser_integration(client):
    """Test that the parsed atoms can be converted to XYZ."""
    atoms = [
        {'element': 'C', 'x': 0.0, 'y': 0.0, 'z': 0.0},
        {'element': 'H', 'x': 1.0, 'y': 0.0, 'z': 0.0},
    ]
    xyz_str = atoms_to_xyz(atoms, "Methane")
    assert "2" in xyz_str
    assert "Methane" in xyz_str
    assert "C      0.000000" in xyz_str