"""Tests for PubChem client functionality."""

import pytest
from unittest.mock import Mock, patch
import pubchempy as pcp

# Add the parent directory to the path for imports
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pubchem.client import PubChemClient, PubChemError
from pubchem.parser import XYZParser


class TestPubChemClient:
    """Test cases for PubChemClient."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.client = PubChemClient(timeout=10)
        self.parser = XYZParser()
    
    def test_init(self):
        """Test client initialization."""
        client = PubChemClient(timeout=20)
        assert client.timeout == 20
    
    @patch('pubchempy.get_compounds')
    def test_search_compound_by_name_success(self, mock_get_compounds):
        """Test successful compound search by name."""
        # Mock compound object
        mock_compound = Mock()
        mock_compound.cid = 962
        mock_get_compounds.return_value = [mock_compound]
        
        result = self.client.search_compound("water", "name")
        
        assert result == mock_compound
        mock_get_compounds.assert_called_once_with("water", "name", record_type="3d")
    
    @patch('pubchempy.get_compounds')
    def test_search_compound_by_cid_success(self, mock_get_compounds):
        """Test successful compound search by CID."""
        mock_compound = Mock()
        mock_compound.cid = 962
        mock_get_compounds.return_value = [mock_compound]
        
        result = self.client.search_compound("962", "cid")
        
        assert result == mock_compound
        mock_get_compounds.assert_called_once_with(962, "cid", record_type="3d")
    
    def test_search_compound_invalid_cid(self):
        """Test compound search with invalid CID."""
        with pytest.raises(PubChemError) as exc_info:
            self.client.search_compound("invalid_cid", "cid")
        
        assert "Invalid CID format" in str(exc_info.value)
    
    @patch('pubchempy.get_compounds')
    def test_search_compound_not_found(self, mock_get_compounds):
        """Test compound search when no results found."""
        mock_get_compounds.return_value = []
        
        result = self.client.search_compound("nonexistent", "name")
        
        assert result is None
    
    @patch('pubchempy.get_compounds')
    def test_search_compound_fallback_to_2d(self, mock_get_compounds):
        """Test fallback to 2D when 3D not available."""
        mock_compound = Mock()
        mock_compound.cid = 123
        
        # First call returns empty, second call returns compound
        mock_get_compounds.side_effect = [[], [mock_compound]]
        
        result = self.client.search_compound("test", "name")
        
        assert result == mock_compound
        assert mock_get_compounds.call_count == 2
    
    def test_get_3d_structure_with_atoms(self):
        """Test getting 3D structure from compound with atoms."""
        # Mock compound with atoms
        mock_atom1 = Mock()
        mock_atom1.element = "O"
        mock_atom1.x = 0.0
        mock_atom1.y = 0.0
        mock_atom1.z = 0.0
        
        mock_atom2 = Mock()
        mock_atom2.element = "H"
        mock_atom2.x = 1.0
        mock_atom2.y = 0.0
        mock_atom2.z = 0.0
        
        mock_compound = Mock()
        mock_compound.atoms = [mock_atom1, mock_atom2]
        mock_compound.cid = 962
        
        result = self.client.get_3d_structure(mock_compound)
        
        expected = [
            {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.0},
            {'element': 'H', 'x': 1.0, 'y': 0.0, 'z': 0.0}
        ]
        assert result == expected
    
    def test_get_3d_structure_no_atoms(self):
        """Test getting 3D structure from compound without atoms."""
        mock_compound = Mock()
        mock_compound.atoms = []
        mock_compound.cid = 123
        
        # Mock the SDF retrieval method
        with patch.object(self.client, '_get_3d_from_sdf', return_value=None):
            result = self.client.get_3d_structure(mock_compound)
            assert result is None
    
    def test_get_compound_info(self):
        """Test getting compound information."""
        mock_compound = Mock()
        mock_compound.cid = 962
        mock_compound.molecular_formula = "H2O"
        mock_compound.molecular_weight = 18.015
        mock_compound.iupac_name = "water"
        mock_compound.synonyms = ["water", "dihydrogen monoxide", "H2O"]
        
        result = self.client.get_compound_info(mock_compound)
        
        expected = {
            'cid': 962,
            'molecular_formula': 'H2O',
            'molecular_weight': 18.015,
            'iupac_name': 'water',
            'synonyms': ['water', 'dihydrogen monoxide', 'H2O']
        }
        assert result == expected
    
    @patch('requests.get')
    def test_get_3d_from_sdf(self, mock_get):
        """Test getting 3D structure from SDF."""
        # Mock SDF response
        sdf_data = """

  -OEChem-01012412483D

  3  2  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1193 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7632   -0.4770 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7632   -0.4770 H   0  0  0  0  0  0  0  0  0  0  0  0
"""
        mock_response = Mock()
        mock_response.text = sdf_data
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        result = self.client._get_3d_from_sdf(962)
        
        assert result is not None
        assert len(result) == 3
        assert result[0]['element'] == 'O'
        assert result[1]['element'] == 'H'
        assert result[2]['element'] == 'H'


class TestXYZParser:
    """Test cases for XYZParser."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.parser = XYZParser()
    
    def test_atoms_to_xyz(self):
        """Test converting atoms to XYZ format."""
        atoms = [
            {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.1193},
            {'element': 'H', 'x': 0.0, 'y': 0.7632, 'z': -0.4770},
            {'element': 'H', 'x': 0.0, 'y': -0.7632, 'z': -0.4770}
        ]
        
        result = self.parser.atoms_to_xyz(atoms, "Water")
        
        lines = result.strip().split('\n')
        assert lines[0] == '3'
        assert lines[1] == 'Water'
        assert 'O ' in lines[2]
        assert 'H ' in lines[3]
        assert 'H ' in lines[4]
    
    def test_atoms_to_xyz_empty(self):
        """Test converting empty atoms list."""
        with pytest.raises(ValueError) as exc_info:
            self.parser.atoms_to_xyz([], "Empty")
        
        assert "No atoms provided" in str(exc_info.value)
    
    def test_atoms_to_xyz_invalid_atom(self):
        """Test converting atoms with missing data."""
        atoms = [
            {'element': 'O', 'x': 0.0, 'y': 0.0}  # Missing 'z'
        ]
        
        with pytest.raises(ValueError) as exc_info:
            self.parser.atoms_to_xyz(atoms, "Invalid")
        
        assert "missing keys" in str(exc_info.value)
    
    def test_validate_xyz_valid(self):
        """Test validating valid XYZ string."""
        xyz_string = """3
Water molecule
O    0.000000    0.000000    0.119262
H    0.000000    0.763239   -0.477047
H    0.000000   -0.763239   -0.477047"""
        
        result = self.parser.validate_xyz(xyz_string)
        
        assert result['valid'] is True
        assert result['num_atoms'] == 3
        assert len(result['atoms']) == 3
        assert result['atoms'][0]['element'] == 'O'
    
    def test_validate_xyz_invalid_count(self):
        """Test validating XYZ with invalid atom count."""
        xyz_string = """invalid
Water molecule
O    0.000000    0.000000    0.119262"""
        
        result = self.parser.validate_xyz(xyz_string)
        
        assert result['valid'] is False
        assert "First line must be number of atoms" in result['error']
    
    def test_validate_xyz_too_short(self):
        """Test validating XYZ string that is too short."""
        xyz_string = "3"
        
        result = self.parser.validate_xyz(xyz_string)
        
        assert result['valid'] is False
        assert "XYZ string too short" in result['error']
    
    def test_format_compound_title(self):
        """Test formatting compound title."""
        compound_info = {
            'cid': 962,
            'molecular_formula': 'H2O',
            'iupac_name': 'water'
        }
        
        result = self.parser.format_compound_title(compound_info, "water")
        
        assert "water (H2O) - CID:962" == result
    
    def test_format_compound_title_long_name(self):
        """Test formatting compound title with long name."""
        compound_info = {
            'cid': 123,
            'molecular_formula': 'C20H30',
            'iupac_name': 'a' * 60  # Very long name
        }
        
        result = self.parser.format_compound_title(compound_info, "long_name")
        
        assert result.endswith("... (C20H30) - CID:123")
        assert len(result.split(" - CID:")[0]) <= 60  # Name part should be truncated