"""
Coordinate Format Compatibility Tests
Tests for handling multiple coordinate formats in MoleculeBuilder
This ensures the molecule builder can handle various input formats flexibly.
"""

import unittest
import sys
import os

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from utils.molecule_builder import MoleculeBuilder
from pyscf import gto


class TestCoordinateFormatCompatibility(unittest.TestCase):
    """Test multiple coordinate format handling"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.builder = MoleculeBuilder()
    
    def test_standard_format_with_elements(self):
        """Test standard [element, x, y, z] format"""
        data = {
            'type': 'coordinates',
            'coordinates': [
                ['O', 0.0, 0.0, 0.0],
                ['H', 0.757, 0.586, 0.0],
                ['H', -0.757, 0.586, 0.0]
            ],
            'charge': 0,
            'spin': 0
        }
        
        mol = self.builder.build_from_data(data)
        
        self.assertEqual(mol.natm, 3)
        self.assertEqual([mol.atom_symbol(i) for i in range(mol.natm)], ['O', 'H', 'H'])
    
    def test_separated_coordinates_and_elements(self):
        """Test separated coordinates [[x,y,z], ...] and elements ['O', 'H', ...]"""
        # This format is what build-molecule currently returns
        data = {
            'type': 'coordinates',
            'coordinates': [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
            'elements': ['O', 'H', 'H'],
            'charge': 0,
            'spin': 0
        }
        
        # This should work after our fix
        mol = self.builder.build_from_data(data)
        
        self.assertEqual(mol.natm, 3)
        self.assertEqual([mol.atom_symbol(i) for i in range(mol.natm)], ['O', 'H', 'H'])
    
    def test_mixed_format_detection(self):
        """Test detection of different coordinate formats"""
        # Standard format
        standard_data = {
            'coordinates': [['O', 0.0, 0.0, 0.0], ['H', 0.757, 0.586, 0.0]]
        }
        
        # Separated format  
        separated_data = {
            'coordinates': [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]],
            'elements': ['O', 'H']
        }
        
        # Test format detection logic
        # Standard format: coordinates[0] is a list with 4 elements, first is string
        is_standard = (len(standard_data['coordinates'][0]) == 4 and 
                      isinstance(standard_data['coordinates'][0][0], str))
        self.assertTrue(is_standard)
        
        # Separated format: coordinates[0] has 3 elements (no element), and elements key exists
        is_separated = (len(separated_data['coordinates'][0]) == 3 and 
                       'elements' in separated_data)
        self.assertTrue(is_separated)
    
    def test_coordinate_conversion_utilities(self):
        """Test coordinate format conversion utilities"""
        # Convert from separated to standard format
        separated_coords = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]]
        elements = ['O', 'H']
        
        # This conversion should be done internally
        standard_coords = []
        for i, (x, y, z) in enumerate(separated_coords):
            standard_coords.append([elements[i], x, y, z])
        
        expected = [['O', 0.0, 0.0, 0.0], ['H', 0.757, 0.586, 0.0]]
        self.assertEqual(standard_coords, expected)
    
    def test_bohr_angstrom_handling(self):
        """Test proper handling of Bohr vs Angstrom units"""
        # Water molecule coordinates in Angstrom
        angstrom_data = {
            'type': 'coordinates',
            'coordinates': [
                ['O', 0.0, 0.0, 0.0],
                ['H', 0.757, 0.586, 0.0],
                ['H', -0.757, 0.586, 0.0]
            ],
            'charge': 0,
            'spin': 0,
            'unit': 'Angstrom'
        }
        
        mol = self.builder.build_from_data(angstrom_data)
        
        # Get coordinates back in Angstrom
        coords_angstrom = self.builder.get_coordinates_in_angstrom(mol)
        
        # Should preserve the original values (approximately)
        self.assertAlmostEqual(coords_angstrom[1][1], 0.757, places=3)  # H x-coord
        self.assertAlmostEqual(coords_angstrom[1][2], 0.586, places=3)  # H y-coord
    
    def test_error_handling_invalid_formats(self):
        """Test error handling for invalid coordinate formats"""
        # Missing elements in separated format
        invalid_separated = {
            'type': 'coordinates',
            'coordinates': [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]],
            # Missing 'elements' key
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(invalid_separated)
        
        # Wrong coordinate length
        invalid_length = {
            'type': 'coordinates',
            'coordinates': [
                ['O', 0.0, 0.0],  # Missing z-coordinate
                ['H', 0.757, 0.586, 0.0]
            ],
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(invalid_length)
    
    def test_element_coordinate_count_mismatch(self):
        """Test handling of mismatched element and coordinate counts"""
        mismatched_data = {
            'type': 'coordinates',
            'coordinates': [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]],  # 2 coordinates
            'elements': ['O', 'H', 'H'],  # 3 elements
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(mismatched_data)
    
    def test_empty_coordinates_handling(self):
        """Test handling of empty coordinates"""
        empty_data = {
            'type': 'coordinates',
            'coordinates': [],
            'elements': [],
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(empty_data)
    
    def test_coordinate_format_auto_detection(self):
        """Test automatic detection and handling of coordinate formats"""
        # This test will pass after we implement the auto-detection feature
        
        # Standard format (use valid molecule - water)
        standard_format = {
            'type': 'coordinates',
            'coordinates': [['O', 0.0, 0.0, 0.0], ['H', 0.757, 0.586, 0.0], ['H', -0.757, 0.586, 0.0]],
            'charge': 0,
            'spin': 0
        }
        
        mol1 = self.builder.build_from_data(standard_format)
        
        # Separated format (what build-molecule returns)
        separated_format = {
            'type': 'coordinates',
            'coordinates': [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]],
            'elements': ['O', 'H', 'H'],
            'charge': 0,
            'spin': 0
        }
        
        # This should work after our enhancement
        mol2 = self.builder.build_from_data(separated_format)
        
        # Both should produce equivalent molecules
        self.assertEqual(mol1.natm, mol2.natm)
        self.assertEqual([mol1.atom_symbol(i) for i in range(mol1.natm)],
                        [mol2.atom_symbol(i) for i in range(mol2.natm)])


if __name__ == '__main__':
    unittest.main(verbosity=2)