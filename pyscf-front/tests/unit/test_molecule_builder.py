"""
Unit tests for MoleculeBuilder class
Tests molecule building from various input formats
"""

import unittest
import sys
import os
from typing import Dict, Any

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from utils.molecule_builder import MoleculeBuilder
from pyscf import gto


class TestMoleculeBuilder(unittest.TestCase):
    """Test cases for MoleculeBuilder class"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.builder = MoleculeBuilder()
    
    def test_build_from_coordinates_water(self):
        """Test building water molecule from coordinates"""
        water_data = {
            'type': 'coordinates',
            'coordinates': [
                ['O', 0.0, 0.0, 0.0],
                ['H', 0.757, 0.586, 0.0],
                ['H', -0.757, 0.586, 0.0]
            ],
            'charge': 0,
            'spin': 0
        }
        
        mol = self.builder.build_from_data(water_data)
        
        # Verify molecule properties
        self.assertIsInstance(mol, gto.Mole)
        self.assertEqual(mol.natm, 3)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.spin, 0)
        self.assertEqual(mol.nelectron, 10)
        
        # Check atom symbols
        symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
        self.assertEqual(symbols, ['O', 'H', 'H'])
    
    def test_build_from_coordinates_methane(self):
        """Test building methane molecule from coordinates"""
        methane_data = {
            'type': 'coordinates',
            'coordinates': [
                ['C', 0.0, 0.0, 0.0],
                ['H', 1.089, 0.0, 0.0],
                ['H', -0.363, 1.027, 0.0],
                ['H', -0.363, -0.513, 0.889],
                ['H', -0.363, -0.513, -0.889]
            ],
            'charge': 0,
            'spin': 0
        }
        
        mol = self.builder.build_from_data(methane_data)
        
        self.assertEqual(mol.natm, 5)
        self.assertEqual(mol.nelectron, 10)
        
        # Check atom symbols
        symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
        self.assertEqual(symbols, ['C', 'H', 'H', 'H', 'H'])
    
    def test_build_from_coordinates_charged(self):
        """Test building charged molecule"""
        data = {
            'type': 'coordinates',
            'coordinates': [
                ['H', 0.0, 0.0, 0.0],
                ['H', 0.74, 0.0, 0.0]
            ],
            'charge': 1,
            'spin': 1
        }
        
        mol = self.builder.build_from_data(data)
        
        self.assertEqual(mol.natm, 2)
        self.assertEqual(mol.charge, 1)
        self.assertEqual(mol.spin, 1)
        self.assertEqual(mol.nelectron, 1)  # H2+ has 1 electron
    
    def test_build_from_xyz_format(self):
        """Test building molecule from XYZ format string"""
        xyz_string = """3
Water molecule
O 0.000000 0.000000 0.000000
H 0.757000 0.586000 0.000000
H -0.757000 0.586000 0.000000"""
        
        data = {
            'type': 'xyz',
            'xyz': xyz_string,
            'charge': 0,
            'spin': 0
        }
        
        mol = self.builder.build_from_data(data)
        
        self.assertEqual(mol.natm, 3)
        self.assertEqual(mol.charge, 0)
        self.assertEqual(mol.spin, 0)
        
        # Verify coordinates are approximately correct
        # Note: mol.atom_coords() returns coordinates in Bohr units
        # 0.757 Angstrom = 0.757 / 0.5291772083 = 1.430 Bohr (approximately)
        coords = mol.atom_coords()
        self.assertAlmostEqual(coords[0][0], 0.0, places=5)  # O x-coord
        self.assertAlmostEqual(coords[1][0], 1.4305, places=3)  # H x-coord in Bohr
        
        # Also test that we can get coordinates back in Angstrom
        coords_angstrom = self.builder.get_coordinates_in_angstrom(mol)
        self.assertAlmostEqual(coords_angstrom[1][1], 0.757, places=3)  # H x-coord in Angstrom
    
    def test_get_test_molecules(self):
        """Test getting predefined test molecules"""
        test_molecules = self.builder.get_test_molecules()
        
        # Should contain expected molecules
        self.assertIn('water', test_molecules)
        self.assertIn('methane', test_molecules)
        self.assertIn('hydrogen', test_molecules)
        
        # Each molecule should have required fields
        for name, mol_data in test_molecules.items():
            self.assertIn('type', mol_data)
            self.assertIn('coordinates', mol_data)
            self.assertIn('charge', mol_data)
            self.assertIn('spin', mol_data)
            
            # Verify we can build each test molecule
            mol = self.builder.build_from_data(mol_data)
            self.assertIsInstance(mol, gto.Mole)
            self.assertGreater(mol.natm, 0)
    
    def test_invalid_input_type(self):
        """Test handling of invalid input type"""
        data = {
            'type': 'invalid_type',
            'data': 'some data'
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(data)
    
    def test_empty_coordinates(self):
        """Test handling of empty coordinates"""
        data = {
            'type': 'coordinates',
            'coordinates': [],
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(data)
    
    def test_invalid_coordinate_format(self):
        """Test handling of invalid coordinate format"""
        data = {
            'type': 'coordinates',
            'coordinates': [
                ['O', 0.0, 0.0],  # Missing z-coordinate
                ['H', 0.757, 0.586, 0.0],
                ['H', -0.757, 0.586, 0.0]
            ],
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(data)
    
    def test_invalid_xyz_format(self):
        """Test handling of invalid XYZ format"""
        invalid_xyz = """2
Comment line
O 0.0 0.0"""  # Missing z-coordinate
        
        data = {
            'type': 'xyz',
            'xyz': invalid_xyz,
            'charge': 0,
            'spin': 0
        }
        
        with self.assertRaises(ValueError):
            self.builder.build_from_data(data)


if __name__ == '__main__':
    unittest.main()