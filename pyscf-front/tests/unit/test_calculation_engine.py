"""
Unit tests for CalculationEngine class
Tests quantum chemistry calculations using PySCF
"""

import unittest
import sys
import os
from typing import Dict, Any

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from calculations.engine import CalculationEngine
from utils.molecule_builder import MoleculeBuilder
from pyscf import gto


class TestCalculationEngine(unittest.TestCase):
    """Test cases for CalculationEngine class"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.engine = CalculationEngine()
        self.builder = MoleculeBuilder()
        
        # Create test molecules
        self.water_mol = self.builder.build_from_data(
            self.builder.get_test_molecules()['water']
        )
        self.hydrogen_mol = self.builder.build_from_data(
            self.builder.get_test_molecules()['hydrogen']
        )
    
    def test_hf_calculation_hydrogen(self):
        """Test Hartree-Fock calculation on H2"""
        result = self.engine.calculate(
            mol=self.hydrogen_mol,
            method='HF',
            basis='STO-3G'
        )
        
        # Verify result structure
        self.assertIsInstance(result, dict)
        self.assertIn('method', result)
        self.assertIn('energy', result)
        self.assertIn('converged', result)
        self.assertIn('num_electrons', result)
        self.assertIn('multiplicity', result)
        
        # Verify values
        self.assertEqual(result['method'], 'HF')
        self.assertTrue(result['converged'])
        self.assertEqual(result['num_electrons'], 2)
        self.assertEqual(result['multiplicity'], 1)
        
        # Energy should be reasonable for H2/STO-3G
        self.assertLess(result['energy'], -1.0)  # Should be negative
        self.assertGreater(result['energy'], -2.0)  # Not too negative
    
    def test_hf_calculation_water(self):
        """Test Hartree-Fock calculation on water"""
        result = self.engine.calculate(
            mol=self.water_mol,
            method='HF',
            basis='STO-3G'
        )
        
        # Verify basic structure
        self.assertEqual(result['method'], 'HF')
        self.assertTrue(result['converged'])
        self.assertEqual(result['num_electrons'], 10)
        self.assertEqual(result['multiplicity'], 1)
        
        # Water HF/STO-3G energy should be around -74.96 Hartree
        self.assertLess(result['energy'], -70.0)
        self.assertGreater(result['energy'], -80.0)
        
        # Should have orbital information
        self.assertIn('homo_energy', result)
        self.assertIn('orbital_energies', result)
    
    def test_dft_calculation_water(self):
        """Test DFT calculation on water"""
        result = self.engine.calculate(
            mol=self.water_mol,
            method='DFT',
            functional='B3LYP',
            basis='STO-3G'
        )
        
        # Verify DFT-specific fields
        self.assertEqual(result['method'], 'DFT')
        self.assertEqual(result['functional'], 'B3LYP')
        self.assertTrue(result['converged'])
        
        # Should have dipole moment for water
        self.assertIn('dipole_moment', result)
        self.assertIsInstance(result['dipole_moment'], list)
        self.assertEqual(len(result['dipole_moment']), 3)
        
        # Water should have significant dipole moment
        if 'dipole_magnitude' in result:
            self.assertGreater(result['dipole_magnitude'], 1.0)
    
    def test_mp2_calculation_hydrogen(self):
        """Test MP2 calculation on H2"""
        result = self.engine.calculate(
            mol=self.hydrogen_mol,
            method='MP2',
            basis='STO-3G'
        )
        
        # Verify MP2-specific fields
        self.assertEqual(result['method'], 'MP2')
        self.assertIn('hf_energy', result)
        self.assertIn('mp2_correction', result)
        self.assertIn('total_energy', result)
        
        # Total energy should be sum of HF + MP2 correction
        expected_total = result['hf_energy'] + result['mp2_correction']
        self.assertAlmostEqual(result['total_energy'], expected_total, places=8)
        
        # MP2 correction should be negative (correlation energy)
        self.assertLess(result['mp2_correction'], 0.0)
    
    def test_different_basis_sets(self):
        """Test calculations with different basis sets"""
        basis_sets = ['STO-3G', '6-31G']
        
        for basis in basis_sets:
            with self.subTest(basis=basis):
                result = self.engine.calculate(
                    mol=self.hydrogen_mol,
                    method='HF',
                    basis=basis
                )
                
                self.assertTrue(result['converged'])
                self.assertLess(result['energy'], 0.0)
    
    def test_different_functionals(self):
        """Test DFT with different functionals"""
        functionals = ['B3LYP', 'PBE']
        
        for functional in functionals:
            with self.subTest(functional=functional):
                result = self.engine.calculate(
                    mol=self.water_mol,
                    method='DFT',
                    functional=functional,
                    basis='STO-3G'
                )
                
                self.assertEqual(result['functional'], functional)
                self.assertTrue(result['converged'])
                self.assertLess(result['energy'], 0.0)
    
    def test_invalid_method(self):
        """Test handling of invalid calculation method"""
        with self.assertRaises(ValueError):
            self.engine.calculate(
                mol=self.water_mol,
                method='INVALID_METHOD',
                basis='STO-3G'
            )
    
    def test_convergence_options(self):
        """Test custom convergence criteria"""
        result = self.engine.calculate(
            mol=self.hydrogen_mol,
            method='HF',
            basis='STO-3G',
            conv_tol=1e-10,
            max_cycle=100
        )
        
        self.assertTrue(result['converged'])
        self.assertLess(result['energy'], 0.0)
    
    def test_orbital_energies_structure(self):
        """Test orbital energies output structure"""
        result = self.engine.calculate(
            mol=self.water_mol,
            method='HF',
            basis='STO-3G'
        )
        
        if 'orbital_energies' in result:
            self.assertIsInstance(result['orbital_energies'], list)
            self.assertGreater(len(result['orbital_energies']), 0)
            
            # All orbital energies should be numbers
            for energy in result['orbital_energies']:
                self.assertIsInstance(energy, (int, float))
        
        # HOMO energy should be reasonable
        if 'homo_energy' in result:
            self.assertIsInstance(result['homo_energy'], (int, float))
            self.assertLess(result['homo_energy'], 0.0)  # Should be negative
        
        # LUMO-HOMO gap should be positive if both exist
        if 'homo_lumo_gap' in result:
            self.assertGreater(result['homo_lumo_gap'], 0.0)


if __name__ == '__main__':
    unittest.main()