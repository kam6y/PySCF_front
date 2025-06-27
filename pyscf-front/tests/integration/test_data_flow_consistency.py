"""
Data Flow Consistency Integration Tests
Tests the complete end-to-end data flow from ProjectPanel to Calculation Results
This ensures data format consistency across the entire application workflow.
"""

import unittest
import sys
import os
import json

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from main import PySCFBackend
from utils.molecule_builder import MoleculeBuilder


class TestDataFlowConsistency(unittest.TestCase):
    """Test complete data flow consistency across the application"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.backend = PySCFBackend()
        self.molecule_builder = MoleculeBuilder()
        
        # Test molecules with various complexity levels
        self.test_molecules = {
            'simple_h2': {
                'type': 'coordinates',
                'coordinates': [
                    ['H', 0.0, 0.0, 0.0],
                    ['H', 0.74, 0.0, 0.0]
                ],
                'charge': 0,
                'spin': 0
            },
            'water': {
                'type': 'coordinates',
                'coordinates': [
                    ['O', 0.0, 0.0, 0.0],
                    ['H', 0.757, 0.586, 0.0],
                    ['H', -0.757, 0.586, 0.0]
                ],
                'charge': 0,
                'spin': 0
            },
            'charged_ion': {
                'type': 'coordinates',
                'coordinates': [
                    ['H', 0.0, 0.0, 0.0],
                    ['H', 0.74, 0.0, 0.0]
                ],
                'charge': 1,
                'spin': 1
            }
        }
    
    def test_complete_workflow_h2_hf(self):
        """Test complete workflow: ProjectPanel → Build → Calculate (H2, HF)"""
        molecule_name = 'simple_h2'
        molecule_data = self.test_molecules[molecule_name]
        
        print(f"\n=== Testing complete workflow with {molecule_name} ===")
        
        # Step 1: User selects molecule in ProjectPanel
        print("Step 1: ProjectPanel molecule selection")
        print(f"  Input format: {molecule_data}")
        
        # Step 2: Build molecule (equivalent to handleMoleculeLoad in App.tsx)
        print("Step 2: Building molecule...")
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': molecule_data,
            'id': f'build_{molecule_name}'
        })
        
        print(f"  Build status: {build_response['status']}")
        self.assertEqual(build_response['status'], 'success')
        
        built_molecule = build_response['molecule']
        print(f"  Built molecule format: {built_molecule}")
        
        # Step 3: Perform calculation (equivalent to handleCalculation in App.tsx)
        print("Step 3: Performing calculation...")
        calc_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,
                'method': 'HF',
                'basis': 'STO-3G'
            },
            'id': f'calc_{molecule_name}'
        })
        
        print(f"  Calculation status: {calc_response['status']}")
        if calc_response['status'] == 'error':
            print(f"  ERROR: {calc_response['message']}")
            
        # This is our main assertion - the workflow should work seamlessly
        self.assertEqual(calc_response['status'], 'success', 
                        "Complete workflow should succeed without coordinate format errors")
        
        results = calc_response['results']
        print(f"  Results: energy={results.get('energy')}, converged={results.get('converged')}")
        
        # Verify results make sense
        self.assertIn('energy', results)
        self.assertTrue(results['converged'])
        self.assertEqual(results['method'], 'HF')
    
    def test_complete_workflow_water_dft(self):
        """Test complete workflow: ProjectPanel → Build → Calculate (Water, DFT)"""
        molecule_data = self.test_molecules['water']
        
        print(f"\n=== Testing water molecule DFT workflow ===")
        
        # Build molecule
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': molecule_data,
            'id': 'build_water'
        })
        
        self.assertEqual(build_response['status'], 'success')
        
        # Calculate with DFT
        calc_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': build_response['molecule'],
                'method': 'DFT',
                'functional': 'B3LYP',
                'basis': 'STO-3G'
            },
            'id': 'calc_water_dft'
        })
        
        print(f"  DFT Calculation status: {calc_response['status']}")
        if calc_response['status'] == 'error':
            print(f"  ERROR: {calc_response['message']}")
            
        self.assertEqual(calc_response['status'], 'success')
        
        results = calc_response['results']
        self.assertEqual(results['method'], 'DFT')
        self.assertEqual(results['functional'], 'B3LYP')
        self.assertIn('dipole_moment', results)
    
    def test_data_format_preservation(self):
        """Test that essential data is preserved throughout the workflow"""
        original_data = self.test_molecules['water']
        
        # Build molecule
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': original_data,
            'id': 'test_preservation'
        })
        
        built_molecule = build_response['molecule']
        
        # Check that key properties are preserved
        self.assertEqual(built_molecule['natoms'], len(original_data['coordinates']))
        self.assertEqual(built_molecule['charge'], original_data['charge'])
        self.assertEqual(built_molecule['spin'], original_data['spin'])
        
        # Check that elements are correctly identified (extract from coordinates)
        original_elements = [coord[0] for coord in original_data['coordinates']]
        built_elements = [coord[0] for coord in built_molecule['coordinates']]
        self.assertEqual(built_elements, original_elements)
        
        # Check that coordinate count matches
        self.assertEqual(len(built_molecule['coordinates']), len(original_data['coordinates']))
    
    def test_coordinate_format_conversion(self):
        """Test that coordinate formats are properly converted for calculations"""
        molecule_data = self.test_molecules['simple_h2']
        
        # Build molecule
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': molecule_data,
            'id': 'test_conversion'
        })
        
        built_molecule = build_response['molecule']
        
        # The built molecule should have coordinates in [element, x, y, z] format  
        # This should work seamlessly for calculations
        print(f"Built molecule coordinates format: {built_molecule['coordinates']}")
        
        # Key test: This conversion should work internally
        calc_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,
                'method': 'HF',
                'basis': 'STO-3G'
            },
            'id': 'test_conversion_calc'
        })
        
        if calc_response['status'] == 'error':
            print(f"Conversion test FAILED: {calc_response['message']}")
            
        self.assertEqual(calc_response['status'], 'success', 
                        "Coordinate format conversion should work seamlessly")
    
    def test_multiple_calculation_types(self):
        """Test that all calculation types work with the same built molecule"""
        molecule_data = self.test_molecules['simple_h2']
        
        # Build molecule once
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': molecule_data,
            'id': 'test_multi'
        })
        
        built_molecule = build_response['molecule']
        
        # Test HF
        hf_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,
                'method': 'HF',
                'basis': 'STO-3G'
            },
            'id': 'test_multi_hf'
        })
        
        # Test DFT
        dft_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,
                'method': 'DFT',
                'functional': 'B3LYP',
                'basis': 'STO-3G'
            },
            'id': 'test_multi_dft'
        })
        
        # Test MP2
        mp2_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,
                'method': 'MP2',
                'basis': 'STO-3G'
            },
            'id': 'test_multi_mp2'
        })
        
        # All should succeed
        for method, response in [('HF', hf_response), ('DFT', dft_response), ('MP2', mp2_response)]:
            print(f"{method} calculation status: {response['status']}")
            if response['status'] == 'error':
                print(f"  {method} ERROR: {response['message']}")
            self.assertEqual(response['status'], 'success', f"{method} calculation should succeed")
    
    def test_charged_molecule_workflow(self):
        """Test workflow with charged molecules"""
        molecule_data = self.test_molecules['charged_ion']
        
        # Build charged molecule
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': molecule_data,
            'id': 'build_charged'
        })
        
        self.assertEqual(build_response['status'], 'success')
        
        built_molecule = build_response['molecule']
        self.assertEqual(built_molecule['charge'], 1)
        
        # Calculate charged molecule
        calc_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,
                'method': 'HF',
                'basis': 'STO-3G'
            },
            'id': 'calc_charged'
        })
        
        self.assertEqual(calc_response['status'], 'success')
        
        results = calc_response['results']
        # H2+ should have 1 electron
        self.assertEqual(results['num_electrons'], 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)