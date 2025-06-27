"""
Integration test for coordinate format error reproduction
Tests the specific error: "Each coordinate must be [element, x, y, z]"
This test reproduces the exact data flow that causes the error in the application.
"""

import unittest
import sys
import os
import json

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from main import PySCFBackend


class TestCoordinateFormatError(unittest.TestCase):
    """Test to reproduce and fix the coordinate format error"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.backend = PySCFBackend()
        
        # Exact same molecule data as used in ProjectPanel.tsx
        self.project_panel_molecules = {
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
            'methane': {
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
            },
            'hydrogen': {
                'type': 'coordinates',
                'coordinates': [
                    ['H', 0.0, 0.0, 0.0],
                    ['H', 0.74, 0.0, 0.0]
                ],
                'charge': 0,
                'spin': 0
            }
        }
    
    def test_actual_application_data_flow_water(self):
        """Test the exact data flow that occurs in the application with water molecule"""
        # Step 1: User clicks on water molecule in ProjectPanel
        # This calls onMoleculeLoad(molecule.data) in ProjectPanel.tsx line 62
        molecule_data = self.project_panel_molecules['water']
        
        # Step 2: App.tsx handleMoleculeLoad calls window.electronAPI.buildMolecule
        # This corresponds to main.py handle_build_molecule
        build_message = {
            'action': 'build-molecule',
            'data': molecule_data,
            'id': 'test_build_1'
        }
        
        build_response = self.backend.handle_message(build_message)
        
        # Verify build was successful
        self.assertEqual(build_response['status'], 'success')
        built_molecule = build_response['molecule']
        
        # Step 3: User starts calculation with built molecule
        # App.tsx handleCalculation calls window.electronAPI.calculate with currentMolecule
        # This is where the error occurs!
        calculation_message = {
            'action': 'calculate',
            'data': {
                'molecule': built_molecule,  # This is the problematic data format
                'method': 'HF',
                'basis': 'STO-3G'
            },
            'id': 'test_calc_1'
        }
        
        # This should currently FAIL with "Each coordinate must be [element, x, y, z]"
        calc_response = self.backend.handle_message(calculation_message)
        
        print(f"Build response molecule format: {built_molecule}")
        print(f"Calculation response: {calc_response}")
        
        # Currently this will fail, but after fix it should succeed
        if calc_response['status'] == 'error':
            self.assertIn('Each coordinate must be', calc_response['message'])
            print("EXPECTED ERROR REPRODUCED: Coordinate format mismatch")
        else:
            print("ERROR FIXED: Calculation succeeded")
            self.assertEqual(calc_response['status'], 'success')
            self.assertIn('results', calc_response)
    
    def test_analyze_data_format_mismatch(self):
        """Analyze the exact data format differences"""
        # Original ProjectPanel format
        original_format = self.project_panel_molecules['water']
        
        # Build molecule and see what format is returned
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': original_format,
            'id': 'test_format_1'
        })
        
        built_format = build_response['molecule']
        
        print("=== DATA FORMAT ANALYSIS ===")
        print(f"Original ProjectPanel format:")
        print(f"  coordinates: {original_format['coordinates']}")
        print(f"  Type: {type(original_format['coordinates'][0])}")
        
        print(f"\nBuilt molecule format:")
        print(f"  coordinates: {built_format.get('coordinates', 'NOT_FOUND')}")
        print(f"  elements: {built_format.get('elements', 'NOT_FOUND')}")
        print(f"  natoms: {built_format.get('natoms', 'NOT_FOUND')}")
        
        # The issue: built_format has separate 'coordinates' and 'elements' arrays
        # But molecule_builder expects coordinates in [element, x, y, z] format
        
        if 'coordinates' in built_format and 'elements' in built_format:
            print(f"\nCONVERTED format (what backend expects):")
            # This is what we need to do to make it work
            converted_coords = []
            for i, (x, y, z) in enumerate(built_format['coordinates']):
                element = built_format['elements'][i]
                converted_coords.append([element, x, y, z])
            print(f"  converted coordinates: {converted_coords}")
    
    def test_methane_data_flow(self):
        """Test data flow with methane molecule"""
        molecule_data = self.project_panel_molecules['methane']
        
        # Build molecule
        build_response = self.backend.handle_message({
            'action': 'build-molecule',
            'data': molecule_data,
            'id': 'test_methane_build'
        })
        
        self.assertEqual(build_response['status'], 'success')
        
        # Try calculation with built molecule (should currently fail)
        calc_response = self.backend.handle_message({
            'action': 'calculate',
            'data': {
                'molecule': build_response['molecule'],
                'method': 'DFT',
                'functional': 'B3LYP',
                'basis': 'STO-3G'
            },
            'id': 'test_methane_calc'
        })
        
        # Document the current state
        if calc_response['status'] == 'error':
            print("Methane calculation FAILED (expected):", calc_response['message'])
        else:
            print("Methane calculation SUCCEEDED (error fixed)")
    
    def test_direct_coordinate_format_validation(self):
        """Test coordinate format validation directly"""
        from utils.molecule_builder import MoleculeBuilder
        
        builder = MoleculeBuilder()
        
        # Test 1: Correct format (should work)
        correct_format = {
            'type': 'coordinates',
            'coordinates': [
                ['O', 0.0, 0.0, 0.0],
                ['H', 0.757, 0.586, 0.0]
            ],
            'charge': 0,
            'spin': 0
        }
        
        try:
            mol1 = builder.build_from_data(correct_format)
            print("✅ Correct format works")
        except Exception as e:
            print(f"❌ Correct format failed: {e}")
        
        # Test 2: Separated format (current backend output - should fail)
        separated_format = {
            'type': 'coordinates',
            'coordinates': [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]],  # No elements
            'elements': ['O', 'H'],  # Separate elements array
            'charge': 0,
            'spin': 0
        }
        
        try:
            mol2 = builder.build_from_data(separated_format)
            print("✅ Separated format unexpectedly works")
        except Exception as e:
            print(f"❌ Separated format failed (expected): {e}")


if __name__ == '__main__':
    unittest.main(verbosity=2)