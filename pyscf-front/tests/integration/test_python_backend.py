"""
Integration tests for Python backend
Tests the complete Python backend workflow including JSON communication
"""

import unittest
import sys
import os
import json
import subprocess
import time
import threading
from typing import Dict, Any

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from main import PySCFBackend


class TestPythonBackendIntegration(unittest.TestCase):
    """Integration tests for Python backend"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.backend = PySCFBackend()
        
        # Test molecules
        self.test_molecules = {
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
    
    def test_molecule_building_message(self):
        """Test molecule building through message interface"""
        message = {
            'action': 'build-molecule',
            'data': self.test_molecules['water'],
            'id': 'test_msg_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['id'], 'test_msg_1')
        self.assertEqual(response['action'], 'build-molecule')
        self.assertEqual(response['status'], 'success')
        
        molecule = response['molecule']
        self.assertEqual(molecule['natoms'], 3)
        self.assertEqual(molecule['charge'], 0)
        self.assertEqual(molecule['spin'], 0)
        self.assertEqual(len(molecule['elements']), 3)
        self.assertEqual(molecule['elements'], ['O', 'H', 'H'])
    
    def test_hf_calculation_message(self):
        """Test HF calculation through message interface"""
        message = {
            'action': 'calculate',
            'data': {
                'molecule': self.test_molecules['hydrogen'],
                'method': 'HF',
                'basis': 'STO-3G'
            },
            'id': 'test_calc_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['id'], 'test_calc_1')
        self.assertEqual(response['action'], 'calculate')
        self.assertEqual(response['status'], 'success')
        
        results = response['results']
        self.assertEqual(results['method'], 'HF')
        self.assertTrue(results['converged'])
        self.assertLess(results['energy'], 0.0)
        self.assertEqual(results['num_electrons'], 2)
    
    def test_dft_calculation_message(self):
        """Test DFT calculation through message interface"""
        message = {
            'action': 'calculate',
            'data': {
                'molecule': self.test_molecules['water'],
                'method': 'DFT',
                'functional': 'B3LYP',
                'basis': 'STO-3G'
            },
            'id': 'test_dft_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'success')
        
        results = response['results']
        self.assertEqual(results['method'], 'DFT')
        self.assertEqual(results['functional'], 'B3LYP')
        self.assertTrue(results['converged'])
        self.assertIn('dipole_moment', results)
    
    def test_status_message(self):
        """Test status query through message interface"""
        message = {
            'action': 'get-status',
            'data': {},
            'id': 'test_status_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'success')
        self.assertIn('backend_status', response)
        self.assertEqual(response['backend_status'], 'running')
        self.assertIn('jobs', response)
    
    def test_invalid_action_message(self):
        """Test handling of invalid action"""
        message = {
            'action': 'invalid_action',
            'data': {},
            'id': 'test_invalid_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'error')
        self.assertIn('Unknown action', response['message'])
    
    def test_malformed_message(self):
        """Test handling of malformed message"""
        message = {
            'invalid_field': 'value'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'error')
    
    def test_calculation_error_handling(self):
        """Test error handling in calculations"""
        message = {
            'action': 'calculate',
            'data': {
                'molecule': self.test_molecules['water'],
                'method': 'INVALID_METHOD',
                'basis': 'STO-3G'
            },
            'id': 'test_error_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'error')
        self.assertIn('message', response)
    
    def test_molecule_building_error(self):
        """Test error handling in molecule building"""
        invalid_molecule = {
            'type': 'coordinates',
            'coordinates': [],  # Empty coordinates
            'charge': 0,
            'spin': 0
        }
        
        message = {
            'action': 'build-molecule',
            'data': invalid_molecule,
            'id': 'test_mol_error_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'error')
    
    def test_xyz_format_molecule(self):
        """Test building molecule from XYZ format"""
        xyz_data = {
            'type': 'xyz',
            'xyz': '''2
Hydrogen molecule
H 0.000000 0.000000 0.000000
H 0.740000 0.000000 0.000000''',
            'charge': 0,
            'spin': 0
        }
        
        message = {
            'action': 'build-molecule',
            'data': xyz_data,
            'id': 'test_xyz_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'success')
        molecule = response['molecule']
        self.assertEqual(molecule['natoms'], 2)
        self.assertEqual(molecule['elements'], ['H', 'H'])
    
    def test_mp2_calculation(self):
        """Test MP2 calculation"""
        message = {
            'action': 'calculate',
            'data': {
                'molecule': self.test_molecules['hydrogen'],
                'method': 'MP2',
                'basis': 'STO-3G'
            },
            'id': 'test_mp2_1'
        }
        
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['status'], 'success')
        
        results = response['results']
        self.assertEqual(results['method'], 'MP2')
        self.assertIn('hf_energy', results)
        self.assertIn('mp2_correction', results)
        self.assertIn('total_energy', results)
        
        # MP2 correction should be negative
        self.assertLess(results['mp2_correction'], 0.0)


class TestBackendProcessCommunication(unittest.TestCase):
    """Test communication with backend as separate process"""
    
    def setUp(self):
        """Set up test process communication"""
        self.python_executable = sys.executable
        self.backend_script = os.path.join(
            os.path.dirname(__file__), 
            '../../src/python/main.py'
        )
    
    def test_subprocess_communication(self):
        """Test basic subprocess communication"""
        # Start backend process
        process = subprocess.Popen(
            [self.python_executable, self.backend_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            # Send a simple status request
            message = {
                'action': 'get-status',
                'data': {},
                'id': 'subprocess_test_1'
            }
            
            message_str = json.dumps(message) + '\n'
            process.stdin.write(message_str)
            process.stdin.flush()
            
            # Read response
            response_line = process.stdout.readline()
            response = json.loads(response_line)
            
            self.assertEqual(response['id'], 'subprocess_test_1')
            self.assertEqual(response['status'], 'success')
            self.assertEqual(response['backend_status'], 'running')
            
        finally:
            process.terminate()
            process.wait()
    
    def test_subprocess_calculation(self):
        """Test calculation through subprocess"""
        process = subprocess.Popen(
            [self.python_executable, self.backend_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            # Send calculation request
            message = {
                'action': 'calculate',
                'data': {
                    'molecule': {
                        'type': 'coordinates',
                        'coordinates': [
                            ['H', 0.0, 0.0, 0.0],
                            ['H', 0.74, 0.0, 0.0]
                        ],
                        'charge': 0,
                        'spin': 0
                    },
                    'method': 'HF',
                    'basis': 'STO-3G'
                },
                'id': 'subprocess_calc_1'
            }
            
            message_str = json.dumps(message) + '\n'
            process.stdin.write(message_str)
            process.stdin.flush()
            
            # Read response (may take a few seconds)
            response_line = process.stdout.readline()
            response = json.loads(response_line)
            
            self.assertEqual(response['id'], 'subprocess_calc_1')
            self.assertEqual(response['status'], 'success')
            
            results = response['results']
            self.assertEqual(results['method'], 'HF')
            self.assertTrue(results['converged'])
            
        finally:
            process.terminate()
            process.wait()


if __name__ == '__main__':
    unittest.main()