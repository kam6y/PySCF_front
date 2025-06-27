#!/usr/bin/env python3
"""
PySCF Frontend Python Backend
Main entry point for the Python backend that handles quantum chemistry calculations.
"""

import sys
import json
import logging
from typing import Dict, Any, Optional
from calculations.engine import CalculationEngine
from utils.molecule_builder import MoleculeBuilder
from utils.job_manager import JobManager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PySCFBackend:
    """Main backend class that handles communication with Electron frontend."""
    
    def __init__(self):
        self.calculation_engine = CalculationEngine()
        self.molecule_builder = MoleculeBuilder()
        self.job_manager = JobManager()
        self.running = True
        
        logger.info("PySCF Backend initialized successfully")
    
    def handle_message(self, message: Dict[str, Any]) -> Dict[str, Any]:
        """Handle incoming messages from Electron frontend."""
        try:
            action = message.get('action')
            data = message.get('data', {})
            message_id = message.get('id')
            
            logger.info(f"Received action: {action}")
            
            if action == 'calculate':
                result = self.handle_calculate(data)
            elif action == 'build-molecule':
                result = self.handle_build_molecule(data)
            elif action == 'get-status':
                result = self.handle_get_status()
            elif action == 'stop':
                result = self.handle_stop()
            else:
                result = {
                    'status': 'error',
                    'message': f'Unknown action: {action}'
                }
            
            return {
                'id': message_id,
                'action': action,
                **result
            }
            
        except Exception as e:
            logger.error(f"Error handling message: {e}")
            return {
                'id': message.get('id'),
                'status': 'error',
                'message': str(e)
            }
    
    def handle_calculate(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Handle calculation requests."""
        try:
            # Extract calculation parameters
            molecule_data = data.get('molecule')
            method = data.get('method', 'DFT')
            functional = data.get('functional', 'B3LYP')
            basis = data.get('basis', '6-31G*')
            
            # Build molecule
            mol = self.molecule_builder.build_from_data(molecule_data)
            
            # Perform calculation
            result = self.calculation_engine.calculate(
                mol=mol,
                method=method,
                functional=functional,
                basis=basis
            )
            
            return {
                'status': 'success',
                'results': result
            }
            
        except Exception as e:
            logger.error(f"Calculation error: {e}")
            return {
                'status': 'error',
                'message': str(e)
            }
    
    def handle_build_molecule(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Handle molecule building requests."""
        try:
            mol = self.molecule_builder.build_from_data(data)
            
            return {
                'status': 'success',
                'molecule': {
                    'natoms': mol.natm,
                    'charge': mol.charge,
                    'spin': mol.spin,
                    'coordinates': mol.atom_coords().tolist(),
                    'elements': [mol.atom_symbol(i) for i in range(mol.natm)]
                }
            }
            
        except Exception as e:
            logger.error(f"Molecule building error: {e}")
            return {
                'status': 'error',
                'message': str(e)
            }
    
    def handle_get_status(self) -> Dict[str, Any]:
        """Handle status requests."""
        return {
            'status': 'success',
            'backend_status': 'running',
            'jobs': self.job_manager.get_status()
        }
    
    def handle_stop(self) -> Dict[str, Any]:
        """Handle stop requests."""
        self.running = False
        return {
            'status': 'success',
            'message': 'Backend stopping'
        }
    
    def run(self):
        """Main loop for handling messages from stdin."""
        logger.info("Backend ready, waiting for messages...")
        
        try:
            while self.running:
                # Read from stdin
                line = sys.stdin.readline()
                if not line:
                    break
                
                try:
                    line_stripped = line.strip()
                    if not line_stripped:
                        continue  # Skip empty lines
                    
                    # Parse JSON with simplified error handling
                    message = json.loads(line_stripped)
                    
                    # Basic validation - must be a dictionary
                    if not isinstance(message, dict):
                        logger.warning(f"Received non-dict message: {type(message)}")
                        continue
                    
                    response = self.handle_message(message)
                    
                    # Send response to stdout with proper formatting
                    response_json = json.dumps(response, separators=(',', ':'))
                    print(response_json, flush=True)
                    
                except json.JSONDecodeError as e:
                    logger.error(f"JSON decode error: {e}")
                    # Send simplified error response
                    error_response = {
                        'status': 'error',
                        'message': 'Invalid JSON format'
                    }
                    print(json.dumps(error_response), flush=True)
                except Exception as e:
                    logger.error(f"Error processing message: {e}")
                    # Send general error response
                    error_response = {
                        'status': 'error',
                        'message': 'Message processing error'
                    }
                    print(json.dumps(error_response), flush=True)
                    
        except KeyboardInterrupt:
            logger.info("Backend interrupted by user")
        except Exception as e:
            logger.error(f"Unexpected error in main loop: {e}")
        finally:
            logger.info("Backend shutting down")


if __name__ == '__main__':
    backend = PySCFBackend()
    backend.run()