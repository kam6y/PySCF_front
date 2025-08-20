import logging
import os
import sys
import json
import time
import signal
import atexit
import multiprocessing
import threading
from datetime import datetime
from flask import Flask, request, jsonify
from flask_cors import CORS
from flask_sock import Sock
from flask_pydantic import validate
from pydantic import ValidationError
import socket
import shutil
from typing import Dict
import gevent
from gevent.pywsgi import WSGIServer

from pubchem.client import PubChemClient, PubChemError, PubChemNotFoundError
from pubchem import parser as xyz_parser
from SMILES.smiles_converter import smiles_to_xyz, SMILESError
from quantum_calc import DFTCalculator, HFCalculator, MP2Calculator, CCSDCalculator, TDDFTCalculator, MolecularOrbitalGenerator, CalculationError, ConvergenceError, InputError, GeometryError
from quantum_calc.exceptions import XYZValidationError, FileManagerError, ProcessManagerError, WebSocketError
from quantum_calc.file_manager import CalculationFileManager
from quantum_calc import get_process_manager, shutdown_process_manager, get_websocket_watcher, shutdown_websocket_watcher
from generated_models import (
    PubChemSearchRequest, SMILESConvertRequest, XYZValidateRequest,
    QuantumCalculationRequest, CalculationUpdateRequest
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize Flask app and extensions
app = Flask(__name__)
CORS(app)  # Enable CORS for cross-origin requests
sock = Sock(app)  # Sockインスタンスを作成

# Initialize PubChem client
pubchem_client = PubChemClient(timeout=30)

# Global WebSocket connection registry for immediate notifications
# calculation_id -> set of websocket connections
active_websockets = {}
websocket_lock = threading.Lock()


def send_immediate_websocket_notification(calculation_id: str, status: str, error_message: str = None):
    """Send immediate WebSocket notification to all connected clients for a calculation."""
    with websocket_lock:
        if calculation_id not in active_websockets:
            return
        
        connections = active_websockets[calculation_id].copy()
    
    if not connections:
        return
        
    # Build complete calculation instance
    try:
        from quantum_calc.file_manager import CalculationFileManager
        file_manager = CalculationFileManager()
        calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        if os.path.exists(calc_dir):
            # Read current data from files
            parameters = file_manager.read_calculation_parameters(calc_dir) or {}
            results = file_manager.read_calculation_results(calc_dir)
            display_name = file_manager._get_display_name(calculation_id, parameters)
            
            # Build calculation instance
            calculation_instance = {
                'id': calculation_id,
                'name': display_name,
                'status': status,
                'createdAt': parameters.get('created_at', datetime.now().isoformat()),
                'updatedAt': datetime.now().isoformat(),
                'parameters': parameters,
                'results': results,
                'workingDirectory': calc_dir,
            }
            
            # Add error if provided
            if error_message:
                calculation_instance['error'] = error_message
            
            message = json.dumps(calculation_instance)
            
            # Send to all connected WebSocket clients
            failed_connections = set()
            for ws in connections:
                try:
                    ws.send(message)
                except Exception as e:
                    logger.warning(f"Failed to send immediate notification to WebSocket client: {e}")
                    failed_connections.add(ws)
            
            # Clean up failed connections
            if failed_connections:
                with websocket_lock:
                    if calculation_id in active_websockets:
                        active_websockets[calculation_id] -= failed_connections
                        if not active_websockets[calculation_id]:
                            del active_websockets[calculation_id]
        else:
            logger.warning(f"Calculation directory not found for immediate notification: {calc_dir}")
            
    except Exception as e:
        logger.error(f"Error sending immediate WebSocket notification for {calculation_id}: {e}")


def cleanup_resources():
    """Clean up resources including the process pool and file watcher."""
    logger.info("Cleaning up resources...")
    try:
        shutdown_process_manager()
        logger.info("Process manager shut down successfully")
    except ProcessManagerError as e:
        logger.warning(f"Process manager already shut down or unavailable: {e}")
    except Exception as e:
        logger.error(f"Unexpected error shutting down process manager: {e}", exc_info=True)
    
    try:
        shutdown_websocket_watcher()
        logger.info("WebSocket file watcher shut down successfully")
    except Exception as e:
        logger.error(f"Unexpected error shutting down WebSocket file watcher: {e}", exc_info=True)


def signal_handler(signum, frame):
    """Handle shutdown signals."""
    logger.info(f"Received signal {signum}, shutting down gracefully...")
    cleanup_resources()
    sys.exit(0)


# Register cleanup functions
atexit.register(cleanup_resources)
signal.signal(signal.SIGTERM, signal_handler)
signal.signal(signal.SIGINT, signal_handler)





@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'ok',
        'service': 'pyscf-front-api',
        'version': '0.3.0'
    })


@app.route('/api/pubchem/search', methods=['POST'])
@validate()
def search_pubchem(body: PubChemSearchRequest):
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    try:
        query = body.query
        search_type = body.search_type.value
        
        logger.info(f"Searching PubChem for '{query}' (type: {search_type})")
        
        compound_data = pubchem_client.search_compound(query, search_type)
        if not compound_data or not compound_data.atoms:
            return jsonify({'success': False, 'error': f'No compound with a 3D structure found for query: {query}'}), 404
        
        logger.info(f"Found CID {compound_data.cid} with {len(compound_data.atoms)} atoms.")
        
        title = xyz_parser.format_compound_title(compound_data, query)
        xyz_string = xyz_parser.atoms_to_xyz(compound_data.atoms, title)
        
        logger.info(f"Successfully generated XYZ for CID {compound_data.cid}")
        
        return jsonify({
            'success': True,
            'data': {
                'xyz': xyz_string,
                'compound_info': {
                    'cid': compound_data.cid,
                    'iupac_name': compound_data.iupac_name,
                    'molecular_formula': compound_data.molecular_formula,
                    'molecular_weight': compound_data.molecular_weight,
                    'synonyms': compound_data.synonyms,
                },
                'atom_count': len(compound_data.atoms)
            }
        })
            
    except PubChemNotFoundError as e:
        logger.warning(f"PubChem search failed (Not Found): {e}")
        return jsonify({'success': False, 'error': str(e)}), 404
    except PubChemError as e:
        logger.error(f"A PubChem API error occurred: {e}", exc_info=True)
        return jsonify({'success': False, 'error': str(e)}), e.status_code or 500
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/smiles/convert', methods=['POST'])
@validate()
def convert_smiles(body: SMILESConvertRequest):
    """Converts a SMILES string to XYZ format."""
    try:
        smiles = body.smiles
        
        logger.info(f"Converting SMILES: {smiles}")

        xyz_string = smiles_to_xyz(smiles, title=f"Molecule from SMILES: {smiles}")

        return jsonify({'success': True, 'data': {'xyz': xyz_string}})

    except SMILESError as e:
        logger.error(f"SMILES conversion failed: {e}")
        return jsonify({'success': False, 'error': str(e)}), 400
    except Exception as e:
        logger.error(f"An unexpected error occurred during SMILES conversion: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/pubchem/validate', methods=['POST'])
@validate()
def validate_xyz_endpoint(body: XYZValidateRequest):
    """Validate an XYZ format string."""
    try:
        xyz_string = body.xyz
        
        # Check for empty or None input
        if not xyz_string or not xyz_string.strip():
            return jsonify({'success': False, 'error': 'XYZ string cannot be empty.'}), 400
        
        validation_result = xyz_parser.validate_xyz(xyz_string)
        
        return jsonify({'success': True, 'data': validation_result})
        
    except (TypeError, AttributeError) as e:
        logger.warning(f"Invalid input data for XYZ validation: {e}")
        return jsonify({'success': False, 'error': 'Invalid input format for XYZ validation.'}), 400
    except MemoryError as e:
        logger.error(f"Memory error during XYZ validation: {e}")
        return jsonify({'success': False, 'error': 'XYZ string too large to process.'}), 413
    except Exception as e:
        logger.error(f"Unexpected error during XYZ validation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculate', methods=['POST'])
@validate()
def quantum_calculate(body: QuantumCalculationRequest):
    """
    Starts a quantum chemistry calculation in the background.
    Immediately returns a calculation ID to track the job.
    """
    try:
        # Prepare parameters using validated data from Pydantic model
        parameters = {
            'calculation_method': body.calculation_method.value,
            'basis_function': body.basis_function,
            'exchange_correlation': body.exchange_correlation,
            'charges': body.charges,
            'spin_multiplicity': body.spin_multiplicity,
            'solvent_method': body.solvent_method.value,
            'solvent': body.solvent,
            'xyz': body.xyz,
            'name': body.name,
            'cpu_cores': body.cpu_cores,
            'memory_mb': body.memory_mb,
            'created_at': datetime.now().isoformat()
        }
        
        # Add TDDFT-specific parameters if the method is TDDFT
        if body.calculation_method == 'TDDFT':
            parameters['tddft_nstates'] = body.tddft_nstates
            parameters['tddft_method'] = body.tddft_method.value
            parameters['tddft_analyze_nto'] = body.tddft_analyze_nto
        
        # Initialize file manager and create directory
        file_manager = CalculationFileManager()
        calc_dir = file_manager.create_calculation_dir(parameters['name'])
        calculation_id = os.path.basename(calc_dir)
        
        # Save initial parameters and status
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, 'running')

        # Submit calculation to process pool
        process_manager = get_process_manager()
        success = process_manager.submit_calculation(calculation_id, parameters)
        
        if not success:
            # If submission failed, update status to error
            file_manager.save_calculation_status(calc_dir, 'error')
            file_manager.save_calculation_results(calc_dir, {'error': 'Failed to submit calculation to process pool.'})
            logger.error(f"Failed to submit calculation {calculation_id} to process pool")
            return jsonify({'success': False, 'error': 'Failed to queue calculation.'}), 500
        
        # Register completion callback for immediate WebSocket notification
        def on_calculation_complete(calc_id: str, success: bool, error_message: str):
            """Callback function to send immediate WebSocket notification when calculation completes."""
            final_status = 'completed' if success else 'error'
            send_immediate_websocket_notification(calc_id, final_status, error_message)
        
        process_manager.register_completion_callback(calculation_id, on_calculation_complete)
        
        logger.info(f"Queued calculation {calculation_id} for molecule '{parameters['name']}'")

        # Return immediately with the new calculation instance
        initial_instance = {
            'id': calculation_id,
            'name': parameters['name'],
            'status': 'running', # Frontend polls for actual status, so this is fine
            'createdAt': parameters['created_at'],
            'updatedAt': parameters['created_at'],
            'parameters': parameters
        }

        return jsonify({'success': True, 'data': {'calculation': initial_instance}}), 202

    except (InputError, GeometryError) as e:
        logger.warning(f"Invalid calculation parameters: {e}")
        return jsonify({'success': False, 'error': f'Invalid input parameters: {str(e)}'}), 400
    except FileManagerError as e:
        logger.error(f"File management error during calculation setup: {e}")
        return jsonify({'success': False, 'error': 'Failed to set up calculation files.'}), 500
    except ProcessManagerError as e:
        logger.error(f"Process manager error during calculation submission: {e}")
        return jsonify({'success': False, 'error': 'Calculation system is currently unavailable.'}), 503
    except OSError as e:
        logger.error(f"System error during calculation setup: {e}")
        return jsonify({'success': False, 'error': 'Insufficient system resources to start calculation.'}), 507
    except PermissionError as e:
        logger.error(f"Permission error during calculation setup: {e}")
        return jsonify({'success': False, 'error': 'System permission error. Please contact administrator.'}), 500
    except Exception as e:
        logger.error(f"Unexpected error queuing calculation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations', methods=['GET'])
def list_calculations():
    """List all available calculation directories."""
    try:
        file_manager = CalculationFileManager()
        calculations = file_manager.list_calculations()
        
        return jsonify({
            'success': True,
            'data': {
                'base_directory': file_manager.get_base_directory(),
                'calculations': calculations,
                'count': len(calculations)
            }
        })
        
    except FileManagerError as e:
        logger.error(f"File manager error while listing calculations: {e}")
        return jsonify({'success': False, 'error': 'Unable to access calculation directory.'}), 500
    except PermissionError as e:
        logger.error(f"Permission error while listing calculations: {e}")
        return jsonify({'success': False, 'error': 'Permission denied accessing calculation directory.'}), 403
    except OSError as e:
        logger.error(f"System error while listing calculations: {e}")
        return jsonify({'success': False, 'error': 'System error accessing calculation files.'}), 500
    except Exception as e:
        logger.error(f"Unexpected error listing calculations: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/status', methods=['GET'])
def get_calculation_status():
    """Get status information about the calculation system."""
    try:
        process_manager = get_process_manager()
        active_calculations = process_manager.get_active_calculations()
        
        return jsonify({
            'success': True,
            'data': {
                'process_pool': {
                    'max_workers': process_manager.max_workers,
                    'active_calculations': active_calculations,
                    'active_count': len(active_calculations),
                    'is_shutdown': process_manager._shutdown
                },
                'system': {
                    'cpu_count': multiprocessing.cpu_count()
                }
            }
        })
        
    except ProcessManagerError as e:
        logger.error(f"Process manager error getting status: {e}")
        return jsonify({'success': False, 'error': 'Process manager is unavailable.'}), 503
    except AttributeError as e:
        logger.error(f"Process manager attribute error: {e}")
        return jsonify({'success': False, 'error': 'Process manager not properly initialized.'}), 500
    except Exception as e:
        logger.error(f"Unexpected error getting calculation status: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['GET'])
def get_calculation_details(calculation_id):
    """
    Get detailed information about a specific calculation.
    This now exclusively reads from the file system, making it the single source of truth.
    """
    try:
        
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)

        if not os.path.isdir(calc_path):
             return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404

        # Read all calculation data from disk
        parameters = file_manager.read_calculation_parameters(calc_path) or {}
        results = file_manager.read_calculation_results(calc_path)
        status = file_manager.read_calculation_status(calc_path)
        
        display_name = file_manager._get_display_name(calculation_id, parameters)
        creation_date = parameters.get('created_at', datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat())

        calculation_instance = {
            'id': calculation_id,
            'name': display_name,
            'status': status,
            'createdAt': creation_date,
            'updatedAt': datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat(),
            'workingDirectory': calc_path,
            'parameters': parameters,
            'results': results
        }
        
        return jsonify({
            'success': True,
            'data': {
                'calculation': calculation_instance,
                'files': {
                    'checkpoint_exists': file_manager.file_exists(calc_path, 'calculation.chk'),
                    'parameters_file_exists': parameters is not None,
                    'results_file_exists': results is not None,
                }
            }
        })
        
    except FileManagerError as e:
        logger.error(f"File manager error getting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Unable to access calculation data.'}), 500
    except OSError as e:
        logger.error(f"System error getting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'System error accessing calculation files.'}), 500
    except ValueError as e:
        logger.warning(f"Invalid calculation ID {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Invalid calculation ID format.'}), 400
    except Exception as e:
        logger.error(f"Unexpected error getting calculation details for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['PUT'])
@validate()
def update_calculation(calculation_id, body: CalculationUpdateRequest):
    """Update calculation metadata (currently only name)."""
    try:
        file_manager = CalculationFileManager()
        
        new_name = body.name
        
        # Update the display name (calculation_id remains the same)
        result_id = file_manager.rename_calculation(calculation_id, new_name)
        if not result_id:
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        logger.info(f"Updated display name for calculation {calculation_id} to '{new_name}'")
        return jsonify({
            'success': True,
            'data': {
                'message': 'Calculation renamed successfully.',
                'name': new_name
            }
        })

    except FileManagerError as e:
        logger.error(f"File manager error updating calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Unable to update calculation data.'}), 500
    except ValueError as e:
        logger.warning(f"Invalid input for calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Invalid input data.'}), 400
    except OSError as e:
        logger.error(f"System error updating calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'System error updating calculation.'}), 500
    except Exception as e:
        logger.error(f"Unexpected error updating calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>/cancel', methods=['POST'])
def cancel_calculation(calculation_id):
    """Cancel a running calculation."""
    try:
        process_manager = get_process_manager()
        
        # Check if calculation is running
        if not process_manager.is_running(calculation_id):
            return jsonify({
                'success': False, 
                'error': f'Calculation "{calculation_id}" is not currently running.'
            }), 400
        
        # Attempt to cancel the calculation
        cancelled = process_manager.cancel_calculation(calculation_id)
        
        if cancelled:
            logger.info(f"Successfully cancelled calculation {calculation_id}")
            return jsonify({
                'success': True,
                'data': {
                    'message': f'Calculation "{calculation_id}" has been cancelled successfully',
                    'calculation_id': calculation_id
                }
            })
        else:
            return jsonify({
                'success': False,
                'error': f'Failed to cancel calculation "{calculation_id}". It may have already completed.'
            }), 400
            
    except ProcessManagerError as e:
        logger.error(f"Process manager error cancelling calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Process manager unavailable.'}), 503
    except ValueError as e:
        logger.warning(f"Invalid calculation ID for cancellation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Invalid calculation ID.'}), 400
    except Exception as e:
        logger.error(f"Unexpected error cancelling calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['DELETE'])
def delete_calculation(calculation_id):
    """Delete a calculation and its files."""
    try:
        process_manager = get_process_manager()
        
        # Try to cancel the calculation if it's running
        if process_manager.is_running(calculation_id):
            logger.info(f"Cancelling running calculation {calculation_id} before deletion")
            process_manager.cancel_calculation(calculation_id)
            # Wait a moment for the cancellation to take effect
            time.sleep(0.5)
        
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)

        if not os.path.isdir(calc_path):
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404

        # Delete the calculation directory
        shutil.rmtree(calc_path)
        logger.info(f"Deleted calculation directory: {calc_path}")
        

        return jsonify({
            'success': True,
            'data': {
                'message': f'Calculation "{calculation_id}" has been deleted successfully',
                'deleted_id': calculation_id
            }
        })
        
    except ProcessManagerError as e:
        logger.error(f"Process manager error deleting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Process manager unavailable.'}), 503
    except FileManagerError as e:
        logger.error(f"File manager error deleting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Unable to access calculation files.'}), 500
    except PermissionError as e:
        logger.error(f"Permission error deleting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Permission denied deleting calculation.'}), 403
    except OSError as e:
        logger.error(f"System error deleting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'System error deleting calculation files.'}), 500
    except ValueError as e:
        logger.warning(f"Invalid calculation ID for deletion {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Invalid calculation ID.'}), 400
    except Exception as e:
        logger.error(f"Unexpected error deleting calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>/orbitals', methods=['GET'])
def get_orbitals(calculation_id):
    """Get molecular orbital information for a calculation."""
    try:
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        if not os.path.isdir(calc_path):
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        # Check if calculation is completed
        status = file_manager.read_calculation_status(calc_path)
        if status != 'completed':
            return jsonify({
                'success': False, 
                'error': f'Calculation "{calculation_id}" is not completed. Status: {status}'
            }), 400
        
        # Initialize orbital generator
        orbital_generator = MolecularOrbitalGenerator(calc_path)
        
        # Validate calculation data
        if not orbital_generator.validate_calculation():
            return jsonify({
                'success': False,
                'error': 'Orbital data is not available or calculation data is invalid.'
            }), 404
        
        # Get orbital summary
        orbital_summary = orbital_generator.get_orbital_summary()
        
        logger.info(f"Retrieved orbital information for calculation {calculation_id}")
        logger.info(f"Total orbitals: {orbital_summary['total_orbitals']}, "
                   f"HOMO: {orbital_summary['homo_index']}, LUMO: {orbital_summary['lumo_index']}")
        
        return jsonify({
            'success': True,
            'data': orbital_summary
        })
        
    except CalculationError as e:
        logger.error(f"Calculation error getting orbitals for {calculation_id}: {e}")
        return jsonify({'success': False, 'error': str(e)}), 500
    except FileManagerError as e:
        logger.error(f"File manager error getting orbitals for {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Unable to access calculation files.'}), 500
    except Exception as e:
        logger.error(f"Unexpected error getting orbitals for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>/orbitals/<int:orbital_index>/cube', methods=['GET'])
def get_orbital_cube(calculation_id, orbital_index):
    """Generate and return CUBE file for specific molecular orbital."""
    try:
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        if not os.path.isdir(calc_path):
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        # Check if calculation is completed
        status = file_manager.read_calculation_status(calc_path)
        if status != 'completed':
            return jsonify({
                'success': False, 
                'error': f'Calculation "{calculation_id}" is not completed. Status: {status}'
            }), 400
        
        # Get query parameters with defaults
        grid_size = request.args.get('gridSize', 80, type=int)
        isovalue_pos = request.args.get('isovaluePos', 0.02, type=float)
        isovalue_neg = request.args.get('isovalueNeg', -0.02, type=float)
        
        # Validate parameters
        if grid_size < 40 or grid_size > 120:
            return jsonify({'success': False, 'error': 'Grid size must be between 40 and 120.'}), 400
        if isovalue_pos < 0.001 or isovalue_pos > 0.1:
            return jsonify({'success': False, 'error': 'Positive isovalue must be between 0.001 and 0.1.'}), 400
        if isovalue_neg > -0.001 or isovalue_neg < -0.1:
            return jsonify({'success': False, 'error': 'Negative isovalue must be between -0.1 and -0.001.'}), 400
        
        # Initialize orbital generator
        orbital_generator = MolecularOrbitalGenerator(calc_path)
        
        # Validate calculation data
        if not orbital_generator.validate_calculation():
            return jsonify({
                'success': False,
                'error': 'Orbital data is not available or calculation data is invalid.'
            }), 404
        
        # Generate CUBE file
        logger.info(f"Generating CUBE file for calculation {calculation_id}, orbital {orbital_index}")
        logger.info(f"Parameters: grid_size={grid_size}, isovalue_pos={isovalue_pos}, isovalue_neg={isovalue_neg}")
        
        cube_data = orbital_generator.generate_cube_file(
            orbital_index=orbital_index,
            grid_size=grid_size,
            isovalue_pos=isovalue_pos,
            isovalue_neg=isovalue_neg,
            return_content=True,
            save_to_disk=True
        )
        
        if cube_data.get('cached', False):
            logger.info(f"Using cached CUBE file for calculation {calculation_id}, orbital {orbital_index}")
        else:
            logger.info(f"Successfully generated CUBE file for calculation {calculation_id}, orbital {orbital_index}")
        logger.info(f"File size: {cube_data['generation_params']['file_size_kb']:.1f} KB")
        
        return jsonify({
            'success': True,
            'data': cube_data
        })
        
    except CalculationError as e:
        logger.error(f"Calculation error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}")
        return jsonify({'success': False, 'error': str(e)}), 500
    except FileManagerError as e:
        logger.error(f"File manager error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}")
        return jsonify({'success': False, 'error': 'Unable to access calculation files.'}), 500
    except ValueError as e:
        logger.warning(f"Invalid orbital index for calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': 'Invalid orbital index.'}), 400
    except Exception as e:
        logger.error(f"Unexpected error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['GET'])
def list_cube_files(calculation_id):
    """List all CUBE files for a calculation."""
    try:
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        if not os.path.isdir(calc_path):
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        # Get CUBE file information
        cube_files = file_manager.get_cube_files_info(calc_path)
        
        logger.info(f"Found {len(cube_files)} CUBE files for calculation {calculation_id}")
        
        return jsonify({
            'success': True,
            'data': {
                'calculation_id': calculation_id,
                'cube_files': cube_files,
                'total_files': len(cube_files),
                'total_size_kb': sum(f['file_size_kb'] for f in cube_files)
            }
        })
        
    except Exception as e:
        logger.error(f"Unexpected error listing CUBE files for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['DELETE'])
def delete_cube_files(calculation_id):
    """Delete CUBE files for a calculation."""
    try:
        # Get query parameters
        orbital_index = request.args.get('orbital_index', type=int)
        
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        if not os.path.isdir(calc_path):
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        # Delete CUBE files
        deleted_count = file_manager.delete_cube_files(calc_path, orbital_index)
        
        if deleted_count > 0:
            if orbital_index is not None:
                logger.info(f"Deleted {deleted_count} CUBE files for orbital {orbital_index} in calculation {calculation_id}")
                message = f"Deleted {deleted_count} CUBE files for orbital {orbital_index}."
            else:
                logger.info(f"Deleted {deleted_count} CUBE files for calculation {calculation_id}")
                message = f"Deleted {deleted_count} CUBE files."
        else:
            if orbital_index is not None:
                message = f"No CUBE files found for orbital {orbital_index}."
            else:
                message = "No CUBE files found."
        
        return jsonify({
            'success': True,
            'data': {
                'calculation_id': calculation_id,
                'orbital_index': orbital_index,
                'deleted_files': deleted_count,
                'message': message
            }
        })
        
    except Exception as e:
        logger.error(f"Unexpected error deleting CUBE files for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@sock.route('/ws/calculations/<calculation_id>')
def calculation_status_socket(ws, calculation_id):
    """
    WebSocket endpoint for efficient, event-driven calculation monitoring.
    Uses file system event watching instead of polling for better performance.
    """
    # 一時的IDの場合は特別なログメッセージを出力
    if calculation_id.startswith('new-calculation-'):
        logger.info(f"WebSocket connection attempt for temporary calculation ID: {calculation_id}")
    else:
        logger.info(f"WebSocket connection established for calculation {calculation_id}")
    file_manager = CalculationFileManager()
    calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
    
    # State management for the connection
    connection_active = True
    
    # Calculation directory existence check
    if not os.path.isdir(calc_path):
        # 一時的IDの場合は特別なログメッセージを出力
        if calculation_id.startswith('new-calculation-'):
            logger.info(f"WebSocket connection attempted for temporary calculation ID: {calculation_id}. Closing connection.")
            error_message = f'Temporary calculation ID "{calculation_id}" does not exist on server.'
        else:
            logger.warning(f"Calculation directory not found: {calc_path}")
            error_message = f'Calculation "{calculation_id}" not found.'
        
        try:
            ws.send(json.dumps({
                'error': error_message,
                'id': calculation_id,
                'is_temporary': calculation_id.startswith('new-calculation-')
            }))
        except Exception as send_error:
            logger.error(f"Failed to send error message: {send_error}")
        finally:
            _close_websocket_safely(ws, calculation_id)
        return

    def build_calculation_instance(calc_id: str, calc_path: str, file_manager: CalculationFileManager) -> Dict:
        """Build a complete calculation instance from file system data."""
        try:
            parameters = file_manager.read_calculation_parameters(calc_path) or {}
            results = file_manager.read_calculation_results(calc_path)
            status = file_manager.read_calculation_status(calc_path) or 'unknown'
            display_name = file_manager._get_display_name(calc_id, parameters)
            
            # Safe date retrieval
            try:
                creation_date = parameters.get('created_at', 
                    datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat())
                updated_date = datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat()
            except OSError:
                current_time = datetime.now().isoformat()
                creation_date = current_time
                updated_date = current_time

            return {
                'id': calc_id,
                'name': display_name,
                'status': status,
                'createdAt': creation_date,
                'updatedAt': updated_date,
                'parameters': parameters,
                'results': results,
                'workingDirectory': calc_path,
            }
        except Exception as e:
            logger.error(f"Error building calculation instance for {calc_id}: {e}")
            raise

    def on_file_change(file_data: Dict):
        """Callback for file system events - sends updated data to WebSocket client."""
        nonlocal connection_active
        
        if not connection_active:
            return
        
        try:
            # Build complete calculation instance directly (no gevent.spawn needed)
            calculation_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
            
            # Send update to WebSocket client
            try:
                ws.send(json.dumps(calculation_instance))
            except Exception as send_error:
                logger.error(f"Failed to send WebSocket update for {calculation_id}: {send_error}")
                connection_active = False
                return
            
            # Close connection if calculation is finished
            if calculation_instance['status'] in ['completed', 'error']:
                logger.info(f"Calculation {calculation_id} finished with status '{calculation_instance['status']}'. Closing WebSocket after brief delay.")
                connection_active = False
                # Small delay to ensure client receives the final status update
                try:
                    gevent.sleep(0.1)
                    _close_websocket_safely(ws, calculation_id)
                except Exception as close_error:
                    logger.debug(f"Error during delayed WebSocket close for {calculation_id}: {close_error}")
                
        except Exception as e:
            logger.error(f"Error in file change callback for {calculation_id}: {e}")
            try:
                ws.send(json.dumps({
                    'error': 'Failed to read calculation data',
                    'id': calculation_id
                }))
            except:
                pass  # Ignore send errors in error handler
            connection_active = False

    try:
        # Register this WebSocket connection for immediate notifications
        with websocket_lock:
            if calculation_id not in active_websockets:
                active_websockets[calculation_id] = set()
            active_websockets[calculation_id].add(ws)
        
        # Initialize file watcher and register this connection
        watcher = get_websocket_watcher(file_manager.get_base_directory())
        watcher.add_connection(calculation_id, on_file_change)
        
        # Send initial state immediately
        try:
            initial_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
            ws.send(json.dumps(initial_instance))
            logger.info(f"Sent initial state for calculation {calculation_id} (status: {initial_instance['status']})")
            
            # If calculation is already finished, close after sending initial state
            if initial_instance['status'] in ['completed', 'error']:
                logger.info(f"Calculation {calculation_id} already finished. Closing WebSocket after initial state.")
                connection_active = False
                _close_websocket_safely(ws, calculation_id)
                return
                
        except Exception as e:
            logger.error(f"Error sending initial state for {calculation_id}: {e}")
            ws.send(json.dumps({
                'error': 'Failed to read initial calculation data',
                'id': calculation_id
            }))
            connection_active = False
            return
        
        # Keep connection alive until client disconnects or calculation completes
        # The file watcher will trigger on_file_change() when files are modified
        try:
            while connection_active:
                # Use gevent-compatible receive to detect client disconnections
                try:
                    message = ws.receive(timeout=30.0)  # 30 second timeout
                    if message is None:
                        # Client disconnected
                        logger.info(f"Client disconnected from calculation {calculation_id}")
                        break
                except gevent.Timeout:
                    # Expected timeout - continue the loop
                    continue
                except Exception as receive_error:
                    # Handle various disconnection scenarios
                    error_msg = str(receive_error).lower()
                    if any(keyword in error_msg for keyword in ['connection', 'closed', 'broken', 'reset']):
                        logger.info(f"Client disconnected from calculation {calculation_id}: {receive_error}")
                    else:
                        logger.debug(f"WebSocket receive error for {calculation_id}: {receive_error}")
                    break
                    
        except Exception as e:
            logger.error(f"Error in WebSocket receive loop for {calculation_id}: {e}")
            
    except WebSocketError as e:
        logger.error(f"WebSocket communication error for {calculation_id}: {e}")
    except FileManagerError as e:
        logger.error(f"File access error in WebSocket for {calculation_id}: {e}")
    except OSError as e:
        logger.error(f"System error in WebSocket for {calculation_id}: {e}")
    except Exception as e:
        logger.error(f"Unexpected WebSocket error for {calculation_id}: {e}", exc_info=True)
    finally:
        # Always clean up the connection
        connection_active = False
        
        # Remove from immediate notification registry
        try:
            with websocket_lock:
                if calculation_id in active_websockets:
                    active_websockets[calculation_id].discard(ws)
                    if not active_websockets[calculation_id]:
                        del active_websockets[calculation_id]
        except Exception as registry_error:
            logger.debug(f"Error removing WebSocket from registry for {calculation_id}: {registry_error}")
        
        # Remove from file watcher
        try:
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.remove_connection(calculation_id, on_file_change)
        except (AttributeError, RuntimeError) as cleanup_error:
            # These are expected errors during shutdown or when observer is not available
            logger.debug(f"Observer cleanup issue for {calculation_id}: {cleanup_error}")
        except Exception as cleanup_error:
            logger.error(f"Unexpected error cleaning up file watcher for {calculation_id}: {cleanup_error}")
        
        _close_websocket_safely(ws, calculation_id)


def _close_websocket_safely(ws, calculation_id: str):
    """Safely close a WebSocket connection with proper logging."""
    try:
        # Check if WebSocket is still open before attempting to close
        if hasattr(ws, 'closed') and not ws.closed:
            ws.close(code=1000, reason="Normal closure")  # Use proper close code
        elif not hasattr(ws, 'closed'):
            # Fallback for different WebSocket implementations
            ws.close()
        logger.info(f"WebSocket connection closed for calculation {calculation_id}")
    except Exception as close_error:
        # Log different types of close errors appropriately
        error_msg = str(close_error).lower()
        if any(keyword in error_msg for keyword in ['closed', 'broken', 'connection']):
            logger.debug(f"WebSocket already closed for {calculation_id}: {close_error}")
        else:
            logger.warning(f"Error closing WebSocket for {calculation_id}: {close_error}")


@app.errorhandler(ValidationError)
def validation_error_handler(error):
    """Handle Pydantic validation errors with consistent error format."""
    errors = error.errors()
    error_messages = []
    for err in errors:
        field = '.'.join(str(x) for x in err['loc'])
        message = err['msg']
        error_messages.append(f"{field}: {message}")
    
    combined_message = "Validation failed: " + "; ".join(error_messages)
    logger.warning(f"Validation error: {combined_message}")
    return jsonify({'success': False, 'error': combined_message}), 400

@app.errorhandler(400)
def bad_request_handler(error):
    """Handle bad requests, filtering out WebSocket close frame errors."""
    # Check if this is a WebSocket close frame being misinterpreted as HTTP
    error_description = str(error).lower()
    if any(indicator in error_description for indicator in [
        'invalid http method', 'expected get method', 'x88x82', '\\x88\\x82'
    ]):
        # This is likely a WebSocket close frame, log as debug instead of error
        logger.debug(f"WebSocket close frame misinterpreted as HTTP request: {error}")
        return '', 400  # Return empty response for WebSocket frames
    
    # For genuine bad requests, return proper error response
    return jsonify({'success': False, 'error': 'Bad request.'}), 400

@app.errorhandler(404)
def not_found(error):
    return jsonify({'success': False, 'error': 'Endpoint not found.'}), 404

@app.errorhandler(405)
def method_not_allowed(error):
    return jsonify({'success': False, 'error': 'Method not allowed for this endpoint.'}), 405

if __name__ == '__main__':
    host = os.environ.get('FLASK_RUN_HOST', '127.0.0.1')
    port_arg = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    debug = os.environ.get('FLASK_ENV') == 'development'

    # GeventのWSGIServerはポート0を自動割り当てしないため、手動で空きポートを探す
    if port_arg == 0:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind((host, 0))
            actual_port = s.getsockname()[1]
    else:
        actual_port = port_arg
    
    # Electronプロセスにポート番号を通知
    print(f"FLASK_SERVER_PORT:{actual_port}", file=sys.stdout, flush=True)
    
    logger.info(f"Starting API server with gevent on http://{host}:{actual_port} (Debug: {debug})")
    
    # ===== ここから修正 =====
    # handler_class の指定を削除し、gevent と Flask-Sock の標準的な連携に任せる
    http_server = WSGIServer((host, actual_port), app)
    # ===== ここまで修正 =====
    try:
        http_server.serve_forever()
    except KeyboardInterrupt:
        logger.info("Server stopped by user.")
    finally:
        logger.info("Shutting down server...")
        http_server.stop()
        cleanup_resources()