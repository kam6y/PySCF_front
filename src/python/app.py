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
from flask_pydantic import validate
from pydantic import ValidationError
import socket
import shutil
from typing import Dict, Optional
from flask_socketio import SocketIO, emit, disconnect

from pubchem.client import PubChemClient, PubChemError, PubChemNotFoundError
from pubchem import parser as xyz_parser
from SMILES.smiles_converter import smiles_to_xyz, SMILESError
from quantum_calc import DFTCalculator, HFCalculator, MP2Calculator, CCSDCalculator, TDDFTCalculator, MolecularOrbitalGenerator, CalculationError, ConvergenceError, InputError, GeometryError, ProcessManagerError
from quantum_calc.ir_spectrum import create_ir_spectrum_from_calculation_results
from quantum_calc.exceptions import XYZValidationError, FileManagerError, ProcessManagerError, WebSocketError
from quantum_calc.file_manager import CalculationFileManager
from quantum_calc import get_process_manager, shutdown_process_manager, get_websocket_watcher, shutdown_websocket_watcher, get_all_supported_parameters, get_current_settings, update_app_settings
from quantum_calc.resource_manager import get_resource_manager
from generated_models import (
    PubChemSearchRequest, SMILESConvertRequest, XYZValidateRequest,
    QuantumCalculationRequest, CalculationUpdateRequest, AppSettings, SettingsUpdateRequest,
    SystemResourceResponse, SystemResourceSummary, SystemResourceInfo, ResourceConstraints, AllocatedResources
)

# Load server configuration
def load_server_config():
    """Load server configuration from JSON file."""
    try:
        # Try to load from config directory first
        config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config', 'server-config.json')
        if not os.path.exists(config_path):
            # Fallback to bundled config for packaged applications
            config_path = os.path.join(os.path.dirname(__file__), 'config', 'server-config.json')
        
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Server configuration file not found at {config_path}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = json.load(f)
        
        print(f"Loaded server configuration from: {config_path}")
        return config
    except Exception as e:
        print(f"WARNING: Failed to load server configuration: {e}. Using defaults.")
        # Return default configuration
        return {
            "server": {"host": "127.0.0.1", "port": {"default": 5000, "auto_detect": True}},
            "gunicorn": {"workers": 1, "threads": 4, "timeout": 0, "worker_class": "sync"},
            "socketio": {"cors_allowed_origins": "*", "async_mode": "threading"},
            "development": {"debug": False},
            "production": {"use_gunicorn": True},
            "logging": {"level": "INFO"}
        }

# Global configuration
SERVER_CONFIG = load_server_config()

# Configure logging based on configuration
log_level = getattr(logging, SERVER_CONFIG.get('logging', {}).get('level', 'INFO').upper())
log_format = SERVER_CONFIG.get('logging', {}).get('format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logging.basicConfig(
    level=log_level,
    format=log_format
)
logger = logging.getLogger(__name__)

# Initialize Flask app and extensions
app = Flask(__name__)
CORS(app)  # Enable CORS for cross-origin requests

# Initialize SocketIO with configuration-based settings
socketio_config = SERVER_CONFIG.get('socketio', {})
socketio = SocketIO(
    app, 
    cors_allowed_origins=socketio_config.get('cors_allowed_origins', "*"),
    async_mode=socketio_config.get('async_mode', 'threading'),
    ping_timeout=socketio_config.get('ping_timeout', 60),
    ping_interval=socketio_config.get('ping_interval', 25)
)

# Initialize PubChem client
pubchem_client = PubChemClient(timeout=30)

# Global WebSocket connection registry for immediate notifications  
# calculation_id -> set of session IDs
active_websockets = {}
websocket_lock = threading.Lock()


def send_immediate_websocket_notification(calculation_id: str, status: str, error_message: str = None):
    """Send immediate WebSocket notification to all connected clients for a calculation."""
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
            
            # Add error if provided (両方のフィールドに設定してフロントエンドとの整合性を確保)
            if error_message:
                calculation_instance['error'] = error_message
                calculation_instance['errorMessage'] = error_message
            
            # Send to all clients in the calculation room
            socketio.emit('calculation_update', calculation_instance, room=f'calculation_{calculation_id}')
            
            # Also send to global updates room for non-active calculations monitoring
            socketio.emit('calculation_update', calculation_instance, room='global_updates')
            logger.debug(f"Sent immediate notification for calculation {calculation_id} with status {status} to both specific and global rooms")
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


@app.route('/api/quantum/supported-parameters', methods=['GET'])
def get_supported_parameters():
    """Get supported quantum chemistry parameters including basis functions, exchange-correlation functionals, and solvents."""
    try:
        logger.info("Getting supported quantum chemistry parameters")
        
        # Get all supported parameters from our quantum_calc module
        parameters = get_all_supported_parameters()
        
        logger.info("Successfully retrieved supported parameters")
        return jsonify({
            'success': True,
            'data': parameters
        }), 200
        
    except Exception as e:
        logger.error(f"Error getting supported parameters: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve supported parameters: {str(e)}'
        }), 500


@app.route('/api/settings', methods=['GET'])
def get_settings():
    """Get current application settings."""
    try:
        logger.info("Getting application settings")
        
        # Get current settings
        settings = get_current_settings()
        
        logger.info(f"Successfully retrieved settings: {settings}")
        return jsonify({
            'success': True,
            'data': {
                'settings': settings.model_dump()
            }
        })
    except Exception as e:
        logger.error(f"Failed to retrieve settings: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve settings: {str(e)}'
        }), 500


@app.route('/api/settings', methods=['PUT'])
@validate()
def update_settings(body: SettingsUpdateRequest):
    """Update application settings."""
    try:
        logger.info(f"Updating application settings: {body}")
        
        # Extract settings from root model
        new_settings = body.root if hasattr(body, 'root') else body
        
        # Update settings
        updated_settings = update_app_settings(new_settings.model_dump())
        
        # Update process manager with new parallel instance limit
        process_manager = get_process_manager()
        process_manager.set_max_parallel_instances(updated_settings.max_parallel_instances)
        
        # Update resource manager with new resource constraints
        resource_manager = get_resource_manager()
        resource_manager.update_resource_constraints(
            max_cpu_utilization_percent=updated_settings.max_cpu_utilization_percent,
            max_memory_utilization_percent=updated_settings.max_memory_utilization_percent
        )
        
        logger.info(f"Successfully updated settings: {updated_settings}")
        return jsonify({
            'success': True,
            'data': {
                'settings': updated_settings.model_dump()
            }
        })
    except Exception as e:
        logger.error(f"Failed to update settings: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to update settings: {str(e)}'
        }), 400 if isinstance(e, ValueError) else 500


@app.route('/api/system/resource-status', methods=['GET'])
def get_system_resource_status():
    """Get current system resource status including constraints and allocation."""
    try:
        logger.info("Getting system resource status")
        
        # Get resource manager
        resource_manager = get_resource_manager()
        
        # Get resource summary
        resource_summary = resource_manager.get_resource_summary()
        
        # Create response using Pydantic models
        system_info = SystemResourceInfo(
            total_cpu_cores=resource_summary['system_info']['total_cpu_cores'],
            total_memory_mb=resource_summary['system_info']['total_memory_mb'],
            available_memory_mb=resource_summary['system_info']['available_memory_mb'],
            cpu_usage_percent=resource_summary['system_info']['cpu_usage_percent'],
            memory_usage_percent=resource_summary['system_info']['memory_usage_percent'],
            timestamp=datetime.fromisoformat(resource_summary['system_info']['timestamp'].replace('Z', '+00:00'))
        )
        
        resource_constraints = ResourceConstraints(
            max_cpu_utilization_percent=resource_summary['resource_constraints']['max_cpu_utilization_percent'],
            max_memory_utilization_percent=resource_summary['resource_constraints']['max_memory_utilization_percent'],
            max_allowed_cpu_cores=resource_summary['resource_constraints']['max_allowed_cpu_cores'],
            max_allowed_memory_mb=resource_summary['resource_constraints']['max_allowed_memory_mb']
        )
        
        allocated_resources = AllocatedResources(
            total_allocated_cpu_cores=resource_summary['allocated_resources']['total_allocated_cpu_cores'],
            total_allocated_memory_mb=resource_summary['allocated_resources']['total_allocated_memory_mb'],
            available_cpu_cores=resource_summary['allocated_resources']['available_cpu_cores'],
            available_memory_mb=resource_summary['allocated_resources']['available_memory_mb'],
            active_calculations_count=resource_summary['allocated_resources']['active_calculations_count']
        )
        
        summary = SystemResourceSummary(
            system_info=system_info,
            resource_constraints=resource_constraints,
            allocated_resources=allocated_resources
        )
        
        logger.info(f"Successfully retrieved system resource status: {summary}")
        return jsonify({
            'success': True,
            'data': summary.model_dump()
        })
        
    except Exception as e:
        logger.error(f"Failed to retrieve system resource status: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve system resource status: {str(e)}'
        }), 500


def validate_calculation_parameters(body: QuantumCalculationRequest) -> Optional[str]:
    """
    Validate calculation parameters for compatibility and theoretical correctness.

    Returns:
        None if validation passes, error message string if validation fails
    """
    def get_enum_value(field_value):
        if hasattr(field_value, 'value'):
            return field_value.value
        return field_value

    calculation_method = get_enum_value(body.calculation_method)

    # Check HF method specific constraints
    if calculation_method == 'HF':
        # HF method should not use exchange-correlation functionals
        if body.exchange_correlation and body.exchange_correlation != 'B3LYP':
            logger.warning(f"HF calculation with non-default exchange_correlation: {body.exchange_correlation}")

        # TDDFT parameters are not applicable to HF
        if body.tddft_nstates and body.tddft_nstates > 10:
            logger.warning(f"TDDFT parameters specified for HF calculation - these will be ignored")

    # Check DFT method constraints
    if calculation_method == 'DFT':
        if not body.exchange_correlation:
            return "DFT method requires an exchange-correlation functional to be specified"

    # Check TDDFT method constraints
    if calculation_method == 'TDDFT':
        if not body.exchange_correlation:
            return "TDDFT method requires an exchange-correlation functional to be specified"
        if not body.tddft_nstates or body.tddft_nstates < 1:
            return "TDDFT method requires tddft_nstates to be specified and greater than 0"

    # Check CASCI/CASSCF method constraints
    if calculation_method in ['CASCI', 'CASSCF']:
        if not body.ncas or body.ncas < 1:
            return f"{calculation_method} method requires ncas (active space orbitals) to be specified and greater than 0"
        if not body.nelecas or body.nelecas < 1:
            return f"{calculation_method} method requires nelecas (active space electrons) to be specified and greater than 0"
        if body.nelecas > 2 * body.ncas:
            return f"{calculation_method} method: nelecas ({body.nelecas}) cannot exceed 2 * ncas ({2 * body.ncas})"

    # Check general parameter constraints
    if body.spin < 0:
        return "Spin multiplicity (2S) cannot be negative"

    if body.charges is not None and abs(body.charges) > 10:
        logger.warning(f"High molecular charge ({body.charges}) - please verify this is correct")

    return None


@app.route('/api/quantum/calculate', methods=['POST'])
@validate()
def quantum_calculate(body: QuantumCalculationRequest):
    """
    Starts a quantum chemistry calculation in the background.
    Immediately returns a calculation ID to track the job.
    """
    try:
        # Validate calculation parameters for compatibility and theoretical correctness
        validation_error = validate_calculation_parameters(body)
        if validation_error:
            logger.warning(f"Parameter validation failed: {validation_error}")
            return jsonify({
                'success': False,
                'error': f'Invalid parameters: {validation_error}'
            }), 400
        # Helper function to safely extract value from enum or string
        def get_enum_value(field_value):
            if hasattr(field_value, 'value'):
                return field_value.value
            return field_value
        
        # Prepare parameters using validated data from Pydantic model
        calculation_method = get_enum_value(body.calculation_method)

        parameters = {
            'calculation_method': calculation_method,
            'basis_function': body.basis_function,
            'charges': body.charges,
            'spin': body.spin,
            'solvent_method': get_enum_value(body.solvent_method),
            'solvent': body.solvent,
            'xyz': body.xyz,
            'name': body.name,
            'cpu_cores': body.cpu_cores,
            'memory_mb': body.memory_mb,
            'created_at': datetime.now().isoformat(),
        }

        # Add exchange_correlation parameter only for DFT methods
        # HF method does not use exchange-correlation functionals
        if calculation_method != 'HF':
            parameters['exchange_correlation'] = body.exchange_correlation
        else:
            # For HF method, explicitly set to None or omit to avoid confusion
            parameters['exchange_correlation'] = None
            logger.info(f"HF calculation - exchange_correlation parameter ignored (was: {body.exchange_correlation})")

        # Always include TDDFT parameters for consistency
        parameters.update({
            'tddft_nstates': body.tddft_nstates,
            'tddft_method': get_enum_value(body.tddft_method) if body.tddft_method else 'TDDFT',
            'tddft_analyze_nto': body.tddft_analyze_nto,
        })

        # Always include CASSCF・CASCI parameters for consistency
        parameters.update({
            'ncas': body.ncas,
            'nelecas': body.nelecas,
            'max_cycle_macro': body.max_cycle_macro,
            'max_cycle_micro': body.max_cycle_micro,
            'natorb': body.natorb,
            'conv_tol': body.conv_tol,
            'conv_tol_grad': body.conv_tol_grad
        })
        
        # Initialize file manager and create directory with error handling
        try:
            file_manager = CalculationFileManager()
            calc_dir = file_manager.create_calculation_dir(parameters['name'])
            calculation_id = os.path.basename(calc_dir)
            
            # Save initial parameters (status will be set based on submission result)
            file_manager.save_calculation_parameters(calc_dir, parameters)
            logger.info(f"Created calculation directory and saved parameters for calculation {calculation_id}")
            
        except Exception as file_error:
            logger.error(f"Failed to set up calculation files: {file_error}")
            return jsonify({'success': False, 'error': f'Failed to initialize calculation: {str(file_error)}'}), 500

        # Get process manager with enhanced error handling
        try:
            process_manager = get_process_manager()
            
            # Apply current settings to process manager if not already done
            try:
                from quantum_calc import update_process_manager_settings
                update_process_manager_settings()
            except Exception as settings_error:
                logger.warning(f"Failed to update process manager settings: {settings_error}")
                # Continue without settings update
                
        except Exception as pm_error:
            logger.error(f"Failed to initialize process manager: {pm_error}")
            # Clean up created directory on process manager failure
            try:
                import shutil
                shutil.rmtree(calc_dir, ignore_errors=True)
            except Exception:
                pass
            return jsonify({'success': False, 'error': f'System initialization error: Unable to initialize calculation system. Please check system resources and try again.'}), 503

        # Submit calculation to process pool or queue with enhanced error handling
        try:
            success, initial_status, waiting_reason = process_manager.submit_calculation(calculation_id, parameters)
        except Exception as submit_error:
            logger.error(f"Unexpected error during calculation submission: {submit_error}")
            # Update status to error and save error information
            file_manager.save_calculation_status(calc_dir, 'error')
            error_message = f'Failed to submit calculation: {str(submit_error)}'
            file_manager.save_calculation_results(calc_dir, {'error': error_message})
            
            # Send immediate WebSocket notification for error status
            send_immediate_websocket_notification(calculation_id, 'error', error_message)
            
            # Return error instance
            error_instance = {
                'id': calculation_id,
                'name': parameters['name'],
                'status': 'error',
                'createdAt': parameters['created_at'],
                'updatedAt': parameters['created_at'],
                'parameters': parameters,
                'error': error_message
            }
            return jsonify(error_instance), 500
        
        if not success:
            # If submission failed, update status to error
            error_message = waiting_reason if waiting_reason else 'Failed to submit calculation to process pool.'
            file_manager.save_calculation_status(calc_dir, 'error')
            file_manager.save_calculation_results(calc_dir, {'error': error_message})
            
            # Send immediate WebSocket notification for error status
            send_immediate_websocket_notification(calculation_id, 'error', error_message)
            
            logger.error(f"Failed to submit calculation {calculation_id} to process pool: {error_message}")
            
            # Create error instance to return
            error_instance = {
                'id': calculation_id,
                'name': parameters['name'],
                'status': 'error',
                'createdAt': parameters['created_at'],
                'updatedAt': parameters['created_at'],
                'parameters': parameters,
                'error': error_message
            }
            
            # Return error instance with 202 status (calculation was accepted but failed due to resources)
            return jsonify(error_instance), 202

        # Set initial status based on submission result
        file_manager.save_calculation_status(calc_dir, initial_status, waiting_reason)
        
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
            'status': initial_status,  # Use actual status (running or waiting)
            'createdAt': parameters['created_at'],
            'updatedAt': parameters['created_at'],
            'parameters': parameters
        }
        
        # Include waiting reason if available
        if waiting_reason is not None:
            initial_instance['waitingReason'] = waiting_reason

        return jsonify({'success': True, 'data': {'calculation': initial_instance}}), 202

    except (InputError, GeometryError) as e:
        logger.warning(f"Invalid calculation parameters: {e}")
        return jsonify({'success': False, 'error': f'Invalid input parameters: {str(e)}'}), 400
    except FileManagerError as e:
        logger.error(f"File management error during calculation setup: {e}")
        return jsonify({'success': False, 'error': 'Failed to set up calculation files.'}), 500
    except ProcessManagerError as e:
        logger.error(f"Process manager initialization/operation error: {e}")
        # Try to clean up any created files
        try:
            if 'calc_dir' in locals():
                import shutil
                shutil.rmtree(calc_dir, ignore_errors=True)
        except Exception:
            pass
        return jsonify({'success': False, 'error': f'Calculation system error: {str(e)}. The system may need to be restarted.'}), 503
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


@app.route('/api/debug/system-diagnostics', methods=['GET'])
def get_system_diagnostics():
    """Get comprehensive system diagnostics for troubleshooting."""
    try:
        logger.info("Getting comprehensive system diagnostics")
        
        # Collect system information
        diagnostics = {
            'timestamp': datetime.now().isoformat(),
            'service_info': {
                'service': 'pyscf-front-api',
                'version': '0.3.0',
                'pid': os.getpid(),
                'working_directory': os.getcwd()
            },
            'system_info': {
                'cpu_count': multiprocessing.cpu_count(),
                'platform': os.name,
                'python_version': sys.version
            }
        }
        
        # Process manager diagnostics
        try:
            process_manager = get_process_manager()
            diagnostics['process_manager'] = {
                'status': 'available',
                'max_workers': process_manager.max_workers,
                'max_parallel_instances': process_manager.max_parallel_instances,
                'active_calculations': len(process_manager.active_futures),
                'queued_calculations': len(process_manager.calculation_queue),
                'is_shutdown': process_manager._shutdown,
                'queue_status': process_manager.get_queue_status()
            }
        except Exception as pm_error:
            diagnostics['process_manager'] = {
                'status': 'error',
                'error': str(pm_error),
                'error_type': type(pm_error).__name__
            }
        
        # Resource manager diagnostics
        try:
            resource_manager = get_resource_manager()
            diagnostics['resource_manager'] = resource_manager.get_diagnostics()
        except Exception as rm_error:
            diagnostics['resource_manager'] = {
                'status': 'error',
                'error': str(rm_error),
                'error_type': type(rm_error).__name__
            }
        
        # File manager diagnostics
        try:
            from quantum_calc.file_manager import CalculationFileManager
            file_manager = CalculationFileManager()
            base_dir = file_manager.get_base_directory()
            
            diagnostics['file_manager'] = {
                'status': 'available',
                'base_directory': base_dir,
                'base_directory_exists': os.path.exists(base_dir),
                'base_directory_writable': os.access(base_dir, os.W_OK) if os.path.exists(base_dir) else False
            }
            
            # Count calculation directories
            try:
                calculations = file_manager.list_calculations()
                diagnostics['file_manager']['total_calculations'] = len(calculations)
                diagnostics['file_manager']['calculation_statuses'] = {}
                
                # Count by status
                status_counts = {}
                for calc in calculations[:20]:  # Limit to first 20 for performance
                    try:
                        status = file_manager.read_calculation_status(os.path.join(base_dir, calc['id']))
                        status_counts[status] = status_counts.get(status, 0) + 1
                    except Exception:
                        status_counts['unknown'] = status_counts.get('unknown', 0) + 1
                
                diagnostics['file_manager']['calculation_statuses'] = status_counts
            except Exception as calc_error:
                diagnostics['file_manager']['calculations_error'] = str(calc_error)
                
        except Exception as fm_error:
            diagnostics['file_manager'] = {
                'status': 'error',
                'error': str(fm_error),
                'error_type': type(fm_error).__name__
            }
        
        # Settings diagnostics
        try:
            settings = get_current_settings()
            diagnostics['settings'] = {
                'status': 'available',
                'settings': settings.model_dump()
            }
        except Exception as settings_error:
            diagnostics['settings'] = {
                'status': 'error',
                'error': str(settings_error),
                'error_type': type(settings_error).__name__
            }
        
        logger.info("Successfully retrieved comprehensive system diagnostics")
        return jsonify({
            'success': True,
            'data': diagnostics
        })
        
    except Exception as e:
        logger.error(f"Failed to retrieve system diagnostics: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve system diagnostics: {str(e)}'
        }), 500


@app.route('/api/debug/process-manager-diagnostics', methods=['GET'])
def get_process_manager_diagnostics():
    """Get detailed process manager diagnostics."""
    try:
        logger.info("Getting process manager diagnostics")
        
        try:
            process_manager = get_process_manager()
            
            # Get detailed process manager state
            diagnostics = {
                'timestamp': datetime.now().isoformat(),
                'status': 'available',
                'configuration': {
                    'max_workers': process_manager.max_workers,
                    'max_parallel_instances': process_manager.max_parallel_instances,
                    'is_shutdown': process_manager._shutdown
                },
                'current_state': {
                    'active_futures_count': len(process_manager.active_futures),
                    'active_calculation_ids': list(process_manager.active_futures.keys()),
                    'queued_calculations_count': len(process_manager.calculation_queue),
                    'completion_callbacks_count': len(process_manager.completion_callbacks)
                },
                'queue_details': [],
                'resource_monitoring': {
                    'monitoring_active': process_manager._resource_monitor_thread is not None and process_manager._resource_monitor_thread.is_alive(),
                    'monitoring_interval': process_manager._resource_monitor_interval
                }
            }
            
            # Get detailed queue information
            for i, queued_calc in enumerate(process_manager.calculation_queue[:10]):  # Limit to first 10
                queue_item = {
                    'position': i + 1,
                    'calculation_id': queued_calc.calculation_id,
                    'created_at': queued_calc.created_at.isoformat(),
                    'waiting_reason': queued_calc.waiting_reason,
                    'calculation_method': queued_calc.parameters.get('calculation_method', 'unknown'),
                    'cpu_cores': queued_calc.parameters.get('cpu_cores', 'unknown'),
                    'memory_mb': queued_calc.parameters.get('memory_mb', 'unknown')
                }
                diagnostics['queue_details'].append(queue_item)
            
            # Executor status
            if process_manager.executor is not None:
                diagnostics['executor'] = {
                    'available': True,
                    'type': type(process_manager.executor).__name__
                }
            else:
                diagnostics['executor'] = {
                    'available': False,
                    'error': 'ProcessPoolExecutor is None'
                }
                
        except Exception as pm_error:
            diagnostics = {
                'timestamp': datetime.now().isoformat(),
                'status': 'error',
                'error': str(pm_error),
                'error_type': type(pm_error).__name__
            }
        
        logger.info("Successfully retrieved process manager diagnostics")
        return jsonify({
            'success': True,
            'data': diagnostics
        })
        
    except Exception as e:
        logger.error(f"Failed to retrieve process manager diagnostics: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve process manager diagnostics: {str(e)}'
        }), 500


@app.route('/api/debug/resource-manager-diagnostics', methods=['GET'])
def get_resource_manager_diagnostics():
    """Get detailed resource manager diagnostics."""
    try:
        logger.info("Getting resource manager diagnostics")
        
        try:
            resource_manager = get_resource_manager()
            diagnostics = resource_manager.get_diagnostics()
            diagnostics['timestamp'] = datetime.now().isoformat()
            
        except Exception as rm_error:
            diagnostics = {
                'timestamp': datetime.now().isoformat(),
                'status': 'error',
                'error': str(rm_error),
                'error_type': type(rm_error).__name__
            }
        
        logger.info("Successfully retrieved resource manager diagnostics")
        return jsonify({
            'success': True,
            'data': diagnostics
        })
        
    except Exception as e:
        logger.error(f"Failed to retrieve resource manager diagnostics: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'Failed to retrieve resource manager diagnostics: {str(e)}'
        }), 500


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
        status, waiting_reason = file_manager.read_calculation_status_details(calc_path)
        
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
        
        # Include waiting reason if available
        if waiting_reason is not None:
            calculation_instance['waitingReason'] = waiting_reason
        
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


@app.route('/api/quantum/calculations/<calculation_id>/ir-spectrum', methods=['GET'])
def get_ir_spectrum(calculation_id):
    """Generate and return IR spectrum for a calculation."""
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
                'error': f'Calculation is not completed (status: {status}). IR spectrum cannot be generated.'
            }), 400
        
        # Read calculation results
        results = file_manager.read_calculation_results(calc_path)
        if not results:
            return jsonify({
                'success': False,
                'error': 'Calculation results not found.'
            }), 404
        
        # Check for vibrational frequency data
        if not results.get('vibrational_frequencies'):
            return jsonify({
                'success': False,
                'error': 'No vibrational frequency data found. Frequency analysis may not have been performed or failed.'
            }), 400
        
        # Get query parameters for spectrum customization
        broadening_fwhm = request.args.get('broadening_fwhm', default=100.0, type=float)
        x_min = request.args.get('x_min', default=400.0, type=float)
        x_max = request.args.get('x_max', default=4000.0, type=float)
        show_peaks = request.args.get('show_peaks', default=True, type=bool)
        
        # Validate parameters
        if broadening_fwhm <= 0:
            return jsonify({
                'success': False,
                'error': 'Broadening FWHM must be positive.'
            }), 400
        
        if x_min >= x_max:
            return jsonify({
                'success': False,
                'error': 'x_min must be less than x_max.'
            }), 400
        
        logger.info(f"Generating IR spectrum for calculation {calculation_id}")
        logger.info(f"Parameters: FWHM={broadening_fwhm}, range=({x_min}, {x_max}), show_peaks={show_peaks}")
        
        # Generate IR spectrum
        ir_result = create_ir_spectrum_from_calculation_results(
            results,
            broadening_fwhm=broadening_fwhm,
            x_range=(x_min, x_max)
        )
        
        if not ir_result.get('success'):
            error_msg = ir_result.get('error', 'Unknown error occurred during IR spectrum generation')
            logger.error(f"IR spectrum generation failed for {calculation_id}: {error_msg}")
            return jsonify({
                'success': False,
                'error': f'IR spectrum generation failed: {error_msg}'
            }), 500
        
        # Extract spectrum data and plot
        spectrum_data = ir_result.get('spectrum_data', {})
        plot_image = ir_result.get('plot_image_base64')
        
        logger.info(f"IR spectrum generated successfully for {calculation_id}")
        logger.info(f"Spectrum contains {len(spectrum_data.get('peaks', []))} peaks")
        
        return jsonify({
            'success': True,
            'data': {
                'calculation_id': calculation_id,
                'spectrum': {
                    'x_axis': spectrum_data.get('x_axis', []),
                    'y_axis': spectrum_data.get('spectrum', []),
                    'peaks': spectrum_data.get('peaks', []),
                    'metadata': spectrum_data.get('metadata', {})
                },
                'plot_image_base64': plot_image,
                'generation_info': {
                    'broadening_fwhm_cm': broadening_fwhm,
                    'frequency_range_cm': [x_min, x_max],
                    'peaks_marked': show_peaks,
                    'generated_at': datetime.now().isoformat()
                }
            }
        })
        
    except ValueError as e:
        logger.error(f"Invalid parameters for IR spectrum generation ({calculation_id}): {e}")
        return jsonify({'success': False, 'error': f'Invalid parameters: {str(e)}'}), 400
    except Exception as e:
        logger.error(f"Unexpected error generating IR spectrum for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


# Flask-SocketIO WebSocket handlers for calculation monitoring
def build_calculation_instance(calc_id: str, calc_path: str, file_manager: CalculationFileManager) -> Dict:
    """Build a complete calculation instance from file system data."""
    try:
        parameters = file_manager.read_calculation_parameters(calc_path) or {}
        results = file_manager.read_calculation_results(calc_path)
        status = file_manager.read_calculation_status(calc_path) or 'pending'
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


@socketio.on('join_calculation')
def on_join_calculation(data):
    """Join a calculation room to receive real-time updates."""
    calculation_id = data.get('calculation_id')
    if not calculation_id:
        emit('error', {'error': 'calculation_id is required'})
        return
    
    from flask_socketio import join_room
    from flask import session
    
    # 一時的IDの場合は特別なログメッセージを出力
    if calculation_id.startswith('new-calculation-'):
        logger.info(f"SocketIO connection attempt for temporary calculation ID: {calculation_id}")
    else:
        logger.info(f"SocketIO connection established for calculation {calculation_id}")
    
    file_manager = CalculationFileManager()
    calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
    
    # Calculation directory existence check
    if not os.path.isdir(calc_path):
        # 一時的IDの場合は特別なログメッセージを出力
        if calculation_id.startswith('new-calculation-'):
            logger.info(f"SocketIO connection attempted for temporary calculation ID: {calculation_id}. Sending error.")
            error_message = f'Temporary calculation ID "{calculation_id}" does not exist on server.'
        else:
            logger.warning(f"Calculation directory not found: {calc_path}")
            error_message = f'Calculation "{calculation_id}" not found.'
        
        emit('error', {
            'error': error_message,
            'id': calculation_id,
            'is_temporary': calculation_id.startswith('new-calculation-')
        })
        return

    # Join the calculation room
    room = f'calculation_{calculation_id}'
    join_room(room)
    
    # Store calculation_id in session for cleanup
    session['calculation_id'] = calculation_id
    
    # Define file change callback for this connection
    def on_file_change(file_data: Dict):
        """Callback for file system events - sends updated data to SocketIO client."""
        try:
            # Build complete calculation instance
            calculation_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
            
            # Send update to all clients in the room
            socketio.emit('calculation_update', calculation_instance, room=room)
            
            # Disconnect clients if calculation is finished
            if calculation_instance['status'] in ['completed', 'error']:
                logger.info(f"Calculation {calculation_id} finished with status '{calculation_instance['status']}'.")
                # Note: Don't force disconnect - let client handle completion
                
        except Exception as e:
            logger.error(f"Error in file change callback for {calculation_id}: {e}")
            socketio.emit('error', {
                'error': 'Failed to read calculation data',
                'id': calculation_id
            }, room=room)

    try:
        # Initialize file watcher and register this connection
        watcher = get_websocket_watcher(file_manager.get_base_directory())
        watcher.add_connection(calculation_id, on_file_change)
        
        # Send initial state immediately
        initial_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
        emit('calculation_update', initial_instance)
        logger.info(f"Sent initial state for calculation {calculation_id} (status: {initial_instance['status']})")
        
        # Store callback reference for cleanup
        session['file_change_callback'] = on_file_change
        
    except Exception as e:
        logger.error(f"Error setting up SocketIO monitoring for {calculation_id}: {e}")
        emit('error', {
            'error': 'Failed to set up calculation monitoring',
            'id': calculation_id
        })


@socketio.on('leave_calculation') 
def on_leave_calculation(data):
    """Leave a calculation room."""
    calculation_id = data.get('calculation_id')
    if not calculation_id:
        return
        
    from flask_socketio import leave_room
    from flask import session
    
    room = f'calculation_{calculation_id}'
    leave_room(room)
    
    # Clean up file watcher connection
    if 'file_change_callback' in session:
        try:
            file_manager = CalculationFileManager()
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.remove_connection(calculation_id, session['file_change_callback'])
        except Exception as e:
            logger.debug(f"Error cleaning up file watcher for {calculation_id}: {e}")
    
    logger.info(f"Client left calculation {calculation_id}")


@socketio.on('join_global_updates')
def on_join_global_updates():
    """Join global updates room to receive all calculation updates."""
    from flask_socketio import join_room
    
    join_room('global_updates')
    logger.info("Client joined global_updates room for real-time monitoring of all calculations")


@socketio.on('leave_global_updates')
def on_leave_global_updates():
    """Leave global updates room."""
    from flask_socketio import leave_room
    
    leave_room('global_updates')
    logger.info("Client left global_updates room")


@socketio.on('disconnect')
def on_disconnect():
    """Handle client disconnection.""" 
    from flask import session
    
    calculation_id = session.get('calculation_id')
    if calculation_id and 'file_change_callback' in session:
        try:
            file_manager = CalculationFileManager()
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.remove_connection(calculation_id, session['file_change_callback'])
            logger.info(f"Cleaned up file watcher for disconnected client (calculation {calculation_id})")
        except Exception as e:
            logger.debug(f"Error cleaning up file watcher on disconnect for {calculation_id}: {e}")


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
    import re
    
    # Check if this is a WebSocket close frame being misinterpreted as HTTP
    error_description = str(error).lower()
    
    # Expanded list of WebSocket protocol indicators
    websocket_indicators = [
        # HTTP method related
        'invalid http method', 'expected get method', 'invalid method', 
        'unexpected http method', 'unsupported method',
        
        # WebSocket specific
        'websocket', 'connection upgrade', 'upgrade required', 
        'websocket handshake', 'sec-websocket',
        
        # Protocol frames and binary data
        'x88x82', '\\x88\\x82', 'x88', '\\x88',  # Close frame
        'x89', '\\x89',  # Ping frame
        'x8a', '\\x8a',  # Pong frame
        'x81', '\\x81',  # Text frame
        'x82', '\\x82',  # Binary frame
        
        # Parser and protocol errors
        'bad request line', 'malformed request', 'protocol error',
        'invalid request line', 'request line too long',
        'invalid header', 'header too long', 'bad http version',
        
        # Connection related
        'connection reset', 'connection closed', 'connection aborted',
        'broken pipe', 'client disconnected',
        
        # Encoding/parsing issues
        'invalid utf-8', 'decode error', 'unicode error',
        'invalid character', 'unexpected character'
    ]
    
    # Regular expression patterns for more sophisticated matching
    websocket_patterns = [
        r'\\x[0-9a-f]{2}',  # Hexadecimal escape sequences (binary data)
        r'x[0-9a-f]{2}',    # Hex bytes without backslash
        r'invalid.*method.*[^a-zA-Z]',  # Invalid method patterns
        r'malformed.*request.*(?:line|header)',  # Malformed request patterns (specific)
        r'connection.*(?:reset|closed|aborted)',  # Connection issues
        r'websocket.*(?:error|frame|close)',  # WebSocket-specific errors
        r'bad.*(?:request\s+line|header)',  # Bad request line/header (specific)
        r'invalid.*(?:request\s+line|header)',  # Invalid request line/header
        r'protocol.*error',  # Protocol errors
        r'unexpected.*(?:character|data|frame)',  # Unexpected data patterns
    ]
    
    # Check simple string indicators first
    is_websocket_related = any(indicator in error_description for indicator in websocket_indicators)
    
    # If not found, check regex patterns
    if not is_websocket_related:
        for pattern in websocket_patterns:
            if re.search(pattern, error_description, re.IGNORECASE):
                is_websocket_related = True
                break
    
    # Additional check: if error description contains mostly non-printable characters
    # This often indicates binary WebSocket frames
    if not is_websocket_related:
        non_printable_count = sum(1 for c in str(error) if ord(c) < 32 and c not in '\n\r\t')
        if non_printable_count > 3:  # Threshold for binary data detection
            is_websocket_related = True
    
    if is_websocket_related:
        # This is likely a WebSocket close frame or protocol error, log as debug
        logger.debug(f"WebSocket protocol frame misinterpreted as HTTP: {error}")
        return '', 400  # Return empty response for WebSocket frames
    
    # For genuine bad requests, return proper error response
    logger.warning(f"Genuine bad request: {error}")
    return jsonify({'success': False, 'error': 'Bad request.'}), 400

@app.errorhandler(404)
def not_found(error):
    return jsonify({'success': False, 'error': 'Endpoint not found.'}), 404

@app.errorhandler(405)
def method_not_allowed(error):
    return jsonify({'success': False, 'error': 'Method not allowed for this endpoint.'}), 405

def create_app():
    """Application factory for Gunicorn compatibility."""
    return app

def create_socketio():
    """SocketIO factory for Gunicorn compatibility."""
    return socketio

def find_available_port(host, start_port, end_port):
    """Find an available port within the specified range."""
    for port in range(start_port, end_port + 1):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind((host, port))
                return port
        except OSError:
            continue
    raise RuntimeError(f"No available port found in range {start_port}-{end_port}")

def start_development_server():
    """Start the server in development mode using Flask-SocketIO."""
    server_config = SERVER_CONFIG.get('server', {})
    dev_config = SERVER_CONFIG.get('development', {})
    socketio_config = SERVER_CONFIG.get('socketio', {})
    
    host = server_config.get('host', '127.0.0.1')
    port_config = server_config.get('port', {})
    debug = dev_config.get('debug', False)
    
    # Command line port argument takes precedence
    port_arg = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    
    if port_arg > 0:
        actual_port = port_arg
    elif port_config.get('auto_detect', True):
        port_range = port_config.get('range', {'start': 5000, 'end': 5100})
        start_port = port_config.get('default', 5000)
        actual_port = find_available_port(host, start_port, port_range['end'])
    else:
        actual_port = port_config.get('default', 5000)
    
    # Notify Electron process of the port
    print(f"FLASK_SERVER_PORT:{actual_port}", file=sys.stdout, flush=True)
    
    logger.info(f"Starting API server with Flask-SocketIO on http://{host}:{actual_port}")
    logger.info(f"Configuration: Debug={debug}, Async_mode={socketio_config.get('async_mode', 'threading')}")
    
    # Start Flask-SocketIO server with configuration-based settings
    try:
        socketio.run(
            app, 
            host=host, 
            port=actual_port, 
            debug=debug, 
            use_reloader=False, 
            allow_unsafe_werkzeug=socketio_config.get('allow_unsafe_werkzeug', True)
        )
    except KeyboardInterrupt:
        logger.info("Server stopped by user.")
    finally:
        logger.info("Shutting down server...")
        cleanup_resources()

if __name__ == '__main__':
    start_development_server()