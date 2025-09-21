"""
Quantum chemistry calculation API endpoints.
Handles all quantum chemistry calculation operations including job submission,
monitoring, results retrieval, and orbital/spectrum analysis.
"""

import logging
import os
import time
import shutil
from datetime import datetime
from flask import Blueprint, request, jsonify
from flask_pydantic import validate

from quantum_calc import (
    DFTCalculator, HFCalculator, MP2Calculator, CCSDCalculator, TDDFTCalculator,
    MolecularOrbitalGenerator, CalculationError, ConvergenceError, InputError,
    GeometryError, ProcessManagerError, get_process_manager, get_all_supported_parameters
)
from quantum_calc.ir_spectrum import create_ir_spectrum_from_calculation_results
from quantum_calc.exceptions import XYZValidationError, FileManagerError, ProcessManagerError
from quantum_calc.file_manager import CalculationFileManager
from generated_models import QuantumCalculationRequest, CalculationUpdateRequest, OrbitalCubeRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create quantum blueprint
quantum_bp = Blueprint('quantum', __name__)


def validate_calculation_parameters(body: QuantumCalculationRequest) -> str:
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


@quantum_bp.route('/api/quantum/supported-parameters', methods=['GET'])
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


@quantum_bp.route('/api/quantum/calculate', methods=['POST'])
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
        
        # Add geometry optimization parameter
        parameters['optimize_geometry'] = body.optimize_geometry
        
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

        except Exception as pm_error:
            logger.error(f"Failed to initialize process manager: {pm_error}")
            # Clean up created directory on process manager failure
            try:
                shutil.rmtree(calc_dir, ignore_errors=True)
            except Exception:
                pass
            return jsonify({'success': False, 'error': f'System initialization error: Unable to initialize calculation system. Please check system resources and try again.'}), 503

        # Submit calculation to process pool with enhanced error handling
        try:
            logger.info(f"About to submit calculation {calculation_id}")
            success, initial_status, waiting_reason = process_manager.submit_calculation(calculation_id, parameters)
            logger.info(f"Submit result: success={success}, status={initial_status}, reason={waiting_reason}")

        except Exception as submit_error:
            logger.error(f"Unexpected error during calculation submission: {submit_error}")
            # Update status to error and save error information
            file_manager.save_calculation_status(calc_dir, 'error')
            error_message = f'Failed to submit calculation: {str(submit_error)}'
            file_manager.save_calculation_results(calc_dir, {'error': error_message})
            
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
            
            return jsonify(error_instance), 202

        # Set initial status based on submission result
        file_manager.save_calculation_status(calc_dir, initial_status, waiting_reason)
        
        logger.info(f"Queued calculation {calculation_id} for molecule '{parameters['name']}'")

        # Return immediately with the new calculation instance
        initial_instance = {
            'id': calculation_id,
            'name': parameters['name'],
            'status': initial_status,
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


@quantum_bp.route('/api/quantum/calculations', methods=['GET'])
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


@quantum_bp.route('/api/quantum/status', methods=['GET'])
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
                    'cpu_count': os.cpu_count()
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['GET'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['PUT'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/cancel', methods=['POST'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['DELETE'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals', methods=['GET'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/<int:orbital_index>/cube', methods=['GET'])
@validate()
def get_orbital_cube(calculation_id, orbital_index, query: OrbitalCubeRequest):
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
        
        # Get parameters from Pydantic model with validation already applied
        grid_size = query.gridSize
        isovalue_pos = query.isovaluePos
        isovalue_neg = query.isovalueNeg
        
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['GET'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['DELETE'])
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/ir-spectrum', methods=['GET'])
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