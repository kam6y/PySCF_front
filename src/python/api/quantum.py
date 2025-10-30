"""
Quantum chemistry calculation API endpoints.
Handles all quantum chemistry calculation operations including job submission,
monitoring, results retrieval, and orbital/spectrum analysis.
"""

import logging
from datetime import datetime
from flask import Blueprint, request, jsonify
from flask_pydantic import validate

from services import get_quantum_service, ServiceError
from generated_models import QuantumCalculationRequest, CalculationUpdateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create quantum blueprint
quantum_bp = Blueprint('quantum', __name__)




@quantum_bp.route('/api/quantum/supported-parameters', methods=['GET'])
def get_supported_parameters():
    """Get supported quantum chemistry parameters including basis functions, exchange-correlation functionals, and solvents."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        parameters = quantum_service.get_supported_parameters()
        
        return jsonify({
            'success': True,
            'data': parameters
        }), 200
        
    except ServiceError as e:
        logger.error(f"Service error getting supported parameters: {e}")
        return jsonify({
            'success': False,
            'error': e.message
        }), e.status_code
    except Exception as e:
        logger.error(f"Error getting supported parameters: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'An internal server error occurred.'
        }), 500


@quantum_bp.route('/api/quantum/calculate', methods=['POST'])
@validate()
def quantum_calculate(body: QuantumCalculationRequest):
    """
    Starts a quantum chemistry calculation in the background.
    Immediately returns a calculation ID to track the job.
    """
    try:
        quantum_service = get_quantum_service()
        
        # Extract enum values and build parameters
        def get_enum_value(field_value):
            if hasattr(field_value, 'value'):
                return field_value.value
            return field_value
        
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
            'optimize_geometry': body.optimize_geometry,
            'tddft_nstates': body.tddft_nstates,
            'tddft_method': get_enum_value(body.tddft_method) if body.tddft_method else 'TDDFT',
            'tddft_analyze_nto': body.tddft_analyze_nto,
            'ncas': body.ncas,
            'nelecas': body.nelecas,
            'max_cycle_macro': body.max_cycle_macro,
            'max_cycle_micro': body.max_cycle_micro,
            'natorb': body.natorb,
            'conv_tol': body.conv_tol,
            'conv_tol_grad': body.conv_tol_grad
        }
        
        # Add exchange_correlation only for DFT methods
        if calculation_method != 'HF':
            parameters['exchange_correlation'] = body.exchange_correlation
        else:
            parameters['exchange_correlation'] = None
        
        # Call service layer
        result = quantum_service.start_calculation(parameters)
        
        return jsonify({'success': True, 'data': {'calculation': result}}), 202
    
    except ServiceError as e:
        logger.error(f"Service error starting calculation: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error queuing calculation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations', methods=['GET'])
def list_calculations():
    """
    List calculation directories with optional filtering.

    Query Parameters:
        name_query (str, optional): Partial match search in calculation name (case-insensitive)
        status (str, optional): Filter by status ("completed", "running", "error")
        calculation_method (str, optional): Filter by calculation method ("DFT", "HF", "MP2", etc.)
        basis_function (str, optional): Filter by basis set (case-insensitive)
        date_from (str, optional): Start date for date range filtering (ISO format: YYYY-MM-DD)
        date_to (str, optional): End date for date range filtering (ISO format: YYYY-MM-DD)
    """
    try:
        quantum_service = get_quantum_service()

        # Get query parameters for filtering
        name_query = request.args.get('name_query', type=str)
        status = request.args.get('status', type=str)
        calculation_method = request.args.get('calculation_method', type=str)
        basis_function = request.args.get('basis_function', type=str)
        date_from = request.args.get('date_from', type=str)
        date_to = request.args.get('date_to', type=str)

        # Call service layer with filters
        result = quantum_service.list_calculations(
            name_query=name_query,
            status=status,
            calculation_method=calculation_method,
            basis_function=basis_function,
            date_from=date_from,
            date_to=date_to
        )

        return jsonify({
            'success': True,
            'data': result
        })

    except ServiceError as e:
        logger.error(f"Service error listing calculations: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error listing calculations: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/status', methods=['GET'])
def get_calculation_status():
    """Get status information about the calculation system."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        result = quantum_service.get_calculation_status()
        
        return jsonify({
            'success': True,
            'data': result
        })
        
    except ServiceError as e:
        logger.error(f"Service error getting calculation status: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error getting calculation status: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['GET'])
def get_calculation_details(calculation_id):
    """Get detailed information about a specific calculation."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        result = quantum_service.get_calculation_details(calculation_id)
        
        return jsonify({
            'success': True,
            'data': result
        })
        
    except ServiceError as e:
        logger.error(f"Service error getting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error getting calculation details for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['PUT'])
@validate()
def update_calculation(calculation_id, body: CalculationUpdateRequest):
    """Update calculation metadata (currently only name)."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        result = quantum_service.update_calculation(calculation_id, body.name)
        
        return jsonify({
            'success': True,
            'data': result
        })

    except ServiceError as e:
        logger.error(f"Service error updating calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error updating calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/cancel', methods=['POST'])
def cancel_calculation(calculation_id):
    """Cancel a running calculation."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        result = quantum_service.cancel_calculation(calculation_id)
        
        return jsonify({
            'success': True,
            'data': result
        })
            
    except ServiceError as e:
        logger.error(f"Service error cancelling calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error cancelling calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['DELETE'])
def delete_calculation(calculation_id):
    """Delete a calculation and its files."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        result = quantum_service.delete_calculation(calculation_id)
        
        return jsonify({
            'success': True,
            'data': result
        })
        
    except ServiceError as e:
        logger.error(f"Service error deleting calculation {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error deleting calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals', methods=['GET'])
def get_orbitals(calculation_id):
    """Get molecular orbital information for a calculation."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        orbital_summary = quantum_service.get_molecular_orbitals(calculation_id)
        
        return jsonify({
            'success': True,
            'data': orbital_summary
        })
        
    except ServiceError as e:
        logger.error(f"Service error getting orbitals for {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error getting orbitals for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/<int:orbital_index>/cube', methods=['GET'])
def get_orbital_cube(calculation_id, orbital_index):
    """Generate and return CUBE file for specific molecular orbital."""
    try:
        quantum_service = get_quantum_service()
        
        # Get parameters from query parameters with default values
        grid_size = request.args.get('gridSize', default=80, type=int)
        isovalue_pos = request.args.get('isovaluePos', type=float)
        isovalue_neg = request.args.get('isovalueNeg', type=float)
        
        # Call service layer
        cube_data = quantum_service.generate_orbital_cube(
            calculation_id,
            orbital_index,
            grid_size=grid_size,
            isovalue_pos=isovalue_pos,
            isovalue_neg=isovalue_neg
        )
        
        return jsonify({
            'success': True,
            'data': cube_data
        })
        
    except ServiceError as e:
        logger.error(f"Service error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['GET'])
def list_cube_files(calculation_id):
    """List all CUBE files for a calculation."""
    try:
        quantum_service = get_quantum_service()
        
        # Call service layer
        result = quantum_service.list_cube_files(calculation_id)
        
        return jsonify({
            'success': True,
            'data': result
        })
        
    except ServiceError as e:
        logger.error(f"Service error listing CUBE files for {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error listing CUBE files for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['DELETE'])
def delete_cube_files(calculation_id):
    """Delete CUBE files for a calculation."""
    try:
        quantum_service = get_quantum_service()
        
        # Get query parameters
        orbital_index = request.args.get('orbital_index', type=int)
        
        # Call service layer
        result = quantum_service.delete_cube_files(calculation_id, orbital_index)
        
        return jsonify({
            'success': True,
            'data': result
        })
        
    except ServiceError as e:
        logger.error(f"Service error deleting CUBE files for {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error deleting CUBE files for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/ir-spectrum', methods=['GET'])
def get_ir_spectrum(calculation_id):
    """Generate and return IR spectrum for a calculation."""
    try:
        quantum_service = get_quantum_service()
        
        # Get query parameters for spectrum customization
        broadening_fwhm = request.args.get('broadening_fwhm', default=100.0, type=float)
        x_min = request.args.get('x_min', default=400.0, type=float)
        x_max = request.args.get('x_max', default=4000.0, type=float)
        show_peaks = request.args.get('show_peaks', default=True, type=bool)
        
        # Call service layer
        result = quantum_service.generate_ir_spectrum(
            calculation_id,
            broadening_fwhm=broadening_fwhm,
            x_min=x_min,
            x_max=x_max,
            show_peaks=show_peaks
        )
        
        return jsonify({
            'success': True,
            'data': result
        })
        
    except ServiceError as e:
        logger.error(f"Service error generating IR spectrum for {calculation_id}: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error generating IR spectrum for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500