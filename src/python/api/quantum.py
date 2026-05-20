"""
Quantum chemistry calculation API endpoints.
Handles all quantum chemistry calculation operations including job submission,
monitoring, results retrieval, and orbital/spectrum analysis.
"""

import logging
from datetime import datetime
from flask import Blueprint, request, jsonify
from flask_pydantic import validate

from services import get_quantum_service
from generated_models import QuantumCalculationRequest, CalculationUpdateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create quantum blueprint
quantum_bp = Blueprint('quantum', __name__)




@quantum_bp.route('/api/quantum/supported-parameters', methods=['GET'])
def get_supported_parameters():
    """Get supported quantum chemistry parameters including basis functions, exchange-correlation functionals, and solvents."""
    quantum_service = get_quantum_service()

    # Call service layer
    parameters = quantum_service.get_supported_parameters()

    return jsonify({
        'success': True,
        'data': parameters
    }), 200


@quantum_bp.route('/api/quantum/calculate', methods=['POST'])
def quantum_calculate():
    """
    Starts a quantum chemistry calculation in the background.
    Immediately returns a calculation ID to track the job.
    """
    quantum_service = get_quantum_service()

    # Get raw JSON data before Pydantic validation
    raw_data = request.get_json()
    if not raw_data:
        return jsonify({
            'success': False,
            'error': 'Request body is required'
        }), 400

    # Extract calculation method for early validation
    calculation_method = raw_data.get('calculation_method')
    if not calculation_method:
        return jsonify({
            'success': False,
            'error': 'calculation_method is required'
        }), 400

    # Validate parameter applicability BEFORE Pydantic validation
    # This ensures we catch inapplicable parameters that Pydantic would ignore
    from quantum_calc.method_defaults import validate_parameters_for_method
    is_valid, applicability_error = validate_parameters_for_method(
        calculation_method,
        raw_data
    )
    if not is_valid:
        logger.warning(f"Parameter applicability check failed: {applicability_error}")
        return jsonify({
            'success': False,
            'error': applicability_error
        }), 400

    # Now validate with Pydantic
    body = QuantumCalculationRequest.model_validate(raw_data)

    # Extract enum values helper function
    def get_enum_value(field_value):
        if hasattr(field_value, 'value'):
            return field_value.value
        return field_value

    # Handle Pydantic RootModel[Union[...]] structure
    # Access .root attribute if present (discriminated union from OpenAPI)
    validated_model = body.root if hasattr(body, 'root') else body

    # Convert Pydantic model to a plain dictionary
    parameters = validated_model.model_dump(exclude_none=False, mode='python')

    # Convert enum values to strings
    for key, value in list(parameters.items()):
        parameters[key] = get_enum_value(value)

    # Add creation timestamp
    parameters['created_at'] = datetime.now().isoformat()

    # Call service layer (also validates parameters for defense-in-depth and AI agent calls)
    result = quantum_service.start_calculation(parameters)

    return jsonify({'success': True, 'data': {'calculation': result}}), 202


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


@quantum_bp.route('/api/quantum/status', methods=['GET'])
def get_calculation_status():
    """Get status information about the calculation system."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.get_calculation_status()

    return jsonify({
        'success': True,
        'data': result
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['GET'])
def get_calculation_details(calculation_id):
    """Get detailed information about a specific calculation."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.get_calculation_details(calculation_id)

    return jsonify({
        'success': True,
        'data': result
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['PUT'])
@validate()
def update_calculation(calculation_id, body: CalculationUpdateRequest):
    """Update calculation metadata (currently only name)."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.update_calculation(calculation_id, body.name)

    return jsonify({
        'success': True,
        'data': result
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/pause', methods=['POST'])
def pause_calculation(calculation_id):
    """Pause a running calculation."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.pause_calculation(calculation_id)

    return jsonify({
        'success': True,
        'data': result
    }), 202


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/resume', methods=['POST'])
def resume_calculation(calculation_id):
    """Resume a paused calculation."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.resume_calculation(calculation_id)

    return jsonify({
        'success': True,
        'data': result
    }), 202


@quantum_bp.route('/api/quantum/calculations/<calculation_id>', methods=['DELETE'])
def delete_calculation(calculation_id):
    """Delete a calculation and its files."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.delete_calculation(calculation_id)

    return jsonify({
        'success': True,
        'data': result
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals', methods=['GET'])
def get_orbitals(calculation_id):
    """Get molecular orbital information for a calculation."""
    quantum_service = get_quantum_service()

    # Call service layer
    orbital_summary = quantum_service.get_molecular_orbitals(calculation_id)

    return jsonify({
        'success': True,
        'data': orbital_summary
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/<int:orbital_index>/cube', methods=['GET'])
def get_orbital_cube(calculation_id, orbital_index):
    """Generate and return CUBE file for specific molecular orbital."""
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


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['GET'])
def list_cube_files(calculation_id):
    """List all CUBE files for a calculation."""
    quantum_service = get_quantum_service()

    # Call service layer
    result = quantum_service.list_cube_files(calculation_id)

    return jsonify({
        'success': True,
        'data': result
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/orbitals/cube-files', methods=['DELETE'])
def delete_cube_files(calculation_id):
    """Delete CUBE files for a calculation."""
    quantum_service = get_quantum_service()

    # Get query parameters
    orbital_index = request.args.get('orbital_index', type=int)

    # Call service layer
    result = quantum_service.delete_cube_files(calculation_id, orbital_index)

    return jsonify({
        'success': True,
        'data': result
    })


@quantum_bp.route('/api/quantum/calculations/<calculation_id>/ir-spectrum', methods=['GET'])
def get_ir_spectrum(calculation_id):
    """Generate and return IR spectrum for a calculation."""
    quantum_service = get_quantum_service()

    # Get query parameters for spectrum customization
    broadening_fwhm = request.args.get('broadening_fwhm', default=100.0, type=float)
    x_min = request.args.get('x_min', default=400.0, type=float)
    x_max = request.args.get('x_max', default=4000.0, type=float)
    show_peaks_raw = request.args.get('show_peaks')
    show_peaks = show_peaks_raw.lower() in ('true', '1') if show_peaks_raw is not None else True

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
