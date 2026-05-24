"""
Quantum chemistry calculation API endpoints.
Handles job submission, monitoring, results retrieval, and orbital/spectrum analysis.
"""

import json
import logging
from datetime import datetime
from typing import Annotated, Any

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import JSONResponse

from generated_models import CalculationUpdateRequest, QuantumCalculationRequest
from quantum_calc.method_defaults import validate_parameters_for_method
from services import get_quantum_service

logger = logging.getLogger(__name__)

router = APIRouter(prefix='/api/quantum')


@router.get('/supported-parameters')
def get_supported_parameters() -> dict[str, Any]:
    """Get supported quantum chemistry parameters."""
    parameters = get_quantum_service().get_supported_parameters()
    return {'success': True, 'data': parameters}


@router.post('/calculate')
async def quantum_calculate(request: Request) -> JSONResponse:
    """Start a quantum chemistry calculation after method-aware validation."""
    raw_body = await request.body()
    if not raw_body:
        raise HTTPException(status_code=400, detail='Request body is required')

    raw_data = json.loads(raw_body)
    if not raw_data:
        raise HTTPException(status_code=400, detail='Request body is required')
    if not isinstance(raw_data, dict):
        raise HTTPException(status_code=400, detail='Request body must be a JSON object')

    calculation_method = raw_data.get('calculation_method')
    if not calculation_method:
        raise HTTPException(status_code=400, detail='calculation_method is required')

    is_valid, applicability_error = validate_parameters_for_method(
        calculation_method,
        raw_data,
    )
    if not is_valid:
        logger.warning("Parameter applicability check failed: %s", applicability_error)
        return JSONResponse(
            {'success': False, 'error': applicability_error},
            status_code=400,
        )

    body = QuantumCalculationRequest.model_validate(raw_data)
    validated_model = body.root if hasattr(body, 'root') else body
    parameters = validated_model.model_dump(exclude_none=False, mode='python')
    for key, value in list(parameters.items()):
        parameters[key] = value.value if hasattr(value, 'value') else value
    parameters['created_at'] = datetime.now().isoformat()

    result = get_quantum_service().start_calculation(parameters)
    return JSONResponse(
        {'success': True, 'data': {'calculation': result}},
        status_code=202,
    )


@router.get('/calculations')
def list_calculations(
    name_query: str | None = None,
    status: str | None = None,
    calculation_method: str | None = None,
    basis_function: str | None = None,
    date_from: str | None = None,
    date_to: str | None = None,
) -> dict[str, Any]:
    """List calculation directories with optional filtering."""
    result = get_quantum_service().list_calculations(
        name_query=name_query,
        status=status,
        calculation_method=calculation_method,
        basis_function=basis_function,
        date_from=date_from,
        date_to=date_to,
    )
    return {'success': True, 'data': result}


@router.get('/status')
def get_calculation_status() -> dict[str, Any]:
    """Get status information about the calculation system."""
    result = get_quantum_service().get_calculation_status()
    return {'success': True, 'data': result}


@router.get('/calculations/{calculation_id}')
def get_calculation_details(calculation_id: str) -> dict[str, Any]:
    """Get detailed information about a specific calculation."""
    result = get_quantum_service().get_calculation_details(calculation_id)
    return {'success': True, 'data': result}


@router.put('/calculations/{calculation_id}')
def update_calculation(
    calculation_id: str,
    body: CalculationUpdateRequest,
) -> dict[str, Any]:
    """Update calculation metadata."""
    result = get_quantum_service().update_calculation(calculation_id, body.name)
    return {'success': True, 'data': result}


@router.post('/calculations/{calculation_id}/pause')
def pause_calculation(calculation_id: str) -> JSONResponse:
    """Pause a running calculation."""
    result = get_quantum_service().pause_calculation(calculation_id)
    return JSONResponse({'success': True, 'data': result}, status_code=202)


@router.post('/calculations/{calculation_id}/resume')
def resume_calculation(calculation_id: str) -> JSONResponse:
    """Resume a paused calculation."""
    result = get_quantum_service().resume_calculation(calculation_id)
    return JSONResponse({'success': True, 'data': result}, status_code=202)


@router.delete('/calculations/{calculation_id}')
def delete_calculation(calculation_id: str) -> dict[str, Any]:
    """Delete a calculation and its files."""
    result = get_quantum_service().delete_calculation(calculation_id)
    return {'success': True, 'data': result}


@router.get('/calculations/{calculation_id}/orbitals')
def get_orbitals(calculation_id: str) -> dict[str, Any]:
    """Get molecular orbital information for a calculation."""
    orbital_summary = get_quantum_service().get_molecular_orbitals(calculation_id)
    return {'success': True, 'data': orbital_summary}


@router.get('/calculations/{calculation_id}/orbitals/{orbital_index}/cube')
def get_orbital_cube(
    calculation_id: str,
    orbital_index: int,
    gridSize: Annotated[int, Query(alias='gridSize')] = 80,
    isovaluePos: Annotated[float | None, Query(alias='isovaluePos')] = None,
    isovalueNeg: Annotated[float | None, Query(alias='isovalueNeg')] = None,
) -> dict[str, Any]:
    """Generate and return CUBE file for a molecular orbital."""
    cube_data = get_quantum_service().generate_orbital_cube(
        calculation_id,
        orbital_index,
        grid_size=gridSize,
        isovalue_pos=isovaluePos,
        isovalue_neg=isovalueNeg,
    )
    return {'success': True, 'data': cube_data}


@router.get('/calculations/{calculation_id}/orbitals/cube-files')
def list_cube_files(calculation_id: str) -> dict[str, Any]:
    """List all CUBE files for a calculation."""
    result = get_quantum_service().list_cube_files(calculation_id)
    return {'success': True, 'data': result}


@router.delete('/calculations/{calculation_id}/orbitals/cube-files')
def delete_cube_files(
    calculation_id: str,
    orbital_index: int | None = None,
) -> dict[str, Any]:
    """Delete CUBE files for a calculation."""
    result = get_quantum_service().delete_cube_files(calculation_id, orbital_index)
    return {'success': True, 'data': result}


@router.get('/calculations/{calculation_id}/ir-spectrum')
def get_ir_spectrum(
    calculation_id: str,
    broadening_fwhm: float = 100.0,
    x_min: float = 400.0,
    x_max: float = 4000.0,
    show_peaks: bool = True,
) -> dict[str, Any]:
    """Generate and return IR spectrum for a calculation."""
    result = get_quantum_service().generate_ir_spectrum(
        calculation_id,
        broadening_fwhm=broadening_fwhm,
        x_min=x_min,
        x_max=x_max,
        show_peaks=show_peaks,
    )
    return {'success': True, 'data': result}
