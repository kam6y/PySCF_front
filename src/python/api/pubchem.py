"""PubChem API endpoints."""

from fastapi import APIRouter

from generated_models import PubChemSearchRequest, XYZValidateRequest
from services import get_pubchem_service

router = APIRouter(prefix='/api/pubchem')


@router.post('/search')
def search_pubchem(body: PubChemSearchRequest) -> dict:
    search_type_value = body.searchType or "name"
    search_type = (
        search_type_value.value
        if hasattr(search_type_value, "value")
        else str(search_type_value)
    )
    result = get_pubchem_service().search_compound(body.query, search_type)
    return {'success': True, 'data': result}


@router.post('/validate')
def validate_xyz_endpoint(body: XYZValidateRequest) -> dict:
    validation_result = get_pubchem_service().validate_xyz(body.xyz)
    return {'success': True, 'data': validation_result}
