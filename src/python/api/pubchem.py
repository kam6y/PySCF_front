"""PubChem API endpoints."""

from fastapi import APIRouter, HTTPException

from generated_models import PubChemSearchRequest, XYZValidateRequest
from services import get_pubchem_service

router = APIRouter(prefix="/api/pubchem")

# Endpoint-level length caps matching openapi.yaml maxLength constraints.
# These act as defense-in-depth even when the Pydantic model enforces them.
MAX_QUERY_LENGTH = 500
MAX_XYZ_LENGTH = 1_000_000


@router.post("/search")
def search_pubchem(body: PubChemSearchRequest) -> dict:
    if len(body.query) > MAX_QUERY_LENGTH:
        raise HTTPException(
            status_code=400,
            detail=f"query must be at most {MAX_QUERY_LENGTH} characters",
        )
    search_type_value = body.searchType or "name"
    search_type = (
        search_type_value.value
        if hasattr(search_type_value, "value")
        else str(search_type_value)
    )
    result = get_pubchem_service().search_compound(body.query, search_type)
    return {"success": True, "data": result}


@router.post("/validate")
def validate_xyz_endpoint(body: XYZValidateRequest) -> dict:
    if len(body.xyz) > MAX_XYZ_LENGTH:
        raise HTTPException(
            status_code=400,
            detail=f"xyz must be at most {MAX_XYZ_LENGTH} characters",
        )
    validation_result = get_pubchem_service().validate_xyz(body.xyz)
    return {"success": True, "data": validation_result}
