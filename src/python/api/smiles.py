"""SMILES conversion API endpoints."""

from fastapi import APIRouter, HTTPException

from generated_models import SMILESConvertRequest
from services import get_smiles_service

router = APIRouter(prefix="/api/smiles")

# Endpoint-level length cap matching openapi.yaml maxLength constraint.
# Acts as defense-in-depth even when the Pydantic model enforces it.
MAX_SMILES_LENGTH = 10_000


@router.post("/convert")
def convert_smiles(body: SMILESConvertRequest) -> dict:
    if len(body.smiles) > MAX_SMILES_LENGTH:
        raise HTTPException(
            status_code=400,
            detail=f"smiles must be at most {MAX_SMILES_LENGTH} characters",
        )
    result = get_smiles_service().convert_smiles(body.smiles)
    return {"success": True, "data": result}
