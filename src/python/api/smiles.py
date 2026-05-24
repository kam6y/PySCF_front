"""SMILES conversion API endpoints."""

from fastapi import APIRouter

from generated_models import SMILESConvertRequest
from services import get_smiles_service

router = APIRouter(prefix='/api/smiles')


@router.post('/convert')
def convert_smiles(body: SMILESConvertRequest) -> dict:
    result = get_smiles_service().convert_smiles(body.smiles)
    return {'success': True, 'data': result}
