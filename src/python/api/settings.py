"""Settings management API endpoints."""

from fastapi import APIRouter

from generated_models import SettingsUpdateRequest
from services import get_settings_service

router = APIRouter(prefix='/api/settings')


@router.get('')
def get_settings() -> dict:
    settings = get_settings_service().get_settings()
    return {'success': True, 'data': {'settings': settings}}


@router.put('')
def update_settings(body: SettingsUpdateRequest) -> dict:
    new_settings = body.root if hasattr(body, 'root') else body
    updated_settings = get_settings_service().update_settings(new_settings.model_dump())
    return {'success': True, 'data': {'settings': updated_settings}}
