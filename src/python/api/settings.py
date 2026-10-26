"""Settings management API endpoints."""

from typing import Any, Dict

from fastapi import APIRouter

from generated_models import SettingsUpdateRequest
from services import get_settings_service

router = APIRouter(prefix="/api/settings")

# ---------------------------------------------------------------------------
# Key used to indicate whether an API key is stored (without revealing it)
# ---------------------------------------------------------------------------
_API_KEY_FIELD = "gemini_api_key"
_API_KEY_CONFIGURED_FIELD = "gemini_api_key_configured"


def _mask_api_key_for_response(settings: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy of *settings* with the real API key removed.

    Adds a boolean ``gemini_api_key_configured`` flag so the UI can tell
    whether a key has been stored, without ever exposing the secret over
    the HTTP boundary.

    Args:
        settings: Raw settings dict (may contain the plaintext key).

    Returns:
        A new dict safe for HTTP responses.
    """
    masked = {**settings}
    real_key = masked.get(_API_KEY_FIELD)
    masked[_API_KEY_CONFIGURED_FIELD] = bool(real_key)
    masked[_API_KEY_FIELD] = ""
    return masked


@router.get("")
def get_settings() -> dict:
    settings = get_settings_service().get_settings()
    return {"success": True, "data": {"settings": _mask_api_key_for_response(settings)}}


@router.put("")
def update_settings(body: SettingsUpdateRequest) -> dict:
    new_settings = body.root if hasattr(body, "root") else body
    updated_settings = get_settings_service().update_settings(
        new_settings.model_dump(mode="json")
    )
    return {
        "success": True,
        "data": {"settings": _mask_api_key_for_response(updated_settings)},
    }
