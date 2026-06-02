"""
System and debug API endpoints.
Handles system resource monitoring and diagnostic information.
"""

import ipaddress
import json
import logging
import os
from datetime import datetime
from typing import Any

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from generated_models import (
    AllocatedResources,
    ResourceConstraints,
    SystemResourceInfo,
    SystemResourceSummary,
)
from services import get_system_service

logger = logging.getLogger(__name__)

router = APIRouter()


def _debug_endpoints_allowed() -> bool:
    """Allow debug diagnostics only in non-production environments.

    The ``PYSCF_ENV`` env var is checked at request time.  When set to
    ``"production"`` the debug endpoints return 403.  In all other cases
    (absent, ``"development"``, etc.) the endpoints are available —
    they still require the existing auth-token middleware.
    """
    return os.getenv("PYSCF_ENV", "").lower() != "production"


def _get_client_host(request: Request) -> str | None:
    return request.client.host if request.client else None


def _is_loopback_address(address: str | None) -> bool:
    if not address:
        return False
    try:
        ip = ipaddress.ip_address(address)
    except ValueError:
        return False
    if ip.is_loopback:
        return True
    if isinstance(ip, ipaddress.IPv6Address) and ip.ipv4_mapped:
        return ip.ipv4_mapped.is_loopback
    return False


@router.get("/api/system/resource-status")
def get_system_resource_status() -> dict[str, Any]:
    """Get current system resource status including constraints and allocation."""
    resource_summary = get_system_service().get_resource_status()
    system_info = SystemResourceInfo(
        total_cpu_cores=resource_summary["system_info"]["total_cpu_cores"],
        total_memory_mb=resource_summary["system_info"]["total_memory_mb"],
        available_memory_mb=resource_summary["system_info"]["available_memory_mb"],
        cpu_usage_percent=resource_summary["system_info"]["cpu_usage_percent"],
        memory_usage_percent=resource_summary["system_info"]["memory_usage_percent"],
        timestamp=datetime.fromisoformat(
            resource_summary["system_info"]["timestamp"].replace("Z", "+00:00")
        ),
    )
    constraints = ResourceConstraints(**resource_summary["resource_constraints"])
    allocated = AllocatedResources(**resource_summary["allocated_resources"])
    summary = SystemResourceSummary(
        system_info=system_info,
        resource_constraints=constraints,
        allocated_resources=allocated,
    )
    return {"success": True, "data": summary.model_dump(mode="json")}


@router.get("/api/system/gpu4pyscf-status")
def get_gpu4pyscf_status() -> dict[str, Any]:
    """Get CUDA detection and GPU4PySCF installation status."""
    status = get_system_service().get_gpu4pyscf_status()
    return {"success": True, "data": status}


@router.post("/api/system/gpu4pyscf-install")
async def install_gpu4pyscf(request: Request) -> Any:
    """Install GPU4PySCF and recommended cuTENSOR from local requests only."""
    client_host = _get_client_host(request)
    if not _is_loopback_address(client_host):
        logger.warning(
            "Blocked non-local GPU4PySCF install request from %s",
            client_host,
        )
        return JSONResponse(
            {
                "success": False,
                "error": "GPU4PySCF installation is only available from the local machine.",
            },
            status_code=403,
        )

    raw_body = await request.body()
    payload = json.loads(raw_body) if raw_body else {}
    if not isinstance(payload, dict):
        payload = {}
    include_cutensor = payload.get("include_cutensor", True) is True
    force_reinstall = payload.get("force_reinstall") is True
    confirm_install = payload.get("confirm_install") is True
    result = get_system_service().install_gpu4pyscf(
        include_cutensor=include_cutensor,
        force_reinstall=force_reinstall,
        confirm_install=confirm_install,
    )
    return {"success": True, "data": result}


_DEBUG_BLOCKED_STATUS_CODE = 403
_DEBUG_BLOCKED_CONTENT = {
    "success": False,
    "error": "Debug endpoints are disabled in production.",
}


def _debug_blocked_response() -> JSONResponse:
    """Create a fresh JSONResponse for blocked debug requests.

    Returns a new object per call to avoid shared mutable state —
    Starlette Response objects are mutated in-place by middleware
    (e.g. CORS raw_headers), and sync endpoints run in a threadpool.
    """
    return JSONResponse(
        content=_DEBUG_BLOCKED_CONTENT, status_code=_DEBUG_BLOCKED_STATUS_CODE
    )


@router.get("/api/debug/system-diagnostics", response_model=None)
def get_system_diagnostics() -> dict[str, Any] | JSONResponse:
    """Get comprehensive system diagnostics for troubleshooting."""
    if not _debug_endpoints_allowed():
        return _debug_blocked_response()
    diagnostics = get_system_service().get_system_diagnostics()
    return {"success": True, "data": diagnostics}


@router.get("/api/debug/process-manager-diagnostics", response_model=None)
def get_process_manager_diagnostics() -> dict[str, Any] | JSONResponse:
    """Get detailed process manager diagnostics."""
    if not _debug_endpoints_allowed():
        return _debug_blocked_response()
    diagnostics = get_system_service().get_process_manager_diagnostics()
    return {"success": True, "data": diagnostics}


@router.get("/api/debug/resource-manager-diagnostics", response_model=None)
def get_resource_manager_diagnostics() -> dict[str, Any] | JSONResponse:
    """Get detailed resource manager diagnostics."""
    if not _debug_endpoints_allowed():
        return _debug_blocked_response()
    diagnostics = get_system_service().get_resource_manager_diagnostics()
    return {"success": True, "data": diagnostics}
