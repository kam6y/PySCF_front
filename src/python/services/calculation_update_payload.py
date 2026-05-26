"""Build calculation update payloads for SSE clients."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from quantum_calc import CalculationRepository


def build_calculation_instance(
    calculation_id: str,
    calculation_path: str | Path,
    repository: CalculationRepository,
    status_override: str | None = None,
    error_message: str | None = None,
) -> dict[str, Any]:
    """Build the calculation instance shape consumed by the frontend."""
    path = Path(calculation_path)
    calc_dir = str(path)

    parameters = repository.read_calculation_parameters(calc_dir) or {}
    results = repository.read_calculation_results(calc_dir)
    stored_status, waiting_reason = repository.read_calculation_status_details(calc_dir)
    status = status_override or stored_status
    display_name = repository.get_display_name(parameters)

    current_time = datetime.now().isoformat()
    try:
        mtime = datetime.fromtimestamp(path.stat().st_mtime).isoformat()
    except OSError:
        mtime = current_time

    instance = {
        "id": calculation_id,
        "name": display_name,
        "status": status,
        "createdAt": parameters.get("created_at", mtime),
        "updatedAt": current_time if status_override else mtime,
        "parameters": parameters,
        "results": results,
        "workingDirectory": calc_dir,
    }

    if waiting_reason is not None:
        instance["waitingReason"] = waiting_reason

    if error_message:
        instance["error"] = error_message
        instance["errorMessage"] = error_message

    return instance
