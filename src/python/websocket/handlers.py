"""
WebSocket handlers for real-time calculation monitoring.
Handles client connections for calculation status updates and progress monitoring.
"""

import logging
import os
from datetime import datetime
from typing import Any

from quantum_calc import CalculationRepository, get_websocket_watcher
from websocket.event_loop_bridge import schedule_coroutine

logger = logging.getLogger(__name__)

_sid_state: dict[str, dict[str, Any]] = {}


def _state_for(sid: str) -> dict[str, Any]:
    return _sid_state.setdefault(sid, {"callbacks": {}})


def build_calculation_instance(
    calc_id: str,
    calc_path: str,
    file_manager: CalculationRepository,
) -> dict[str, Any]:
    """Build a complete calculation instance from file system data."""
    try:
        parameters = file_manager.read_calculation_parameters(calc_path) or {}
        results = file_manager.read_calculation_results(calc_path)
        status = file_manager.read_calculation_status(calc_path) or "pending"
        display_name = file_manager.get_display_name(calc_id, parameters)

        try:
            creation_date = parameters.get(
                "created_at",
                datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat(),
            )
            updated_date = datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat()
        except OSError:
            current_time = datetime.now().isoformat()
            creation_date = current_time
            updated_date = current_time

        return {
            "id": calc_id,
            "name": display_name,
            "status": status,
            "createdAt": creation_date,
            "updatedAt": updated_date,
            "parameters": parameters,
            "results": results,
            "workingDirectory": calc_path,
        }
    except Exception:
        logger.exception("Error building calculation instance for %s", calc_id)
        raise


def register_websocket_handlers(sio: Any) -> None:
    """Register all WebSocket event handlers with the Socket.IO AsyncServer."""
    if getattr(sio, "_pyscf_handlers_registered", False):
        logger.debug("WebSocket handlers already registered; skipping.")
        return
    setattr(sio, "_pyscf_handlers_registered", True)

    @sio.event
    async def connect(sid: str, environ: dict[str, Any], auth: dict[str, Any] | None):
        expected_token = os.getenv("PYSCF_AUTH_TOKEN")
        if expected_token and (not auth or auth.get("token") != expected_token):
            logger.warning("Unauthorized Socket.IO connection attempt")
            return False
        _state_for(sid)
        return True

    @sio.event
    async def join_global_updates(sid: str) -> None:
        await sio.enter_room(sid, "global_updates")
        logger.info(
            "Client joined global_updates room for real-time monitoring of all calculations"
        )

    @sio.event
    async def leave_global_updates(sid: str) -> None:
        await sio.leave_room(sid, "global_updates")
        logger.info("Client left global_updates room")

    @sio.event
    async def join_calculation(sid: str, data: dict[str, Any] | None) -> None:
        calculation_id = (data or {}).get("calculation_id")
        if not calculation_id:
            await sio.emit("error", {"error": "calculation_id is required"}, to=sid)
            return

        if calculation_id.startswith("new-calculation-"):
            logger.info(
                "Socket.IO connection attempt for temporary calculation ID: %s",
                calculation_id,
            )
        else:
            logger.info(
                "Socket.IO connection established for calculation %s",
                calculation_id,
            )

        from quantum_calc import get_current_settings

        settings = get_current_settings()
        file_manager = CalculationRepository(base_dir=settings.calculations_directory)
        try:
            calc_path = str(file_manager.resolve_calculation_path(calculation_id))
        except ValueError:
            await sio.emit(
                "error",
                {
                    "error": "Invalid calculation ID.",
                    "id": calculation_id,
                    "is_temporary": calculation_id.startswith("new-calculation-"),
                },
                to=sid,
            )
            return

        if not os.path.isdir(calc_path):
            if calculation_id.startswith("new-calculation-"):
                error_message = (
                    f'Temporary calculation ID "{calculation_id}" does not exist on server.'
                )
            else:
                error_message = f'Calculation "{calculation_id}" not found.'
            await sio.emit(
                "error",
                {
                    "error": error_message,
                    "id": calculation_id,
                    "is_temporary": calculation_id.startswith("new-calculation-"),
                },
                to=sid,
            )
            return

        room = f"calculation_{calculation_id}"
        await sio.enter_room(sid, room)
        state = _state_for(sid)

        def on_file_change(file_data: dict[str, Any]) -> None:
            async def emit_update() -> None:
                try:
                    calculation_instance = build_calculation_instance(
                        calculation_id,
                        calc_path,
                        file_manager,
                    )
                    await sio.emit("calculation_update", calculation_instance, room=room)
                    if calculation_instance["status"] in ["completed", "error"]:
                        logger.info(
                            "Calculation %s finished with status '%s'.",
                            calculation_id,
                            calculation_instance["status"],
                        )
                except Exception:
                    logger.exception(
                        "Error in file change callback for %s",
                        calculation_id,
                    )
                    await sio.emit(
                        "error",
                        {"error": "Failed to read calculation data", "id": calculation_id},
                        room=room,
                    )

            schedule_coroutine(emit_update())

        try:
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.add_connection(calculation_id, on_file_change)
            state["callbacks"][calculation_id] = on_file_change
            initial_instance = build_calculation_instance(
                calculation_id,
                calc_path,
                file_manager,
            )
            await sio.emit("calculation_update", initial_instance, to=sid)
        except Exception:
            logger.exception(
                "Error setting up Socket.IO monitoring for %s",
                calculation_id,
            )
            await sio.emit(
                "error",
                {"error": "Failed to set up calculation monitoring", "id": calculation_id},
                to=sid,
            )

    @sio.event
    async def leave_calculation(sid: str, data: dict[str, Any] | None) -> None:
        calculation_id = (data or {}).get("calculation_id")
        if not calculation_id:
            await sio.emit("error", {"error": "calculation_id is required"}, to=sid)
            return

        room = f"calculation_{calculation_id}"
        await sio.leave_room(sid, room)
        state = _state_for(sid)
        callback = state["callbacks"].pop(calculation_id, None)
        if callback is not None:
            from quantum_calc import get_current_settings

            settings = get_current_settings()
            file_manager = CalculationRepository(base_dir=settings.calculations_directory)
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.remove_connection(calculation_id, callback)
        logger.info("Client left calculation %s", calculation_id)

    @sio.event
    async def disconnect(sid: str) -> None:
        state = _sid_state.pop(sid, {"callbacks": {}})
        if not state["callbacks"]:
            return

        from quantum_calc import get_current_settings

        settings = get_current_settings()
        file_manager = CalculationRepository(base_dir=settings.calculations_directory)
        watcher = get_websocket_watcher(file_manager.get_base_directory())
        for calculation_id, callback in state["callbacks"].items():
            watcher.remove_connection(calculation_id, callback)
            logger.info(
                "Cleaned up file watcher for disconnected client (calculation %s)",
                calculation_id,
            )
