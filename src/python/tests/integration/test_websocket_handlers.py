"""
Handler-level tests for Socket.IO calculation monitoring.

Task 9 adds end-to-end ASGI smoke coverage. These tests call the registered
async handlers directly so the migration can validate handler behavior without
Flask-SocketIO's test client.
"""

import asyncio
import json
from pathlib import Path
from typing import Any

import pytest

from generated_models import AppSettings, Timezone
from websocket.event_loop_bridge import bind_event_loop, clear_event_loop
from websocket.handlers import _sid_state, register_websocket_handlers


class FakeAsyncSocketIO:
    def __init__(self) -> None:
        self.handlers: dict[str, Any] = {}
        self.entered_rooms: list[tuple[str, str]] = []
        self.left_rooms: list[tuple[str, str]] = []
        self.emitted: list[dict[str, Any]] = []

    def event(self, handler: Any) -> Any:
        self.handlers[handler.__name__] = handler
        return handler

    async def enter_room(self, sid: str, room: str) -> None:
        self.entered_rooms.append((sid, room))

    async def leave_room(self, sid: str, room: str) -> None:
        self.left_rooms.append((sid, room))

    async def emit(
        self,
        event: str,
        data: Any = None,
        *,
        to: str | None = None,
        room: str | None = None,
    ) -> None:
        self.emitted.append({"name": event, "args": [data], "to": to, "room": room})


class FakeWatcher:
    def __init__(self) -> None:
        self.callbacks: dict[str, list[Any]] = {}
        self.removed: list[tuple[str, Any]] = []

    def add_connection(self, calculation_id: str, callback: Any) -> None:
        self.callbacks.setdefault(calculation_id, []).append(callback)

    def remove_connection(self, calculation_id: str, callback: Any) -> None:
        self.removed.append((calculation_id, callback))
        callbacks = self.callbacks.get(calculation_id, [])
        if callback in callbacks:
            callbacks.remove(callback)


@pytest.fixture
def sio() -> FakeAsyncSocketIO:
    _sid_state.clear()
    fake_sio = FakeAsyncSocketIO()
    register_websocket_handlers(fake_sio)
    yield fake_sio
    _sid_state.clear()


@pytest.fixture
def watcher() -> FakeWatcher:
    return FakeWatcher()


@pytest.fixture
def calculations_dir(app) -> Path:
    return Path(app.state.CALCULATIONS_DIR)


@pytest.fixture(autouse=True)
def patch_current_settings(mocker, calculations_dir) -> AppSettings:
    settings = AppSettings(
        max_parallel_instances=3,
        max_cpu_utilization_percent=95.0,
        max_memory_utilization_percent=95.0,
        system_total_cores=8,
        system_total_memory_mb=16384,
        calculations_directory=str(calculations_dir),
        timezone=Timezone.UTC,
        gemini_api_key=None,
        research_email=None,
    )
    mocker.patch("quantum_calc.get_current_settings", return_value=settings)
    return settings


@pytest.fixture(autouse=True)
def patch_watcher(mocker, watcher) -> FakeWatcher:
    mocker.patch("websocket.handlers.get_websocket_watcher", return_value=watcher)
    return watcher


def write_calculation(
    calculations_dir: Path,
    calc_id: str,
    *,
    status: str = "completed",
    parameters: dict[str, Any] | None = None,
    results: dict[str, Any] | None = None,
) -> Path:
    calc_dir = calculations_dir / calc_id
    calc_dir.mkdir(parents=True, exist_ok=True)
    params = {
        "name": "Test Calculation",
        "calculation_method": "HF",
        "basis_function": "sto-3g",
        "created_at": "2024-01-01T00:00:00",
        **(parameters or {}),
    }
    (calc_dir / "parameters.json").write_text(json.dumps(params))
    (calc_dir / "status.json").write_text(json.dumps({"status": status}))
    if results is not None:
        (calc_dir / "results.json").write_text(json.dumps(results))
    return calc_dir


def run(handler: Any, *args: Any) -> Any:
    return asyncio.run(handler(*args))


class TestWebSocketAuthentication:
    def test_connect_without_token_rejected_when_auth_token_configured(
        self,
        sio,
        monkeypatch,
    ):
        monkeypatch.setenv("PYSCF_AUTH_TOKEN", "socket-secret")

        connected = run(sio.handlers["connect"], "sid-1", {}, None)

        assert connected is False

    def test_connect_with_token_allowed_when_auth_token_configured(
        self,
        sio,
        monkeypatch,
    ):
        monkeypatch.setenv("PYSCF_AUTH_TOKEN", "socket-secret")

        connected = run(
            sio.handlers["connect"],
            "sid-1",
            {},
            {"token": "socket-secret"},
        )

        assert connected is True
        assert "sid-1" in _sid_state


class TestJoinCalculationWebSocket:
    def test_join_calculation_success(self, sio, calculations_dir, watcher):
        calc_id = "test-calc-123"
        write_calculation(calculations_dir, calc_id)

        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})

        assert ("sid-1", f"calculation_{calc_id}") in sio.entered_rooms
        assert calc_id in watcher.callbacks
        update_events = [
            msg for msg in sio.emitted if msg["name"] == "calculation_update"
        ]
        assert len(update_events) == 1
        initial_state = update_events[0]["args"][0]
        assert initial_state["id"] == calc_id
        assert initial_state["name"] == "Test Calculation"
        assert initial_state["status"] == "completed"

    def test_join_calculation_not_found(self, sio):
        calc_id = "nonexistent-calc-999"

        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})

        error_events = [msg for msg in sio.emitted if msg["name"] == "error"]
        assert len(error_events) == 1
        error_data = error_events[0]["args"][0]
        assert calc_id in error_data["error"] or "not found" in error_data["error"].lower()
        assert error_events[0]["to"] == "sid-1"

    def test_join_calculation_missing_id(self, sio):
        run(sio.handlers["join_calculation"], "sid-1", {})

        error_events = [msg for msg in sio.emitted if msg["name"] == "error"]
        assert len(error_events) == 1
        assert error_events[0]["args"][0]["error"] == "calculation_id is required"
        assert error_events[0]["to"] == "sid-1"

    def test_join_calculation_temporary_id(self, sio):
        calc_id = "new-calculation-temp-123"

        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})

        error_events = [msg for msg in sio.emitted if msg["name"] == "error"]
        assert len(error_events) == 1
        assert error_events[0]["args"][0]["is_temporary"] is True


class TestLeaveCalculationWebSocket:
    def test_leave_calculation_success(self, sio, calculations_dir, watcher):
        calc_id = "test-calc-leave"
        write_calculation(calculations_dir, calc_id, status="running")
        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})
        sio.emitted.clear()

        run(sio.handlers["leave_calculation"], "sid-1", {"calculation_id": calc_id})

        assert ("sid-1", f"calculation_{calc_id}") in sio.left_rooms
        assert len(watcher.removed) == 1
        assert watcher.removed[0][0] == calc_id
        assert [msg for msg in sio.emitted if msg["name"] == "error"] == []

    def test_leave_calculation_missing_id_returns_error(self, sio):
        run(sio.handlers["leave_calculation"], "sid-1", {})

        error_events = [msg for msg in sio.emitted if msg["name"] == "error"]
        assert len(error_events) == 1
        assert error_events[0]["args"][0]["error"] == "calculation_id is required"


class TestGlobalUpdatesWebSocket:
    def test_join_global_updates(self, sio):
        run(sio.handlers["join_global_updates"], "sid-1")

        assert ("sid-1", "global_updates") in sio.entered_rooms
        assert [msg for msg in sio.emitted if msg["name"] == "error"] == []

    def test_leave_global_updates(self, sio):
        run(sio.handlers["leave_global_updates"], "sid-1")

        assert ("sid-1", "global_updates") in sio.left_rooms
        assert [msg for msg in sio.emitted if msg["name"] == "error"] == []


class TestWebSocketDisconnection:
    def test_disconnect_cleans_up(self, sio, calculations_dir, watcher):
        calc_id = "test-calc-disconnect"
        write_calculation(calculations_dir, calc_id, status="running")
        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})

        run(sio.handlers["disconnect"], "sid-1")

        assert len(watcher.removed) == 1
        assert watcher.removed[0][0] == calc_id
        assert "sid-1" not in _sid_state


class TestCalculationUpdates:
    def test_receive_update_on_file_change(self, sio, calculations_dir, watcher):
        async def main() -> None:
            bind_event_loop(asyncio.get_running_loop())
            try:
                calc_id = "test-calc-update"
                write_calculation(calculations_dir, calc_id, status="running")
                await sio.handlers["join_calculation"](
                    "sid-1",
                    {"calculation_id": calc_id},
                )
                sio.emitted.clear()

                (calculations_dir / calc_id / "status.json").write_text(
                    json.dumps({"status": "completed"})
                )
                watcher.callbacks[calc_id][0]({"path": "status.json"})
                await asyncio.sleep(0.05)

                update_events = [
                    msg for msg in sio.emitted if msg["name"] == "calculation_update"
                ]
                assert len(update_events) == 1
                calc_data = update_events[0]["args"][0]
                assert calc_data["id"] == calc_id
                assert calc_data["status"] == "completed"
                assert Path(calc_data["workingDirectory"]) == (
                    calculations_dir / calc_id
                ).resolve()
            finally:
                clear_event_loop()

        asyncio.run(main())

    def test_calculation_update_contains_all_fields(self, sio, calculations_dir):
        calc_id = "test-calc-complete"
        write_calculation(
            calculations_dir,
            calc_id,
            status="completed",
            parameters={
                "name": "Complete Calculation",
                "exchange_correlation": "b3lyp",
                "charges": 0,
                "spin": 0,
                "xyz": "H 0 0 0\nH 0 0 0.74",
            },
            results={"energy": -1.06, "dipole": [0, 0, 0]},
        )

        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})

        update_events = [
            msg for msg in sio.emitted if msg["name"] == "calculation_update"
        ]
        assert len(update_events) == 1
        calc_data = update_events[0]["args"][0]
        assert set(calc_data) >= {
            "id",
            "name",
            "status",
            "createdAt",
            "updatedAt",
            "parameters",
            "results",
            "workingDirectory",
        }
        assert calc_data["id"] == calc_id
        assert calc_data["name"] == "Complete Calculation"
        assert calc_data["status"] == "completed"
        assert calc_data["results"]["energy"] == -1.06

    def test_calculation_update_with_error(self, sio, calculations_dir):
        calc_id = "test-calc-error"
        write_calculation(
            calculations_dir,
            calc_id,
            status="error",
            parameters={"name": "Error Calculation"},
        )

        run(sio.handlers["join_calculation"], "sid-1", {"calculation_id": calc_id})

        update_events = [
            msg for msg in sio.emitted if msg["name"] == "calculation_update"
        ]
        assert len(update_events) == 1
        assert update_events[0]["args"][0]["status"] == "error"
