import json
import time
from pathlib import Path
from typing import Callable

import httpx

from quantum_calc import get_current_settings
from services.calculation_update_stream import get_calculation_update_stream_hub
from websocket.event_loop_bridge import schedule_coroutine


AUTH_HEADERS = {"X-Auth-Token": "socket-token"}


def _wait_until(predicate: Callable[[], bool], timeout: float = 1.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.01)
    assert predicate()


def _extract_sse_payload(line: str) -> dict:
    assert line.startswith("data: ")
    return json.loads(line.removeprefix("data: "))


def test_global_stream_returns_event_stream_content_type(asgi_server) -> None:
    with httpx.Client(timeout=2.0) as client:
        with client.stream(
            "GET",
            f"{asgi_server}/api/quantum/calculations/updates/stream",
            headers=AUTH_HEADERS,
        ) as response:
            assert response.status_code == 200
            assert response.headers["content-type"].startswith("text/event-stream")


def test_calculation_stream_emits_initial_payload(asgi_server) -> None:
    calc_id = "calc-sse-initial"
    calc_dir = Path(get_current_settings().calculations_directory) / calc_id
    calc_dir.mkdir(parents=True)
    (calc_dir / "parameters.json").write_text(
        json.dumps(
            {
                "name": "SSE Calculation",
                "calculation_method": "HF",
                "created_at": "2026-05-26T00:00:00",
            }
        ),
        encoding="utf-8",
    )
    (calc_dir / "status.json").write_text(
        json.dumps({"status": "waiting", "waiting_reason": "All slots are busy"}),
        encoding="utf-8",
    )

    with httpx.Client(timeout=2.0) as client:
        with client.stream(
            "GET",
            f"{asgi_server}/api/quantum/calculations/{calc_id}/updates/stream",
            headers=AUTH_HEADERS,
        ) as response:
            assert response.status_code == 200
            first_line = next(response.iter_lines())

    event = _extract_sse_payload(first_line)
    assert event["type"] == "calculation_update"
    assert event["payload"]["calculation"]["id"] == calc_id
    assert event["payload"]["calculation"]["status"] == "waiting"
    assert event["payload"]["calculation"]["waitingReason"] == "All slots are busy"


def test_calculation_stream_missing_calculation_emits_error(asgi_server) -> None:
    with httpx.Client(timeout=2.0) as client:
        with client.stream(
            "GET",
            f"{asgi_server}/api/quantum/calculations/missing-calc/updates/stream",
            headers=AUTH_HEADERS,
        ) as response:
            assert response.status_code == 200
            lines = list(response.iter_lines())

    data_lines = [line for line in lines if line]
    assert len(data_lines) == 1
    event = _extract_sse_payload(data_lines[0])
    assert event["type"] == "error"
    assert event["payload"]["calculation_id"] == "missing-calc"
    assert "not found" in event["payload"]["message"]


def test_global_stream_receives_published_update(asgi_server) -> None:
    hub = get_calculation_update_stream_hub()

    with httpx.Client(timeout=2.0) as client:
        with client.stream(
            "GET",
            f"{asgi_server}/api/quantum/calculations/updates/stream",
            headers=AUTH_HEADERS,
        ) as response:
            assert response.status_code == 200
            future = schedule_coroutine(
                hub.publish_calculation_update(
                    {"id": "calc-global", "status": "completed"}
                )
            )
            assert future is not None
            future.result(timeout=1.0)
            first_line = next(response.iter_lines())

    event = _extract_sse_payload(first_line)
    assert event["type"] == "calculation_update"
    assert event["payload"]["calculation"]["id"] == "calc-global"


def test_global_stream_unsubscribes_on_disconnect(asgi_server) -> None:
    hub = get_calculation_update_stream_hub()

    with httpx.Client(timeout=2.0) as client:
        with client.stream(
            "GET",
            f"{asgi_server}/api/quantum/calculations/updates/stream",
            headers=AUTH_HEADERS,
        ) as response:
            assert response.status_code == 200
            _wait_until(lambda: len(hub._global_subscribers) == 1)

        _wait_until(lambda: len(hub._global_subscribers) == 0)
