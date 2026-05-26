import asyncio
import json

import pytest

from services.calculation_update_stream import get_calculation_update_stream_hub
from services.notification_service import NotificationService
from websocket.event_loop_bridge import bind_event_loop, clear_event_loop


@pytest.fixture
def anyio_backend() -> str:
    return "asyncio"


@pytest.mark.anyio
async def test_send_calculation_update_publishes_to_sse_hub(
    tmp_path,
    mocker,
) -> None:
    calc_id = "calc-notification"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    (calc_dir / "parameters.json").write_text(
        json.dumps({"name": "Notification Calc", "created_at": "2026-05-26T00:00:00"}),
        encoding="utf-8",
    )
    (calc_dir / "status.json").write_text(
        json.dumps({"status": "running"}),
        encoding="utf-8",
    )

    settings = mocker.Mock(calculations_directory=str(tmp_path))
    mocker.patch("quantum_calc.get_current_settings", return_value=settings)

    hub = get_calculation_update_stream_hub()
    subscriber = hub.subscribe_global()
    service = NotificationService()
    bind_event_loop(asyncio.get_running_loop())

    try:
        service.send_calculation_update(calc_id, "completed")
        event = await asyncio.wait_for(subscriber.queue.get(), timeout=1)
        assert event["payload"]["calculation"]["id"] == calc_id
        assert event["payload"]["calculation"]["status"] == "completed"
    finally:
        clear_event_loop()
        hub.unsubscribe(subscriber)
