import asyncio

import pytest

from services.calculation_update_stream import (
    CalculationUpdateStreamHub,
    STREAM_QUEUE_SIZE,
    format_sse_event,
)


@pytest.fixture
def anyio_backend() -> str:
    return "asyncio"


def test_format_sse_event_serializes_json_payload() -> None:
    event = format_sse_event(
        "calculation_update",
        {"calculation": {"id": "calc-1", "status": "running"}},
    )

    assert event.startswith("data: ")
    assert '"type":"calculation_update"' in event
    assert '"id":"calc-1"' in event
    assert event.endswith("\n\n")


@pytest.mark.anyio
async def test_global_subscriber_receives_published_update() -> None:
    hub = CalculationUpdateStreamHub()
    subscriber = hub.subscribe_global()

    await hub.publish_calculation_update(
        {"id": "calc-1", "status": "completed"}
    )

    event = await asyncio.wait_for(subscriber.queue.get(), timeout=1)

    assert event["type"] == "calculation_update"
    assert event["payload"]["calculation"]["id"] == "calc-1"
    hub.unsubscribe(subscriber)


@pytest.mark.anyio
async def test_calculation_subscriber_receives_only_matching_updates() -> None:
    hub = CalculationUpdateStreamHub()
    subscriber = hub.subscribe_calculation("calc-1")

    await hub.publish_calculation_update({"id": "calc-2", "status": "running"})
    await hub.publish_calculation_update({"id": "calc-1", "status": "running"})

    event = await asyncio.wait_for(subscriber.queue.get(), timeout=1)

    assert event["payload"]["calculation"]["id"] == "calc-1"
    assert subscriber.queue.empty()
    hub.unsubscribe(subscriber)


@pytest.mark.anyio
async def test_full_queue_drops_oldest_update_and_keeps_latest() -> None:
    hub = CalculationUpdateStreamHub()
    subscriber = hub.subscribe_global()

    for index in range(STREAM_QUEUE_SIZE + 1):
        await hub.publish_calculation_update(
            {"id": f"calc-{index}", "status": "running"}
        )

    assert subscriber.queue.qsize() == STREAM_QUEUE_SIZE

    retained_calculation_ids = [
        subscriber.queue.get_nowait()["payload"]["calculation"]["id"]
        for _ in range(STREAM_QUEUE_SIZE)
    ]

    assert retained_calculation_ids[0] == "calc-1"
    assert retained_calculation_ids[-1] == f"calc-{STREAM_QUEUE_SIZE}"
    hub.unsubscribe(subscriber)


@pytest.mark.anyio
async def test_unsubscribe_removes_subscriber() -> None:
    hub = CalculationUpdateStreamHub()
    subscriber = hub.subscribe_global()

    hub.unsubscribe(subscriber)
    await hub.publish_calculation_update(
        {"id": "calc-1", "status": "completed"}
    )

    assert subscriber.queue.empty()
