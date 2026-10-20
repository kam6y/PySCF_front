"""SSE routes for calculation update monitoring."""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from pathlib import Path

from fastapi import APIRouter, Request
from fastapi.responses import StreamingResponse

from quantum_calc import CalculationRepository, get_current_settings
from services.calculation_update_payload import build_calculation_instance
from services.calculation_update_stream import (
    StreamSubscriber,
    format_sse_event,
    get_calculation_update_stream_hub,
)

router = APIRouter(prefix="/api/quantum/calculations")

HEARTBEAT_SECONDS = 15.0


async def _event_stream(
    request: Request,
    subscriber: StreamSubscriber,
) -> AsyncIterator[str]:
    hub = get_calculation_update_stream_hub()
    try:
        while True:
            if await request.is_disconnected():
                break
            try:
                event = await asyncio.wait_for(
                    subscriber.queue.get(),
                    timeout=HEARTBEAT_SECONDS,
                )
                yield format_sse_event(event["type"], event["payload"])
            except asyncio.TimeoutError:
                yield format_sse_event("heartbeat", {})
    finally:
        hub.unsubscribe(subscriber)


async def _single_event_stream(
    event_type: str,
    payload: dict[str, object],
) -> AsyncIterator[str]:
    yield format_sse_event(event_type, payload)


@router.get("/updates/stream")
async def stream_calculation_updates(request: Request) -> StreamingResponse:
    subscriber = get_calculation_update_stream_hub().subscribe_global()
    return StreamingResponse(
        _event_stream(request, subscriber),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache"},
    )


@router.get("/{calculation_id}/updates/stream")
async def stream_calculation_updates_for_calculation(
    calculation_id: str,
    request: Request,
) -> StreamingResponse:
    settings = get_current_settings()
    repository = CalculationRepository(base_dir=settings.calculations_directory)
    try:
        calc_path = str(repository.resolve_calculation_path(calculation_id))
    except ValueError:
        return StreamingResponse(
            _single_event_stream(
                "error",
                {
                    "message": "Invalid calculation ID.",
                    "calculation_id": calculation_id,
                },
            ),
            media_type="text/event-stream",
            headers={"Cache-Control": "no-cache"},
        )

    if not Path(calc_path).is_dir():
        return StreamingResponse(
            _single_event_stream(
                "error",
                {
                    "message": f"Calculation {calculation_id} not found.",
                    "calculation_id": calculation_id,
                },
            ),
            media_type="text/event-stream",
            headers={"Cache-Control": "no-cache"},
        )

    initial_calculation = build_calculation_instance(
        calculation_id,
        calc_path,
        repository,
    )
    subscriber = get_calculation_update_stream_hub().subscribe_calculation(
        calculation_id
    )
    subscriber.queue.put_nowait(
        {
            "type": "calculation_update",
            "payload": {"calculation": initial_calculation},
        }
    )

    return StreamingResponse(
        _event_stream(request, subscriber),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache"},
    )
