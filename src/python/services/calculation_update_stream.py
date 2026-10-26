"""SSE stream hub for calculation updates."""

from __future__ import annotations

import asyncio
import json
import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

STREAM_QUEUE_SIZE = 20

StreamEvent = dict[str, Any]


@dataclass(frozen=True)
class StreamSubscriber:
    queue: asyncio.Queue[StreamEvent]
    calculation_id: str | None = None


def format_sse_event(event_type: str, payload: dict[str, Any]) -> str:
    data = json.dumps(
        {"type": event_type, "payload": payload},
        ensure_ascii=False,
        separators=(",", ":"),
    )
    return f"data: {data}\n\n"


class CalculationUpdateStreamHub:
    """In-process publisher for calculation update SSE subscribers."""

    def __init__(self) -> None:
        self._global_subscribers: set[StreamSubscriber] = set()
        self._calculation_subscribers: dict[str, set[StreamSubscriber]] = {}

    def subscribe_global(self) -> StreamSubscriber:
        subscriber = StreamSubscriber(queue=asyncio.Queue(maxsize=STREAM_QUEUE_SIZE))
        self._global_subscribers.add(subscriber)
        return subscriber

    def subscribe_calculation(self, calculation_id: str) -> StreamSubscriber:
        subscriber = StreamSubscriber(
            queue=asyncio.Queue(maxsize=STREAM_QUEUE_SIZE),
            calculation_id=calculation_id,
        )
        self._calculation_subscribers.setdefault(calculation_id, set()).add(subscriber)
        return subscriber

    def unsubscribe(self, subscriber: StreamSubscriber) -> None:
        if subscriber.calculation_id is None:
            self._global_subscribers.discard(subscriber)
            return

        subscribers = self._calculation_subscribers.get(subscriber.calculation_id)
        if subscribers is None:
            return
        subscribers.discard(subscriber)
        if not subscribers:
            self._calculation_subscribers.pop(subscriber.calculation_id, None)

    async def publish_calculation_update(self, calculation: dict[str, Any]) -> None:
        calculation_id = calculation.get("id")
        if not isinstance(calculation_id, str) or not calculation_id:
            logger.warning("Dropped calculation update without id: %s", calculation)
            return

        event = {
            "type": "calculation_update",
            "payload": {"calculation": calculation},
        }
        subscribers = set(self._global_subscribers)
        subscribers.update(self._calculation_subscribers.get(calculation_id, set()))

        for subscriber in subscribers:
            self._put_latest(subscriber.queue, event)

    def _put_latest(
        self,
        queue: asyncio.Queue[StreamEvent],
        event: StreamEvent,
    ) -> None:
        if queue.full():
            try:
                queue.get_nowait()
            except asyncio.QueueEmpty:
                pass
        queue.put_nowait(event)


_stream_hub = CalculationUpdateStreamHub()


def get_calculation_update_stream_hub() -> CalculationUpdateStreamHub:
    return _stream_hub
