# SSE Calculation Updates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Socket.IO calculation monitoring with FastAPI Server-Sent Events while preserving real-time calculation update behavior.

**Architecture:** Add FastAPI `StreamingResponse` endpoints backed by a small in-process calculation update hub. Refactor `NotificationService` to publish calculation payloads to that hub, and replace the frontend Socket.IO transport with `@microsoft/fetch-event-source` streams. Remove Socket.IO dependencies after equivalent SSE tests pass.

**Tech Stack:** FastAPI, Starlette `StreamingResponse`, Python `asyncio`, React hooks, `@microsoft/fetch-event-source`, React Query, OpenAPI codegen.

---

## Source References

- Design spec: `docs/superpowers/specs/2026-05-26-sse-calculation-updates-design.md`
- FastAPI docs checked through Context7: `StreamingResponse` accepts sync or async iterators and streams chunks with `media_type="text/event-stream"`.
- Existing SSE client pattern: `src/web/api/agent.ts`
- Existing Socket.IO update processing: `src/web/websocket/useCalculationSync.ts`
- Existing notification source: `src/python/services/notification_service.py`

## File Structure

- Modify `src/api-spec/openapi.yaml`: add two SSE endpoints and schemas.
- Regenerate `src/python/generated_models.py` and `src/web/types/generated-api.ts` with `npm run codegen`.
- Create `src/python/services/calculation_update_stream.py`: async subscriber hub and SSE formatting helpers.
- Create `src/python/services/calculation_update_payload.py`: shared calculation-instance payload builder for SSE routes and notification publishing.
- Create `src/python/api/calculation_updates.py`: FastAPI streaming route handlers.
- Modify `src/python/api/__init__.py`: register calculation update routes.
- Modify `src/python/services/notification_service.py`: publish to SSE hub instead of Socket.IO.
- Modify `src/python/app.py`: remove Socket.IO composition and expose plain FastAPI ASGI app.
- Modify `src/python/tests/conftest.py`: replace Socket.IO ASGI test fixture with plain FastAPI Uvicorn fixture.
- Create `src/python/tests/unit/test_services/test_calculation_update_stream.py`: hub behavior tests.
- Create `src/python/tests/integration/test_api_endpoints/test_calculation_updates_api.py`: SSE endpoint tests.
- Delete `src/python/tests/integration/test_socketio_asgi.py`.
- Delete `src/python/tests/integration/test_websocket_handlers.py`.
- Delete `src/python/websocket/handlers.py`.
- Modify `src/python/websocket/__init__.py`: keep it as a package marker while `event_loop_bridge.py` remains in the package.
- Keep `src/python/websocket/event_loop_bridge.py` if process-manager callbacks still need to schedule coroutines onto the app event loop.
- Create `src/web/realtime/useCalculationUpdateStream.ts`: SSE connection lifecycle.
- Create `src/web/realtime/useCalculationSync.ts`: transport-independent update handling.
- Move or recreate `src/web/websocket/useCalculationNotifier.ts` as `src/web/realtime/useCalculationNotifier.ts`.
- Move or recreate `src/web/websocket/invalidateQueriesWithRetry.ts` as `src/web/realtime/invalidateQueriesWithRetry.ts`.
- Create `src/web/hooks/useCalculationUpdates.ts`: public hook used by `App.tsx`.
- Modify `src/web/App.tsx`: use `useCalculationUpdates`.
- Delete `src/web/hooks/useUnifiedWebSocket.ts`.
- Delete `src/web/websocket/useSocketTransport.ts`, `src/web/websocket/useWebSocketConnection.ts`, and old `src/web/websocket/useCalculationSync.ts`.
- Modify `package.json` and `package-lock.json`: remove `socket.io-client` and update `verify-build-env`.
- Modify `.github/environment.yml`: remove `python-socketio`.
- Regenerate `.github/pyscf-env.conda-lock.yml` and `.github/conda-locks/*` with `npm run conda-lock:generate`.
- Modify `config/server-config.json`: remove `socketio` config block.
- Modify docs/tests references that mention Socket.IO.

## Task 1: OpenAPI Contract And Codegen

**Files:**
- Modify: `src/api-spec/openapi.yaml`
- Regenerate: `src/python/generated_models.py`
- Regenerate: `src/web/types/generated-api.ts`

- [ ] **Step 1: Add SSE path entries to OpenAPI**

Add these paths near existing `/api/quantum/calculations` paths:

```yaml
  /api/quantum/calculations/updates/stream:
    parameters:
      - $ref: '#/components/parameters/AuthTokenHeader'
    get:
      security:
        - AuthToken: []
      tags:
        - Quantum Calculation
      summary: Stream calculation updates
      operationId: streamCalculationUpdates
      description: |
        Opens a Server-Sent Events stream for all calculation updates.
        Each event line contains a JSON `CalculationUpdateStreamEvent`.
      responses:
        '200':
          description: Successfully opened a calculation update SSE stream.
          content:
            text/event-stream:
              schema:
                type: string
        '401':
          $ref: '#/components/responses/UnauthorizedError'

  /api/quantum/calculations/{calculationId}/updates/stream:
    parameters:
      - $ref: '#/components/parameters/AuthTokenHeader'
    get:
      security:
        - AuthToken: []
      tags:
        - Quantum Calculation
      summary: Stream updates for one calculation
      operationId: streamCalculationUpdatesForCalculation
      parameters:
        - name: calculationId
          in: path
          required: true
          schema:
            type: string
      description: |
        Opens a Server-Sent Events stream for one calculation and emits the
        current calculation payload before future updates.
      responses:
        '200':
          description: Successfully opened a calculation-specific SSE stream.
          content:
            text/event-stream:
              schema:
                type: string
        '401':
          $ref: '#/components/responses/UnauthorizedError'
```

- [ ] **Step 2: Add stream event schemas for documentation**

Add these schemas in `components.schemas`:

```yaml
    CalculationUpdateStreamEvent:
      type: object
      required:
        - type
        - payload
      properties:
        type:
          type: string
          enum:
            - calculation_update
            - heartbeat
            - error
        payload:
          type: object
          additionalProperties: true

    CalculationUpdateStreamPayload:
      type: object
      required:
        - calculation
      properties:
        calculation:
          $ref: '#/components/schemas/CalculationInstance'

    CalculationUpdateStreamErrorPayload:
      type: object
      required:
        - message
      properties:
        message:
          type: string
        calculation_id:
          type: string
```

- [ ] **Step 3: Run codegen**

Run:

```bash
npm run codegen
```

Expected: `src/python/generated_models.py` and `src/web/types/generated-api.ts` are regenerated without manual edits.

- [ ] **Step 4: Verify generated file drift is limited**

Run:

```bash
git diff -- src/api-spec/openapi.yaml src/python/generated_models.py src/web/types/generated-api.ts
```

Expected: only the new SSE path/schema and generated type changes appear.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add src/api-spec/openapi.yaml src/python/generated_models.py src/web/types/generated-api.ts
git commit -m "feat: add calculation update SSE contract"
```

## Task 2: Backend SSE Hub

**Files:**
- Create: `src/python/services/calculation_update_stream.py`
- Test: `src/python/tests/unit/test_services/test_calculation_update_stream.py`

- [ ] **Step 1: Write hub tests**

Create `src/python/tests/unit/test_services/test_calculation_update_stream.py`:

```python
import asyncio

import pytest

from services.calculation_update_stream import (
    CalculationUpdateStreamHub,
    format_sse_event,
)


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
async def test_unsubscribe_removes_subscriber() -> None:
    hub = CalculationUpdateStreamHub()
    subscriber = hub.subscribe_global()

    hub.unsubscribe(subscriber)
    await hub.publish_calculation_update(
        {"id": "calc-1", "status": "completed"}
    )

    assert subscriber.queue.empty()
```

- [ ] **Step 2: Run tests and confirm failure**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_services/test_calculation_update_stream.py -v
```

Expected: fails with `ModuleNotFoundError: No module named 'services.calculation_update_stream'`.

- [ ] **Step 3: Implement the hub**

Create `src/python/services/calculation_update_stream.py`:

```python
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
```

- [ ] **Step 4: Run hub tests**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_services/test_calculation_update_stream.py -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit Task 2**

Run:

```bash
git add src/python/services/calculation_update_stream.py src/python/tests/unit/test_services/test_calculation_update_stream.py
git commit -m "feat: add calculation update SSE hub"
```

## Task 3: Backend SSE Routes

**Files:**
- Create: `src/python/api/calculation_updates.py`
- Create: `src/python/services/calculation_update_payload.py`
- Modify: `src/python/api/__init__.py`
- Test: `src/python/tests/integration/test_api_endpoints/test_calculation_updates_api.py`

- [ ] **Step 1: Write API tests**

Create `src/python/tests/integration/test_api_endpoints/test_calculation_updates_api.py`:

```python
import json
from pathlib import Path

from quantum_calc import get_current_settings
from services.calculation_update_stream import get_calculation_update_stream_hub


def _extract_sse_payload(line: str) -> dict:
    assert line.startswith("data: ")
    return json.loads(line.removeprefix("data: "))


def test_global_stream_returns_event_stream_content_type(client) -> None:
    with client.stream("GET", "/api/quantum/calculations/updates/stream") as response:
        assert response.status_code == 200
        assert response.headers["content-type"].startswith("text/event-stream")


def test_calculation_stream_emits_initial_payload(client) -> None:
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
        json.dumps({"status": "running"}),
        encoding="utf-8",
    )

    with client.stream(
        "GET",
        f"/api/quantum/calculations/{calc_id}/updates/stream",
    ) as response:
        assert response.status_code == 200
        first_line = next(response.iter_lines())

    event = _extract_sse_payload(first_line)
    assert event["type"] == "calculation_update"
    assert event["payload"]["calculation"]["id"] == calc_id
    assert event["payload"]["calculation"]["status"] == "running"


def test_calculation_stream_missing_calculation_emits_error(client) -> None:
    with client.stream(
        "GET",
        "/api/quantum/calculations/missing-calc/updates/stream",
    ) as response:
        assert response.status_code == 200
        first_line = next(response.iter_lines())

    event = _extract_sse_payload(first_line)
    assert event["type"] == "error"
    assert event["payload"]["calculation_id"] == "missing-calc"
    assert "not found" in event["payload"]["message"]


def test_global_stream_receives_published_update(client) -> None:
    hub = get_calculation_update_stream_hub()

    with client.stream("GET", "/api/quantum/calculations/updates/stream") as response:
        assert response.status_code == 200
        client.portal.call(
            hub.publish_calculation_update,
            {"id": "calc-global", "status": "completed"},
        )
        first_line = next(response.iter_lines())

    event = _extract_sse_payload(first_line)
    assert event["type"] == "calculation_update"
    assert event["payload"]["calculation"]["id"] == "calc-global"
```

- [ ] **Step 2: Run tests and confirm failure**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_calculation_updates_api.py -v
```

Expected: route tests fail because `api.calculation_updates` is not registered.

- [ ] **Step 3: Create the shared calculation payload builder**

Create `src/python/services/calculation_update_payload.py`:

```python
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
    status = status_override or repository.read_calculation_status(calc_dir)
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

    if error_message:
        instance["error"] = error_message
        instance["errorMessage"] = error_message

    return instance
```

- [ ] **Step 4: Implement streaming routes**

Create `src/python/api/calculation_updates.py`:

```python
"""SSE routes for calculation update monitoring."""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator
from pathlib import Path

from fastapi import APIRouter, Request
from fastapi.responses import StreamingResponse

from quantum_calc import CalculationRepository, get_current_settings
from services.calculation_update_stream import (
    StreamSubscriber,
    format_sse_event,
    get_calculation_update_stream_hub,
)
from services.calculation_update_payload import build_calculation_instance

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

    subscriber = get_calculation_update_stream_hub().subscribe_calculation(
        calculation_id
    )
    initial_calculation = build_calculation_instance(
        calculation_id,
        calc_path,
        repository,
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
```

- [ ] **Step 5: Register the router**

Modify `src/python/api/__init__.py`:

```python
from .calculation_updates import router as calculation_updates_router
```

Add registration after `quantum_router`:

```python
    app.include_router(quantum_router)
    app.include_router(calculation_updates_router)
```

- [ ] **Step 6: Run route tests**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_calculation_updates_api.py -v
```

Expected: all SSE route tests pass.

- [ ] **Step 7: Commit Task 3**

Run:

```bash
git add src/python/api/calculation_updates.py src/python/services/calculation_update_payload.py src/python/api/__init__.py src/python/tests/integration/test_api_endpoints/test_calculation_updates_api.py
git commit -m "feat: add calculation update SSE routes"
```

## Task 4: Notification Service And FastAPI App Cleanup

**Files:**
- Modify: `src/python/services/notification_service.py`
- Modify: `src/python/app.py`
- Modify: `src/python/tests/conftest.py`
- Create: `src/python/tests/unit/test_services/test_notification_service.py`

- [ ] **Step 1: Write notification service test**

Create `src/python/tests/unit/test_services/test_notification_service.py`:

```python
import asyncio
import json

import pytest

from services.calculation_update_stream import get_calculation_update_stream_hub
from services.notification_service import NotificationService
from websocket.event_loop_bridge import bind_event_loop, clear_event_loop


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
```

- [ ] **Step 2: Run notification test and confirm failure**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_services/test_notification_service.py -v
```

Expected: fails because `NotificationService` still requires Socket.IO binding.

- [ ] **Step 3: Refactor notification service**

Replace `src/python/services/notification_service.py` Socket.IO binding with
SSE publishing:

```python
"""Calculation update notification service."""

import logging
import os
from typing import Optional

from websocket.event_loop_bridge import schedule_coroutine

logger = logging.getLogger(__name__)


class NotificationService:
    """Service for publishing calculation update notifications."""

    def send_calculation_update(
        self,
        calculation_id: str,
        status: str,
        error_message: Optional[str] = None,
    ) -> None:
        """Publish immediate notification for calculation status changes."""
        try:
            from quantum_calc import CalculationRepository, get_current_settings
            from services.calculation_update_stream import (
                get_calculation_update_stream_hub,
            )
            from services.calculation_update_payload import (
                build_calculation_instance,
            )

            settings = get_current_settings()
            file_manager = CalculationRepository(
                base_dir=settings.calculations_directory
            )
            calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)

            if not os.path.exists(calc_dir):
                logger.warning("Calculation directory not found: %s", calc_dir)
                return

            calculation_instance = build_calculation_instance(
                calculation_id,
                calc_dir,
                file_manager,
                status_override=status,
                error_message=error_message,
            )

        except Exception:
            logger.exception("Error building calculation notification payload")
            return

        scheduled_future = schedule_coroutine(
            get_calculation_update_stream_hub().publish_calculation_update(
                calculation_instance
            )
        )
        if scheduled_future is None:
            logger.warning(
                "Dropped calculation notification for %s because scheduling failed",
                calculation_id,
            )
            return

        logger.debug(
            "Scheduled calculation notification for %s with status %s",
            calculation_id,
            status,
        )
```

Keep the existing singleton:

```python
_notification_service: Optional[NotificationService] = None


def get_notification_service() -> NotificationService:
    global _notification_service
    if _notification_service is None:
        _notification_service = NotificationService()
    return _notification_service
```

Remove `bind_socketio` and `bind_notification_service`.

- [ ] **Step 4: Simplify FastAPI app**

Modify `src/python/app.py`:

```python
import asyncio
import json
import logging
from contextlib import asynccontextmanager
from typing import Any

from fastapi import FastAPI, Request
```

Remove:

```python
import socketio
from websocket import register_websocket_handlers
```

Remove `create_socketio_server`, `compose_asgi_app`, `sio`, and Socket.IO
binding. Keep:

```python
def create_app(
    server_port: int | None = None,
    test_config: dict[str, Any] | None = None,
) -> FastAPI:
    return create_fastapi_app(server_port=server_port, test_config=test_config)


app = create_app()
```

- [ ] **Step 5: Update test fixtures**

Modify `src/python/tests/conftest.py`:

- `app` fixture no longer sets a `SOCKETIO` test config block.
- Keep fixture name `asgi_server` for existing integration tests, but update it
  to return a Uvicorn URL for the plain FastAPI `create_app`.
- Remove imports of `services.notification_service`, `websocket.handlers._sid_state`,
  and Socket.IO-specific cleanup.

The server fixture should run:

```python
from app import create_app

test_app = create_app(server_port=port, test_config=test_config)
server = uvicorn.Server(
    uvicorn.Config(test_app, host="127.0.0.1", port=port, log_level="warning")
)
```

- [ ] **Step 6: Run backend app tests**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_app_startup_fastapi.py tests/unit/test_services/test_notification_service.py tests/integration/test_api_endpoints/test_calculation_updates_api.py -v
```

Expected: all tests pass.

- [ ] **Step 7: Commit Task 4**

Run:

```bash
git add src/python/app.py src/python/services/notification_service.py src/python/tests/conftest.py src/python/tests/unit/test_services/test_notification_service.py
git commit -m "feat: publish calculation updates over SSE"
```

## Task 5: Frontend SSE Transport

**Files:**
- Create: `src/web/realtime/useCalculationUpdateStream.ts`
- Create: `src/web/realtime/useCalculationSync.ts`
- Create: `src/web/realtime/useCalculationNotifier.ts`
- Create: `src/web/realtime/invalidateQueriesWithRetry.ts`
- Create: `src/web/hooks/useCalculationUpdates.ts`
- Modify: `src/web/App.tsx`
- Test: `src/web/realtime/useCalculationUpdateStream.test.ts`

- [ ] **Step 1: Write SSE client tests**

Create `src/web/realtime/useCalculationUpdateStream.test.ts`:

```typescript
import assert from 'assert';
import { parseCalculationUpdateStreamEvent } from './useCalculationUpdateStream';

const updateEvent = JSON.stringify({
  type: 'calculation_update',
  payload: {
    calculation: {
      id: 'calc-1',
      name: 'SSE Calc',
      status: 'completed',
      createdAt: '2026-05-26T00:00:00',
      updatedAt: '2026-05-26T00:00:01',
      parameters: {},
      results: {},
    },
  },
});

const parsed = parseCalculationUpdateStreamEvent(updateEvent);
assert.strictEqual(parsed.type, 'calculation_update');
assert.strictEqual(parsed.calculation?.id, 'calc-1');

const heartbeat = parseCalculationUpdateStreamEvent(
  JSON.stringify({ type: 'heartbeat', payload: {} })
);
assert.strictEqual(heartbeat.type, 'heartbeat');
assert.strictEqual(heartbeat.calculation, undefined);

const error = parseCalculationUpdateStreamEvent(
  JSON.stringify({ type: 'error', payload: { message: 'Missing calculation' } })
);
assert.strictEqual(error.type, 'error');
assert.strictEqual(error.errorMessage, 'Missing calculation');

console.log('calculation update SSE parser tests passed');
```

- [ ] **Step 2: Add test script entry**

Modify `package.json` `test:web` to include the new parser test:

```json
"test:web": "ts-node --files src/web/api/agent.test.ts && ts-node --files src/web/runtime-config.test.ts && ts-node --files src/web/realtime/useCalculationUpdateStream.test.ts && ts-node --files src/splash/splash-html.test.ts"
```

- [ ] **Step 3: Run parser test and confirm failure**

Run:

```bash
npm run test:web
```

Expected: fails because `src/web/realtime/useCalculationUpdateStream.ts` does not exist.

- [ ] **Step 4: Create transport-independent sync hook**

Create `src/web/realtime/useCalculationSync.ts` by moving logic from
`src/web/websocket/useCalculationSync.ts` and changing its public return shape:

```typescript
export interface UseCalculationSyncOptions {
  activeCalculationId: string | null;
  onCalculationUpdate?: (
    updated: CalculationInstance,
    previousStatus: string | undefined
  ) => void;
  onStreamError?: (error: string) => void;
}

export interface UseCalculationSyncReturn {
  handleCalculationUpdate: (updated: CalculationInstance) => void;
  handleStreamError: (errorData: unknown) => void;
  handleReconnect: () => Promise<void>;
}
```

Remove imports of `Socket` and `SocketEventHandlers`. Keep:

```typescript
const updateCalculationCache = useCallback(
  (updatedCalculation: CalculationInstance) => {
    const calculationId = updatedCalculation.id;
    queryClient.invalidateQueries({
      queryKey: calculationQueryKeys.detail(calculationId),
    });
    queryClient.invalidateQueries({
      queryKey: calculationQueryKeys.list(),
    });
  },
  [queryClient]
);
```

Keep duplicate suppression by `updatedAt` and existing `onCalculationUpdate`
callback behavior.

- [ ] **Step 5: Move notifier and retry helper**

Copy existing code without behavior changes:

```bash
mkdir -p src/web/realtime
cp src/web/websocket/useCalculationNotifier.ts src/web/realtime/useCalculationNotifier.ts
cp src/web/websocket/invalidateQueriesWithRetry.ts src/web/realtime/invalidateQueriesWithRetry.ts
```

Then update relative imports in `src/web/realtime/useCalculationSync.ts` to use
`./invalidateQueriesWithRetry`. In `src/web/realtime/useCalculationNotifier.ts`,
rename `handleWebSocketError` to `handleStreamError` and change log prefixes
from `[UnifiedWebSocket]` to `[CalculationUpdates]`.

- [ ] **Step 6: Implement SSE hook**

Create `src/web/realtime/useCalculationUpdateStream.ts`:

```typescript
import { useCallback, useEffect, useRef, useState } from 'react';
import type { MutableRefObject } from 'react';
import { fetchEventSource } from '@microsoft/fetch-event-source';
import { getApiBaseUrl } from '../api/core';
import type { CalculationInstance } from '../types/api-types';

const STREAM_RETRY_INTERVAL_MS = 2000;

export type CalculationUpdateStreamEventType =
  | 'calculation_update'
  | 'heartbeat'
  | 'error';

export interface ParsedCalculationUpdateStreamEvent {
  type: CalculationUpdateStreamEventType;
  calculation?: CalculationInstance;
  errorMessage?: string;
}

export function parseCalculationUpdateStreamEvent(
  data: string
): ParsedCalculationUpdateStreamEvent {
  const parsed = JSON.parse(data);
  if (parsed.type === 'calculation_update') {
    return {
      type: 'calculation_update',
      calculation: parsed.payload?.calculation,
    };
  }
  if (parsed.type === 'error') {
    return {
      type: 'error',
      errorMessage:
        parsed.payload?.message || 'A calculation update stream error occurred.',
    };
  }
  return { type: 'heartbeat' };
}

export interface UseCalculationUpdateStreamOptions {
  activeCalculationId: string | null;
  onCalculationUpdate: (calculation: CalculationInstance) => void;
  onStreamError: (error: string) => void;
  onReconnect: () => Promise<void>;
}

export function useCalculationUpdateStream({
  activeCalculationId,
  onCalculationUpdate,
  onStreamError,
  onReconnect,
}: UseCalculationUpdateStreamOptions) {
  const [isConnected, setIsConnected] = useState(false);
  const globalAbortRef = useRef<AbortController | null>(null);
  const detailAbortRef = useRef<AbortController | null>(null);
  const hasGlobalConnectedRef = useRef(false);

  const startStream = useCallback(
    async (
      path: string,
      abortRef: MutableRefObject<AbortController | null>,
      trackGlobalConnection = false
    ) => {
      abortRef.current?.abort();
      const controller = new AbortController();
      abortRef.current = controller;
      const authToken = await window.electronAPI?.getAuthToken?.();
      const headers: HeadersInit = { Accept: 'text/event-stream' };
      if (authToken) {
        (headers as Record<string, string>)['X-Auth-Token'] = authToken;
      }

      await fetchEventSource(`${getApiBaseUrl()}${path}`, {
        method: 'GET',
        headers,
        signal: controller.signal,
        onopen: async response => {
          if (!response.ok) {
            throw new Error(`Stream failed: ${response.status}`);
          }
          if (trackGlobalConnection) {
            const wasConnected = hasGlobalConnectedRef.current;
            hasGlobalConnectedRef.current = true;
            setIsConnected(true);
            if (wasConnected) {
              await onReconnect();
            }
          }
        },
        onmessage: event => {
          const parsed = parseCalculationUpdateStreamEvent(event.data);
          if (parsed.type === 'calculation_update' && parsed.calculation) {
            onCalculationUpdate(parsed.calculation);
          }
          if (parsed.type === 'error') {
            onStreamError(parsed.errorMessage || 'Calculation stream error.');
          }
        },
        onclose: () => {
          if (trackGlobalConnection) {
            setIsConnected(false);
          }
        },
        onerror: error => {
          if (controller.signal.aborted) {
            return;
          }
          if (trackGlobalConnection) {
            setIsConnected(false);
          }
          onStreamError(error instanceof Error ? error.message : String(error));
          return STREAM_RETRY_INTERVAL_MS;
        },
      });
    },
    [onCalculationUpdate, onReconnect, onStreamError]
  );

  useEffect(() => {
    void startStream(
      '/api/quantum/calculations/updates/stream',
      globalAbortRef,
      true
    );
    return () => globalAbortRef.current?.abort();
  }, [startStream]);

  useEffect(() => {
    detailAbortRef.current?.abort();
    if (!activeCalculationId || activeCalculationId.startsWith('new-calculation-')) {
      return;
    }
    void startStream(
      `/api/quantum/calculations/${encodeURIComponent(activeCalculationId)}/updates/stream`,
      detailAbortRef,
      false
    );
    return () => detailAbortRef.current?.abort();
  }, [activeCalculationId, startStream]);

  return {
    isConnected,
    reconnect: () => {
      void startStream(
        '/api/quantum/calculations/updates/stream',
        globalAbortRef,
        true
      );
    },
    disconnect: () => {
      globalAbortRef.current?.abort();
      detailAbortRef.current?.abort();
      hasGlobalConnectedRef.current = false;
      setIsConnected(false);
    },
  };
}
```

- [ ] **Step 7: Create public hook**

Create `src/web/hooks/useCalculationUpdates.ts`:

```typescript
import { useCalculationNotifier } from '../realtime/useCalculationNotifier';
import { useCalculationSync } from '../realtime/useCalculationSync';
import { useCalculationUpdateStream } from '../realtime/useCalculationUpdateStream';

export interface UseCalculationUpdatesOptions {
  activeCalculationId: string | null;
}

export const useCalculationUpdates = ({
  activeCalculationId,
}: UseCalculationUpdatesOptions) => {
  const { notifyCalculationUpdate, handleStreamError: notifyStreamError } =
    useCalculationNotifier(activeCalculationId);
  const { handleCalculationUpdate, handleStreamError, handleReconnect } =
    useCalculationSync({
      activeCalculationId,
      onCalculationUpdate: notifyCalculationUpdate,
      onStreamError: notifyStreamError,
    });

  return useCalculationUpdateStream({
    activeCalculationId,
    onCalculationUpdate: handleCalculationUpdate,
    onStreamError: handleStreamError,
    onReconnect: handleReconnect,
  });
};
```

- [ ] **Step 8: Update App import and usage**

Modify `src/web/App.tsx`:

```typescript
import { useCalculationUpdates } from './hooks/useCalculationUpdates';
```

Replace:

```typescript
useUnifiedWebSocket({
  activeCalculationId,
});
```

With:

```typescript
useCalculationUpdates({
  activeCalculationId,
});
```

- [ ] **Step 9: Run frontend tests**

Run:

```bash
npm run test:web
npm run typecheck
```

Expected: parser tests pass and TypeScript has no Socket.IO type references.

- [ ] **Step 10: Commit Task 5**

Run:

```bash
git add package.json src/web/App.tsx src/web/hooks/useCalculationUpdates.ts src/web/realtime
git commit -m "feat: stream calculation updates with SSE"
```

## Task 6: Remove Socket.IO Code And Dependencies

**Files:**
- Delete: `src/python/websocket/handlers.py`
- Modify: `src/python/websocket/__init__.py`
- Delete: `src/python/tests/integration/test_socketio_asgi.py`
- Delete: `src/python/tests/integration/test_websocket_handlers.py`
- Delete: `src/web/hooks/useUnifiedWebSocket.ts`
- Delete: `src/web/websocket/useSocketTransport.ts`
- Delete: `src/web/websocket/useWebSocketConnection.ts`
- Delete: `src/web/websocket/useCalculationSync.ts`
- Delete after move: `src/web/websocket/useCalculationNotifier.ts`
- Delete after move: `src/web/websocket/invalidateQueriesWithRetry.ts`
- Modify: `package.json`
- Modify: `package-lock.json`
- Modify: `.github/environment.yml`
- Modify: `.github/pyscf-env.conda-lock.yml`
- Modify: `.github/conda-locks/pyscf-env-linux-64.lock`
- Modify: `.github/conda-locks/pyscf-env-osx-64.lock`
- Modify: `.github/conda-locks/pyscf-env-osx-arm64.lock`
- Modify: `config/server-config.json`
- Modify: `scripts/verify-environment.py`
- Modify: `scripts/test-python-standalone.js`
- Modify: `scripts/build-python-linux.sh`
- Modify: `src/python/pytest.ini`
- Modify: `src/python/tests/README.md`

- [ ] **Step 1: Remove Socket.IO source files**

Run:

```bash
git rm src/python/websocket/handlers.py
git rm src/python/tests/integration/test_socketio_asgi.py src/python/tests/integration/test_websocket_handlers.py
git rm src/web/hooks/useUnifiedWebSocket.ts
git rm src/web/websocket/useSocketTransport.ts src/web/websocket/useWebSocketConnection.ts src/web/websocket/useCalculationSync.ts
git rm src/web/websocket/useCalculationNotifier.ts src/web/websocket/invalidateQueriesWithRetry.ts
```

Replace `src/python/websocket/__init__.py` with a package marker while
`event_loop_bridge.py` is still imported:

```python
"""Event-loop bridge package for backend realtime publishing."""
```

Keep `src/python/websocket/event_loop_bridge.py` unless `rg "event_loop_bridge|schedule_coroutine|bind_event_loop|clear_event_loop" src/python` proves it is unused after notification refactor.

- [ ] **Step 2: Remove JS dependency**

Run:

```bash
npm uninstall socket.io-client
```

Expected: `socket.io-client` is removed from `package.json` and `package-lock.json`.

- [ ] **Step 3: Remove Python Socket.IO dependency**

Modify `.github/environment.yml` and remove:

```yaml
  - python-socketio=5.16.1
```

Modify `package.json` `verify-build-env` to remove `socketio` import:

```json
"verify-build-env": "bash -c 'CONDA_BASE=$(conda info --base 2>/dev/null || echo \"$HOME/miniforge3\") && source \"$CONDA_BASE/etc/profile.d/conda.sh\" && conda activate pyscf-env && echo \"=== Build Environment Verification ===\" && $(conda info --base)/envs/pyscf-env/bin/python -c \"import conda_pack; print(f\\\"conda-pack: {conda_pack.__version__}\\\")\" && $(conda info --base)/envs/pyscf-env/bin/python -c \"import gunicorn; print(f\\\"Gunicorn: {gunicorn.__version__}\\\")\" && $(conda info --base)/envs/pyscf-env/bin/python -c \"import fastapi, uvicorn; print(\\\"FastAPI ASGI stack verified\\\")\" && echo \"✓ All build tools verified\"'"
```

- [ ] **Step 4: Regenerate conda lock files**

Run:

```bash
npm run conda-lock:generate
```

Expected: `python-socketio` and `python-engineio` disappear from the generated lock outputs unless another dependency still requires them.

- [ ] **Step 5: Remove Socket.IO config and verification references**

Remove the `socketio` object from `config/server-config.json`.

Update `scripts/verify-environment.py`:

- Remove `socketio` from required packages.
- Replace "FastAPI ASGI + Socket.IO" test with a FastAPI `TestClient` smoke test only.

Update `scripts/test-python-standalone.js` FastAPI import command:

```javascript
'import fastapi, uvicorn; from uvicorn.workers import UvicornWorker; print("FastAPI ASGI import successful")'
```

Update `scripts/build-python-linux.sh` and remove the `python-socketio` check.

Update `src/python/pytest.ini` and keep the Flask plugin disable because ASE
still pulls Flask into the environment:

```ini
    # FastAPI tests define their own app/client fixtures; disable pytest-flask
    # when it is present indirectly through ASE in the conda environment.
    -p no:flask
```

Rename the realtime marker description:

```ini
    realtime: Realtime update stream tests
```

Update `src/python/tests/README.md` to describe SSE endpoint tests instead of Socket.IO tests.

- [ ] **Step 6: Search for stale references**

Run:

```bash
rg -n "socketio|Socket.IO|socket.io-client|join_calculation|join_global_updates|leave_calculation|calculation_update room|useUnifiedWebSocket|UnifiedWebSocket" src package.json package-lock.json .github config scripts
```

Expected: no references remain. Rename remaining frontend log prefixes from `[UnifiedWebSocket]` to `[CalculationUpdates]`.

- [ ] **Step 7: Run dependency and type verification**

Run:

```bash
npm run typecheck
npm run test:web
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_calculation_updates_api.py tests/unit/test_services/test_calculation_update_stream.py -v
```

Expected: all pass.

- [ ] **Step 8: Commit Task 6**

Run:

```bash
git add package.json package-lock.json .github/environment.yml .github/pyscf-env.conda-lock.yml .github/conda-locks config/server-config.json scripts src
git commit -m "chore: remove Socket.IO calculation monitoring"
```

## Task 7: Full Verification And Manual Smoke

**Files:**
- No planned file edits.

- [ ] **Step 1: Run codegen drift check**

Run:

```bash
npm run codegen:check
```

Expected: no generated-file drift.

- [ ] **Step 2: Run frontend checks**

Run:

```bash
npm run typecheck
npm run test:web
npm run test:main
```

Expected: all commands exit 0.

- [ ] **Step 3: Run backend focused tests**

Run:

```bash
cd src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_services/test_calculation_update_stream.py tests/unit/test_services/test_notification_service.py tests/integration/test_api_endpoints/test_calculation_updates_api.py tests/integration/test_app_startup_fastapi.py -v
```

Expected: all tests pass.

- [ ] **Step 4: Run environment verification**

Run:

```bash
npm run verify-env
npm run verify-build-env
```

Expected: both commands exit 0 and no longer import `socketio`.

- [ ] **Step 5: Run stale reference search**

Run:

```bash
rg -n "socketio|Socket.IO|socket.io-client|python-socketio|python-engineio|join_calculation|join_global_updates|leave_calculation|useUnifiedWebSocket|UnifiedWebSocket" src package.json package-lock.json .github config scripts
```

Expected: no implementation references remain.

- [ ] **Step 6: Manual dev smoke**

Run:

```bash
npm run dev
```

Manual checks:

- Open the app.
- Start a small HF or DFT calculation.
- Confirm the calculation list updates without page refresh.
- Open the calculation results page while it is running.
- Confirm the active calculation detail updates.
- Stop the dev server cleanly.

- [ ] **Step 7: Confirm working tree is clean**

Run:

```bash
git status --short
```

Expected: no uncommitted implementation changes remain after Task 6. If Step 5
found stale references, edit the exact files reported by `rg`, rerun Step 5,
and commit those concrete files with `git commit -m "test: verify SSE calculation monitoring"`.

## Self-Review Checklist

- Spec coverage: Tasks cover OpenAPI, generated artifacts, backend hub, backend routes, notification publishing, FastAPI app simplification, frontend SSE transport, Socket.IO deletion, dependency cleanup, tests, and manual smoke.
- Placeholder scan: This plan contains no unresolved placeholders.
- Type consistency: Backend event type is `calculation_update`; frontend parser expects `payload.calculation`; OpenAPI documents `text/event-stream` responses as strings plus schemas for event payload documentation.
- Risk note: `src/python/websocket/event_loop_bridge.py` is intentionally retained unless later search proves it unused, because process-manager callbacks may still need the app event loop bridge.
