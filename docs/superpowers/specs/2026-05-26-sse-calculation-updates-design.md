# SSE Calculation Updates Design

## Goal

Replace Socket.IO-based calculation monitoring with FastAPI Server-Sent Events
while preserving the user-visible behavior of real-time calculation updates.

## Current State

The backend currently composes a FastAPI app with `socketio.ASGIApp` in
`src/python/app.py`. Socket.IO handlers in `src/python/websocket/handlers.py`
provide `join_global_updates`, `join_calculation`, `leave_calculation`, and
`calculation_update`. The frontend uses `socket.io-client` through
`src/web/websocket/useSocketTransport.ts` and `src/web/websocket/useCalculationSync.ts`.

The agent chat already uses SSE through `@microsoft/fetch-event-source`, so the
frontend has a proven streaming client dependency.

## Scope

In scope:

- Add FastAPI SSE endpoints for calculation update streams.
- Replace the frontend Socket.IO transport with an SSE transport.
- Keep existing calculation cache invalidation and notification behavior.
- Remove Socket.IO backend composition, handlers, tests, config, and frontend
  dependency after the SSE path is verified.
- Update OpenAPI and generated Python/TypeScript API artifacts.

Out of scope:

- Changing calculation status semantics.
- Changing REST calculation APIs.
- Changing AI chat SSE behavior.
- Adding native FastAPI WebSocket endpoints.
- Keeping a Socket.IO compatibility layer.

## API Design

Add two SSE endpoints:

- `GET /api/quantum/calculations/updates/stream`
  Streams all calculation updates that were previously emitted to the
  `global_updates` Socket.IO room.

- `GET /api/quantum/calculations/{calculationId}/updates/stream`
  Streams updates for one calculation and sends an initial calculation payload
  when the calculation exists.

Both endpoints return `text/event-stream`.

Event payloads use the same envelope shape as agent chat SSE:

```json
{"type":"calculation_update","payload":{"calculation":{"id":"...","status":"running"}}}
```

Additional event types:

- `heartbeat`: keeps long-lived connections alive.
- `error`: reports stream-level errors such as invalid calculation IDs.

Authentication uses the same `X-Auth-Token` header as existing HTTP APIs. The
frontend uses `fetchEventSource`, not browser `EventSource`, because auth
headers are required in packaged Electron.

## Backend Design

Create a small SSE hub responsible for subscriber management:

- `src/python/services/calculation_update_stream.py`
  Owns async queues for global subscribers and per-calculation subscribers.
  Provides `publish_calculation_update`, `subscribe_global`, and
  `subscribe_calculation`.

Refactor notification delivery:

- `src/python/services/notification_service.py`
  Builds the existing calculation instance payload, then publishes to the SSE
  hub instead of calling `Socket.IO.emit`.

Add FastAPI route handlers:

- `src/python/api/calculation_updates.py`
  Exposes the two streaming endpoints with `StreamingResponse`.
  Handles disconnects by unregistering subscribers.
  Emits heartbeats periodically.

Simplify app startup:

- `src/python/app.py`
  Returns a plain FastAPI app as `app`.
  Removes `socketio.ASGIApp`, `create_socketio_server`,
  `compose_asgi_app`, and Socket.IO binding.

## Frontend Design

Replace the Socket.IO transport with an SSE monitor:

- `src/web/websocket/useCalculationSse.ts`
  Opens the global stream once while the app is mounted.
  Opens a calculation-specific stream when `activeCalculationId` is a persisted
  calculation ID.
  Aborts and recreates the calculation-specific stream when active calculation
  changes.

Keep update processing:

- Preserve the existing duplicate suppression by `updatedAt`.
- Preserve React Query invalidation for calculation detail and list queries.
- Preserve user notifications for status changes and stream errors.

Remove Socket.IO-specific files after replacement:

- `src/web/websocket/useSocketTransport.ts`
- Socket.IO event-handler types and `emit` room management
- `socket.io-client` dependency

## Data Flow

1. A calculation changes status or result data.
2. `NotificationService.send_calculation_update` builds a calculation payload.
3. The SSE hub publishes the payload to global subscribers and matching
   calculation-specific subscribers.
4. The frontend receives `calculation_update`.
5. Existing cache invalidation and notification logic updates the UI.
6. If an SSE stream reconnects, the frontend invalidates the calculation list
   and active detail query to recover missed events.

## Error Handling

- Invalid calculation IDs return an SSE `error` event and close the stream.
- Missing calculations return an SSE `error` event and close the stream.
- Client disconnect removes its queue from the hub.
- Slow subscribers use bounded queues; if a queue fills, the newest update
  replaces stale pending data for that subscriber.
- Heartbeats prevent idle connections from looking dead to intermediaries.

## Testing

Backend tests:

- SSE endpoint rejects missing or invalid auth in production mode.
- Global stream receives a published `calculation_update`.
- Calculation-specific stream emits an initial payload for existing
  calculations.
- Calculation-specific stream reports not found for missing calculations.
- Subscriber cleanup happens on disconnect.

Frontend tests:

- SSE parser routes `calculation_update` to the existing update handler.
- Aborting a stream suppresses expected abort errors.
- Active calculation changes close the previous detail stream and open the next.

Migration tests:

- Remove Socket.IO ASGI smoke tests and handler tests after equivalent SSE tests
  exist.
- Run `npm run codegen`, `npm run typecheck`, `npm run test:web`,
  `npm run test:main`, and focused Python tests.

## Migration Order

1. Add OpenAPI paths for both SSE endpoints and regenerate API artifacts.
2. Add backend SSE hub and streaming routes with tests.
3. Refactor `NotificationService` to publish to the SSE hub.
4. Add frontend SSE transport while keeping update-processing logic.
5. Remove Socket.IO backend composition, handlers, tests, config, and frontend
   dependency.
6. Run focused verification and a manual app smoke test with an active
   calculation.

## Accepted Tradeoffs

SSE is one-way. That is acceptable because calculation monitoring only needs
server-to-client update delivery. Client commands such as starting, pausing,
resuming, and deleting calculations remain REST APIs.

SSE does not provide Socket.IO rooms. URL-level stream selection replaces rooms:
one global stream for list updates and one optional detail stream for the active
calculation.
