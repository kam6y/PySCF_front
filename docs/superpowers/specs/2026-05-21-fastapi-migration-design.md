# FastAPI Migration Design

Date: 2026-05-21
Branch/worktree: `migration-fastapi`

## Goal

Replace the Python backend framework from Flask/Flask-SocketIO to a simple
FastAPI/ASGI stack while preserving the existing frontend API contract, SSE wire
format, Socket.IO event names, quantum-calculation service behavior, and Electron
startup flow from the user's perspective.

## Approved Approach

Use a full ASGI migration:

- FastAPI for HTTP APIs.
- FastAPI `APIRouter` modules instead of Flask `Blueprint` modules.
- FastAPI `StreamingResponse` for `/api/agent/chat` SSE.
- `python-socketio.AsyncServer(async_mode="asgi")` plus `socketio.ASGIApp` for
  Socket.IO traffic.
- Gunicorn with `uvicorn.workers.UvicornWorker` for Electron-managed production
  startup.
- No Flask compatibility adapter and no mixed Flask/FastAPI runtime.

This was chosen over a staged HTTP-only migration because the project is still in
development and the repo policy prefers simple breaking changes over long-lived
compatibility layers.

## Current Evidence

Baseline checks in the `migration-fastapi` worktree:

- `npm run typecheck`: passed.
- `~/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests/integration/test_api_endpoints/test_health_api.py src/python/tests/integration/test_auth_security.py src/python/tests/integration/test_auth_production.py src/python/tests/integration/test_websocket_handlers.py -v`: 24 passed, 3 skipped.

Existing dependency state before migration:

- `fastapi`: missing.
- `uvicorn`: missing.
- `starlette`: missing.
- `pydantic`: installed, version 2.11.7.
- `python-socketio`: installed and can construct ASGI apps.
- `flask`, `flask-cors`, `flask-socketio`, `flask-pydantic`: installed.

Primary references checked:

- FastAPI documentation through Context7 for `APIRouter`, dependencies,
  exception handlers, `StreamingResponse`, lifespan, and `TestClient`.
- python-socketio documentation through Context7 for `AsyncServer`,
  `ASGIApp`, ASGI deployment, rooms, and emits.

## Architecture

`src/python/app.py` becomes the ASGI entry point. It should expose:

- `fastapi_app`: the FastAPI application, useful for HTTP tests.
- `sio`: the `socketio.AsyncServer` instance, useful for Socket.IO handlers.
- `app`: the top-level ASGI application, built with `socketio.ASGIApp(sio, fastapi_app)`.

The numerical thread-control block that currently starts `src/python/app.py`
must remain the first executable code in the ASGI entry point. Set defaults for
`OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`BLIS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, and `NUMEXPR_NUM_THREADS` before
importing FastAPI routers, services, PySCF, NumPy, or any module that can import
numerical libraries. This preserves the current BLAS/PySCF startup invariant.

The application factory should:

1. Load `ServerConfig`.
2. Determine the server port from `PYSCF_SERVER_PORT`, CLI argument, or config.
3. Store framework-neutral settings on `fastapi_app.state`.
4. Register CORS/auth middleware and exception handlers.
5. Include all API routers.
6. Register Socket.IO handlers.
7. Bind the notification service.
8. Initialize the process manager callback.
9. Register lifespan cleanup for process manager and websocket watcher shutdown.

Create the FastAPI app with default documentation routes disabled:

```python
FastAPI(docs_url=None, redoc_url=None, openapi_url=None, lifespan=lifespan)
```

The existing development-only `/api-docs` and `/api-docs/spec.json` endpoints
remain the only documentation surface.

## HTTP API Design

All existing public paths remain unchanged:

- `/health`
- `/api/pubchem/*`
- `/api/smiles/*`
- `/api/settings`
- `/api/system/*`
- `/api/debug/system-diagnostics`
- `/api/debug/process-manager-diagnostics`
- `/api/debug/resource-manager-diagnostics`
- `/api/quantum/*`
- `/api/agent/chat`
- `/api/chat-history/*`
- Development-only `/api-docs` and `/api-docs/spec.json`

Response envelopes stay unchanged where existing frontend code expects them:

```json
{"success": true, "data": {}}
```

```json
{"success": false, "error": "message"}
```

FastAPI route functions should use generated Pydantic models from
`src/python/generated_models.py` for request body validation where the Flask
routes currently use `flask_pydantic`.

The existing OpenAPI source of truth remains `src/api-spec/openapi.yaml`.
FastAPI's auto-generated OpenAPI is not the contract source for this app.

FastAPI validation and parse errors must be normalized to the existing envelope
contract. Add exception handlers for:

- `fastapi.exceptions.RequestValidationError`
- malformed JSON/body parse errors
- Pydantic `ValidationError` raised by manual model validation

These handlers should return HTTP `400` with:

```json
{"success": false, "error": "Validation failed: ..."}
```

Do not allow FastAPI's default `422` validation response to leak through public
API endpoints. Add regression tests for missing required fields, malformed JSON,
and invalid query parameter types.

`POST /api/quantum/calculate` must preserve its current raw-body validation
order:

1. Read the raw JSON body.
2. Return `400` if the body is missing.
3. Extract and require `calculation_method`.
4. Run `validate_parameters_for_method(calculation_method, raw_data)` before
   Pydantic model validation so inapplicable method parameters cannot be silently
   ignored.
5. Validate with `QuantumCalculationRequest`.

Cover this order with regression tests for missing `calculation_method` and an
inapplicable parameter that Pydantic would otherwise ignore.

## Authentication

Replace Flask `before_request` with a single FastAPI HTTP middleware.

Behavior must match the current contract:

- In development, `/api-docs` and `/api-docs/*` are accessible without
  `X-Auth-Token`.
- `OPTIONS` requests are allowed without a token.
- If `PYSCF_AUTH_TOKEN` is set, all protected HTTP requests must provide an
  exact `X-Auth-Token` match.
- If no token is set and `PYSCF_ENV=production`, protected HTTP requests fail
  with `401`.
- In test/development mode without a token, protected requests are allowed with a
  warning-level log.

HTTP CORS behavior must remain compatible with the Electron renderer calling the
local backend from `file://` and `http://127.0.0.1:<port>`. Add
`CORSMiddleware` with:

- allowed localhost/loopback origins for HTTP and WebSocket development traffic,
- support for `file://` origins where Starlette can represent them,
- all current HTTP methods,
- headers including `Content-Type` and `X-Auth-Token`,
- unauthenticated `OPTIONS` handling.

Add CORS preflight tests that assert `OPTIONS` succeeds without a token and that
`X-Auth-Token` is allowed.

Socket.IO authentication remains auth-payload based:

- If `PYSCF_AUTH_TOKEN` is unset, connection is allowed.
- If set, `auth.token` must match.
- Rejected connections must preserve the current frontend behavior.

Keep the endpoint-specific GPU4PySCF install guard. `POST
/api/system/gpu4pyscf-install` must reject non-loopback clients with the existing
`403` response envelope even when token auth succeeds. Use `Request.client.host`
with the current IPv4/IPv6 loopback helper logic, including IPv4-mapped IPv6
addresses. Add a regression test for the non-loopback `403` envelope.

## SSE Design

`/api/agent/chat` continues returning `text/event-stream` with the current event
format:

```text
data: {"type": "...", "payload": {...}}

```

The existing `_format_sse_event` and chat stream logic can remain synchronous at
first. FastAPI can serve it with `StreamingResponse`. The migration should avoid
rewriting Gemini chat behavior unless required for framework integration.

## Socket.IO Design

Use `python-socketio.AsyncServer(async_mode="asgi")`.

Keep event names and payloads stable:

- `connect`
- `join_calculation`
- `leave_calculation`
- `join_global_updates`
- `leave_global_updates`
- `disconnect`
- emitted `calculation_update`
- emitted `error`

Flask `session` must be replaced by explicit per-SID state in a small helper
owned by `src/python/websocket/handlers.py`. Room operations should use
`await sio.enter_room(sid, room)`, `await sio.leave_room(sid, room)`, and
`await sio.emit(...)`.

File watcher callbacks are synchronous today. When they need to emit from a
non-async callback, use an event-loop-safe bridge rather than blocking in the
watcher thread. Keep the bridge small and covered by tests.

The bridge must be a single shared implementation used by websocket handlers and
`NotificationService`. It should:

1. Capture the running event loop during FastAPI lifespan startup.
2. Schedule async emits from synchronous watcher/process callbacks with
   `asyncio.run_coroutine_threadsafe`.
3. Log failures from the returned future.
4. Clear the stored loop during shutdown.

Do not make notification delivery depend on an opportunistic "loop is available"
check that silently drops calculation updates during normal ASGI runtime.

## Notification Service

`NotificationService` should bind to the Socket.IO async server and support
calculation-update notifications from synchronous process-manager callbacks.

The service should expose a simple synchronous method to callers:

```python
def send_calculation_update(
    self,
    calculation_id: str,
    status: str,
    error_message: str | None = None,
) -> None:
    ...
```

Internally it should submit the async emit through the shared ASGI event-loop
bridge. During normal app lifespan, the loop must be bound. If notification is
called before startup or after shutdown, it should log and return, matching
current behavior. Tests should cover notification scheduling from a background
thread.

## Configuration

`ServerConfig` remains the source for `config/server-config.json`.

Replace `configure_flask_app` with
`configure_fastapi_app(app, config, server_port) -> None`. This function writes
the same logical settings to `app.state` instead of Flask `app.config`.

Existing code that reads Flask `current_app.config` must be decoupled. The known
callers are:

- `src/python/api/health.py`
- `src/python/quantum_calc/config_manager.py`

`quantum_calc/config_manager.py` should read the centralized config directly, or
through a tiny framework-neutral provider. It should no longer import Flask.

## Electron Startup

`src/main/python-server.ts` should still spawn the backend from `src/python`, but
Gunicorn arguments change from WSGI sync worker to ASGI worker:

- target remains `app:app`
- worker class becomes `uvicorn.workers.UvicornWorker`
- keep host, port, timeout, access log, log level, and preload behavior from
  `config/server-config.json` where meaningful.
- keep `workers` at `1` unless an external Socket.IO client manager/message queue
  is added, because in-memory Socket.IO rooms and file-watcher state are not
  safe across multiple workers.

Rename user-facing diagnostics from Flask-specific wording to Python/FastAPI
backend wording. The frontend can continue using `window.flaskPort` for now only
if changing preload/global names would cause unnecessary blast radius; however,
new backend code and diagnostics should avoid adding more Flask naming.

Update all local/package smoke scripts that currently assume Flask/Gunicorn sync
startup:

- `package.json` `test:gunicorn-local`
- `scripts/test-python-standalone.js`
- `scripts/validate-build-completeness.py`
- `scripts/build-python-linux.sh`

These checks should import FastAPI/uvicorn, verify the Uvicorn worker is
available, and start `app:app` through Gunicorn with
`uvicorn.workers.UvicornWorker`.

## Dependencies

Update `.github/environment.yml`:

- Add `fastapi`.
- Add `uvicorn`.
- Keep `gunicorn`.
- Declare `python-socketio` directly.
- Remove Flask runtime dependencies after all imports/tests are migrated:
  `flask`, `flask-cors`, `flask-pydantic`, `flask-socketio`, `flask-sock`,
  `pytest-flask`.

Update scripts and verification messages that currently check Flask imports so
they check FastAPI/uvicorn instead.

## Testing Strategy

Use test-driven migration for behavior changes:

1. Add or convert tests to define desired FastAPI/ASGI behavior.
2. Run the targeted test and confirm it fails for the expected reason.
3. Implement the minimal migration for that behavior.
4. Run the targeted test again.
5. Run the nearby regression set.

Test fixture migration:

- Replace Flask `app.test_client()` with FastAPI `TestClient(fastapi_app)`.
- Replace `app.app_context()` usage with explicit fixture setup.
- Replace Flask-SocketIO test client usage with at least one real ASGI Socket.IO
  smoke test using `socketio.AsyncClient`.
- Update OpenAPI contract tests to parse FastAPI router definitions instead of
  Flask `@blueprint.route` decorators.

The Socket.IO smoke test must cover:

- connection rejected with missing/wrong token when `PYSCF_AUTH_TOKEN` is set,
- connection accepted with the correct token,
- `join_global_updates`,
- `join_calculation`,
- receipt of a `calculation_update` event.

Handler-level tests may still be added for edge cases, but they are not a
substitute for the ASGI protocol smoke test.

Minimum verification before completion:

- `npm run typecheck`
- `npm run codegen`
- `~/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests/integration/test_api_endpoints -v`
- `~/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests/integration/test_auth_security.py src/python/tests/integration/test_auth_production.py -v`
- `~/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests/integration/test_websocket_handlers.py -v`
- `npm run test:gunicorn-local`
- `npm run test:python-standalone`
- `npm run verify-env`
- `npm run verify-build-env`

## Files Expected To Change

Backend:

- `src/python/app.py`
- `src/python/config.py`
- `src/python/api/__init__.py`
- `src/python/api/agent.py`
- `src/python/api/chat_history.py`
- `src/python/api/health.py`
- `src/python/api/pubchem.py`
- `src/python/api/quantum.py`
- `src/python/api/settings.py`
- `src/python/api/smiles.py`
- `src/python/api/swagger_ui.py`
- `src/python/api/system.py`
- `src/python/websocket/handlers.py`
- `src/python/services/notification_service.py`
- `src/python/quantum_calc/config_manager.py`
- `src/python/tests/conftest.py`
- affected integration tests under `src/python/tests/integration/`

Electron/scripts/config:

- `src/main/python-server.ts`
- `src/main/config.ts`
- `config/server-config.json`
- `.github/environment.yml`
- `scripts/verify-environment.py`
- `scripts/validate-build-completeness.py`
- `scripts/build-python-linux.sh`
- `scripts/test-python-standalone.js`
- `package.json`

Documentation/comments may also need wording updates where they explicitly
describe the backend as Flask.

## Risks And Mitigations

Socket.IO migration is the highest-risk area. Mitigate by keeping event names and
payloads unchanged, migrating handler tests early, and preserving room semantics
before touching calculation workflow code.

Electron startup is high risk because packaged startup currently assumes
Gunicorn WSGI. Mitigate with a focused local Gunicorn/Uvicorn worker smoke test
before broad cleanup.

Configuration coupling is medium risk because `current_app.config` leaks into
non-API modules. Mitigate by creating one framework-neutral config provider and
removing Flask imports from service/quantum code.

SSE is lower risk because the current formatter is framework-independent.
Mitigate by testing the raw `data: ...\n\n` stream output.

Dependency cleanup is medium risk. Remove Flask packages only after `rg "flask|Flask|flask_socketio|flask_pydantic"` shows no production imports.

## Out Of Scope

- Changing frontend API paths or response envelopes.
- Replacing Socket.IO with raw WebSocket.
- Changing quantum calculation behavior.
- Changing AI chat provider or prompt behavior.
- Making FastAPI's generated OpenAPI the contract source.
- Renaming `window.flaskPort` everywhere unless it becomes necessary for clarity
  after backend migration is stable.

## Success Criteria

- No production Flask imports remain.
- Backend starts as ASGI through Electron's Python process manager.
- `/health` returns the same status payload shape.
- Existing API endpoints keep their paths, methods, query parameters, and response
  envelopes.
- `/api/agent/chat` streams the same SSE events.
- Socket.IO clients can connect, join/leave calculation rooms, join/leave global
  updates, and receive `calculation_update`.
- Authentication behavior matches current development, test, and production
  expectations.
- Codegen, TypeScript typecheck, targeted backend tests, and environment
  verification pass.
