# FastAPI Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Python backend runtime from Flask/Flask-SocketIO to FastAPI/ASGI while preserving the existing HTTP paths, response envelopes, SSE wire format, Socket.IO event names, authentication behavior, and Electron startup flow.

**Architecture:** `src/python/app.py` becomes the ASGI entry point and exposes `fastapi_app`, `sio`, and `app`. HTTP routes move to FastAPI `APIRouter` modules, SSE uses `StreamingResponse`, Socket.IO uses `python-socketio.AsyncServer(async_mode="asgi")`, and Gunicorn starts `app:app` with `uvicorn.workers.UvicornWorker`.

**Tech Stack:** Python 3.12, FastAPI, Starlette `TestClient`, Pydantic v2, python-socketio ASGI, Uvicorn worker, Gunicorn, Electron main process TypeScript, pytest.

---

## References

- Design spec: `docs/superpowers/specs/2026-05-21-fastapi-migration-design.md`
- API contract: `src/api-spec/openapi.yaml`
- Generated Python models: `src/python/generated_models.py`
- Generated TypeScript API types: `src/web/types/generated-api.ts`
- Server config: `config/server-config.json`

## File Structure

- `src/python/app.py`: ASGI app factory, thread-control invariant, auth/CORS middleware registration, exception handlers, Socket.IO ASGI composition, lifespan startup/shutdown.
- `src/python/config.py`: framework-neutral `configure_fastapi_app(app, config, server_port)` that writes settings to `app.state`.
- `src/python/api/__init__.py`: FastAPI router registration.
- `src/python/api/health.py`: `/health` router reading `request.app.state`.
- `src/python/api/swagger_ui.py`: development-only `/api-docs` and `/api-docs/spec.json` routers.
- `src/python/api/settings.py`, `src/python/api/pubchem.py`, `src/python/api/smiles.py`, `src/python/api/chat_history.py`: service-backed HTTP routers with preserved response envelopes.
- `src/python/api/quantum.py`: quantum HTTP router preserving raw-body validation order for `/api/quantum/calculate`.
- `src/python/api/system.py`: system HTTP router preserving diagnostics and loopback-only GPU4PySCF install guard.
- `src/python/api/agent.py`: `/api/agent/chat` `StreamingResponse` SSE route.
- `src/python/websocket/event_loop_bridge.py`: shared sync-to-async scheduling bridge for callback-driven emits.
- `src/python/websocket/handlers.py`: async Socket.IO handlers and per-SID room state.
- `src/python/services/notification_service.py`: synchronous notification API that schedules async Socket.IO emits through the bridge.
- `src/python/quantum_calc/config_manager.py`: framework-neutral quantum config access.
- `src/python/tests/conftest.py`: FastAPI `TestClient`, shared app fixtures, Socket.IO ASGI server fixtures.
- `src/python/tests/integration/test_app_startup_fastapi.py`: ASGI startup and docs disabling regressions.
- `src/python/tests/integration/test_validation_fastapi.py`: validation envelope, malformed JSON, query validation, CORS preflight.
- `src/python/tests/integration/test_socketio_asgi.py`: real Socket.IO ASGI smoke tests.
- `src/python/tests/integration/test_websocket_handlers.py`: handler-level behavior converted away from Flask-SocketIO test client.
- `src/python/tests/integration/test_api_endpoints/test_openapi_contract.py`: contract scanner updated from Flask decorators to FastAPI decorators.
- `.github/environment.yml`: FastAPI/Uvicorn/python-socketio dependencies and Flask dependency removal.
- `package.json`: verification scripts and Gunicorn smoke command.
- `scripts/verify-environment.py`: FastAPI/Uvicorn import and app smoke checks.
- `scripts/validate-build-completeness.py`: packaged dependency validation.
- `scripts/build-python-linux.sh`: Linux package dependency check.
- `scripts/test-python-standalone.js`: standalone backend import and Gunicorn/Uvicorn worker smoke checks.
- `src/main/python-server.ts`: Electron-managed Gunicorn arguments and diagnostics.
- `src/main/config.ts`: server config type update for ASGI worker settings.
- `config/server-config.json`: worker class update and worker count invariant.

## Commit Cadence

Commit after each task passes its targeted tests. Use short commits such as:

```bash
git add <changed files>
git commit -m "feat: add FastAPI app shell"
```

Run `git status --short` before each commit and keep unrelated changes out of the commit.

---

### Task 1: Dependencies And Smoke Scripts

**Files:**
- Modify: `.github/environment.yml`
- Modify: `package.json`
- Modify: `scripts/verify-environment.py`
- Modify: `scripts/validate-build-completeness.py`
- Modify: `scripts/build-python-linux.sh`
- Modify: `scripts/test-python-standalone.js`

- [ ] **Step 1: Add a failing dependency smoke check**

Run this before editing dependencies:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python - <<'PY'
from fastapi import FastAPI
from fastapi.testclient import TestClient
from uvicorn.workers import UvicornWorker
import socketio

server = socketio.AsyncServer(async_mode="asgi")
app = FastAPI()
TestClient(app)
print(FastAPI.__name__, TestClient.__name__, UvicornWorker.__name__, server.async_mode)
PY
```

Expected: FAIL with `ModuleNotFoundError: No module named 'fastapi'` or `No module named 'uvicorn'` in the current conda environment.

- [ ] **Step 2: Update `.github/environment.yml`**

Add FastAPI ASGI dependencies first. Keep Flask dependencies in this task because the migration is not implemented yet; remove them in Task 11 after the runtime import search is clean. Keep the exact FastAPI versions aligned with the latest conda resolution used during implementation.

```yaml
  - fastapi=0.136.1
  - uvicorn=0.47.0
  - python-socketio=5.16.1
  - gunicorn=23.0.0
```

Leave these Flask-specific entries in place until Task 11:

```yaml
  - flask=3.1.2
  - flask-cors=6.0.1
  - pip:
    - flask-pydantic==0.13.2
    - flask-socketio==5.5.1
    - flask-sock==0.7.0
    - pytest-flask==1.3.0
```

- [ ] **Step 3: Update `package.json` verification scripts**

Change `test:gunicorn-local` to start the ASGI app through Uvicorn worker:

```json
"test:gunicorn-local": "bash -c 'CONDA_BASE=$(conda info --base 2>/dev/null || echo \"$HOME/miniforge3\") && source \"$CONDA_BASE/etc/profile.d/conda.sh\" && conda activate pyscf-env && cd src/python && echo \"Testing Gunicorn/Uvicorn locally...\" && gunicorn --bind 127.0.0.1:5000 --workers 1 --worker-class uvicorn.workers.UvicornWorker --timeout 0 --log-level info app:app'"
```

Change `test:conda-env` so it checks FastAPI and Uvicorn instead of Flask:

```json
"test:conda-env": "bash -c 'CONDA_BASE=$(conda info --base 2>/dev/null || echo \"$HOME/miniforge3\") && export PATH=\"$CONDA_BASE/envs/pyscf-env/bin:$PATH\" && cd dist/mac-arm64/Pyscf_front.app/Contents/Resources/conda_env/bin && echo \"Testing conda environment...\" && ./python --version && ./python -c \"import gunicorn; print('\"'\"'Gunicorn:'\"'\"', gunicorn.__version__)\" && ./python -c \"import fastapi, uvicorn, socketio; print('\"'\"'FastAPI ASGI stack: Available'\"'\"')\" && echo \"Conda environment test completed.\"'"
```

Extend `verify-build-env` with FastAPI/Uvicorn imports:

```json
"verify-build-env": "bash -c 'CONDA_BASE=$(conda info --base 2>/dev/null || echo \"$HOME/miniforge3\") && source \"$CONDA_BASE/etc/profile.d/conda.sh\" && conda activate pyscf-env && echo \"=== Build Environment Verification ===\" && $(conda info --base)/envs/pyscf-env/bin/python -c \"import conda_pack; print(f\\\"conda-pack: {conda_pack.__version__}\\\")\" && $(conda info --base)/envs/pyscf-env/bin/python -c \"import gunicorn; print(f\\\"Gunicorn: {gunicorn.__version__}\\\")\" && $(conda info --base)/envs/pyscf-env/bin/python -c \"import fastapi, uvicorn, socketio; print(\\\"FastAPI ASGI stack verified\\\")\" && echo \"✓ All build tools verified\"'"
```

- [ ] **Step 4: Update Python packaging verification scripts**

In `scripts/verify-environment.py`, replace Flask dependency keys with:

```python
dist_name_map = {
    'pyscf': 'pyscf',
    'rdkit': 'rdkit',
    'geometric': 'geometric',
    'fastapi': 'fastapi',
    'uvicorn': 'uvicorn',
    'socketio': 'python-socketio',
    'pydantic': 'pydantic',
    'gunicorn': 'gunicorn',
    'requests': 'requests',
}

required_packages = [
    ('pyscf', 'PySCF - 量子化学計算'),
    ('rdkit', 'RDKit - 化学情報学'),
    ('geometric', 'geometric - 分子幾何最適化'),
    ('fastapi', 'FastAPI - ASGI Web フレームワーク'),
    ('uvicorn', 'Uvicorn - ASGI サーバー'),
    ('socketio', 'python-socketio - Socket.IO ASGI対応'),
    ('pydantic', 'Pydantic - データバリデーション'),
    ('gunicorn', 'Gunicorn - 本番プロセスマネージャ'),
    ('requests', 'Requests - HTTP クライアント'),
]
```

Rename the Flask functionality check to `check_fastapi_functionality` and use this exact smoke app:

```python
def check_fastapi_functionality() -> bool:
    """FastAPI の基本機能をテスト"""
    log_info("FastAPI の基本機能をテスト中...")
    try:
        from fastapi import FastAPI
        from fastapi.testclient import TestClient
        from uvicorn.workers import UvicornWorker
        import socketio

        app = FastAPI(docs_url=None, redoc_url=None, openapi_url=None)

        @app.get('/test')
        def test_endpoint():
            return {'message': 'FastAPI test successful'}

        client = TestClient(app)
        response = client.get('/test')
        server = socketio.AsyncServer(async_mode='asgi')

        if response.status_code == 200 and server.async_mode == 'asgi' and UvicornWorker is not None:
            log_success("FastAPI ASGI 機能テスト成功 ✓")
            return True

        log_error(f"FastAPI テスト失敗: ステータスコード {response.status_code}")
        return False
    except Exception as e:
        log_error(f"FastAPI テストに失敗: {e}")
        return False
```

Update `main()` so the test list calls the renamed function:

```python
tests = [
    ("Python バージョン", check_python_version),
    ("必須パッケージ", check_required_packages),
    ("PySCF 機能", check_pyscf_functionality),
    ("RDKit 機能", check_rdkit_functionality),
    ("FastAPI ASGI 機能", check_fastapi_functionality),
    ("conda 環境", check_conda_environment),
    ("プロジェクト構造", check_project_structure),
]
```

In `scripts/validate-build-completeness.py`, replace the packaged import list:

```python
required_modules = [
    "fastapi",
    "uvicorn",
    "socketio",
    "gunicorn",
    "pydantic",
]
```

In `scripts/build-python-linux.sh`, replace the Flask import smoke with:

```bash
python -c "import fastapi; print(f\"FastAPI: {fastapi.__version__}\")"
python -c "import uvicorn; print(f\"Uvicorn: {uvicorn.__version__}\")"
python -c "import socketio; print('python-socketio: available')"
```

In `scripts/test-python-standalone.js`, replace the Flask import test with:

```javascript
console.log('\n=== Test 4: FastAPI ASGI Import Test ===');
await testPythonCommand([
  pythonExecutablePath,
  '-c',
  'import fastapi, uvicorn, socketio; from uvicorn.workers import UvicornWorker; print("FastAPI ASGI import successful")'
]);
```

Add `--worker-class uvicorn.workers.UvicornWorker` to the standalone Gunicorn args:

```javascript
const gunicornArgs = [
  '-m', 'gunicorn',
  '--bind', '127.0.0.1:0',
  '--workers', '1',
  '--worker-class', 'uvicorn.workers.UvicornWorker',
  '--timeout', '5',
  '--log-level', 'warning',
  'app:app',
];
```

- [ ] **Step 5: Install or update the conda environment**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
conda env update -n pyscf-env -f .github/environment.yml --prune
```

Expected: conda resolves FastAPI, Uvicorn, python-socketio, and existing PySCF/RDKit packages without conflicts.

- [ ] **Step 6: Verify dependency smoke passes**

Run the smoke command from Step 1 again.

Expected output includes:

```text
FastAPI TestClient UvicornWorker asgi
```

- [ ] **Step 7: Commit**

```bash
git add .github/environment.yml package.json scripts/verify-environment.py scripts/validate-build-completeness.py scripts/build-python-linux.sh scripts/test-python-standalone.js
git commit -m "build: add FastAPI ASGI dependencies"
```

---

### Task 2: FastAPI App Shell, Config State, And Test Fixtures

**Files:**
- Modify: `src/python/app.py`
- Modify: `src/python/config.py`
- Modify: `src/python/api/__init__.py`
- Modify: `src/python/api/health.py`
- Modify: `src/python/api/swagger_ui.py`
- Create: `src/python/websocket/event_loop_bridge.py`
- Modify: `src/python/tests/conftest.py`
- Modify: `src/python/tests/test_fixtures.py`
- Create: `src/python/tests/integration/test_app_startup_fastapi.py`

- [ ] **Step 1: Write failing FastAPI startup tests**

Create `src/python/tests/integration/test_app_startup_fastapi.py`:

```python
import importlib
import os

from fastapi import FastAPI


def test_thread_control_vars_are_set_before_app_import(monkeypatch):
    thread_vars = [
        'OMP_NUM_THREADS',
        'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS',
        'BLIS_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS',
        'NUMEXPR_NUM_THREADS',
    ]
    for name in thread_vars:
        monkeypatch.delenv(name, raising=False)

    module = importlib.reload(importlib.import_module('app'))

    for name in thread_vars:
        assert os.environ[name] == '1'
    assert hasattr(module, 'fastapi_app')
    assert hasattr(module, 'sio')
    assert hasattr(module, 'app')


def test_create_fastapi_app_disables_default_docs():
    from app import create_fastapi_app

    app = create_fastapi_app(server_port=5000, test_config={'TESTING': True})

    assert isinstance(app, FastAPI)
    routes = {route.path for route in app.routes}
    assert '/docs' not in routes
    assert '/redoc' not in routes
    assert '/openapi.json' not in routes
    assert '/api-docs' in routes
    assert '/api-docs/spec.json' in routes


def test_config_is_stored_on_app_state():
    from app import create_fastapi_app

    app = create_fastapi_app(
        server_port=5123,
        test_config={
            'TESTING': True,
            'CALCULATIONS_DIR': '/tmp/pyscf-fastapi-test',
            'WEBSOCKET_WATCHER_ENABLED': False,
        },
    )

    assert app.state.TESTING is True
    assert app.state.SERVER_PORT == 5123
    assert app.state.CALCULATIONS_DIR == '/tmp/pyscf-fastapi-test'
    assert app.state.WEBSOCKET_WATCHER_ENABLED is False
```

- [ ] **Step 2: Run startup tests to verify they fail**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_app_startup_fastapi.py -v
```

Expected: FAIL because `create_fastapi_app` and `fastapi_app` do not exist and the module still creates a Flask app.

- [ ] **Step 3: Add framework-neutral FastAPI config**

In `src/python/config.py`, keep `ServerConfig`, `get_server_config`, and `determine_server_port`. Replace `configure_flask_app` with:

```python
def configure_fastapi_app(app, config: ServerConfig, server_port: int) -> None:
    """
    Configure FastAPI application state from ServerConfig.

    FastAPI has no app.config dict, so application settings are stored on
    app.state with the same logical keys the Flask runtime used.
    """
    app.state.SERVER_CONFIG = config.to_dict()
    app.state.SERVER_CONFIG_OBJECT = config
    app.state.SERVER_HOST = config.get_server_host()
    app.state.SERVER_PORT = server_port
    app.state.DEBUG = config.get('server.debug', False)
    app.state.TESTING = False

    app.state.GUNICORN = config.get('gunicorn', {})
    app.state.SOCKETIO = config.get('socketio', {})
    app.state.DEVELOPMENT = config.get('development', {})
    app.state.PRODUCTION = config.get('production', {})
    app.state.QUANTUM_CALCULATIONS = config.get('quantum_calculations', {})
    app.state.QUANTUM_CALCULATION_DEFAULTS = config.get('quantum_calculation_defaults', {})
    app.state.APP_INFO = config.get('app_info', {})
    app.state.APP_VERSION = config.get('app_info.version', 'unknown')
    app.state.EXTERNAL_SERVICES = config.get('external_services', {})
    app.state.LOGGING = config.get('logging', {})
    app.state.AI_AGENT = config.get('ai_agent', {})
    app.state.VERSION = config.get('app_info.version', 'unknown')

    logger.info(f"FastAPI app configured with server port: {server_port}")
```

Preserve `configure_flask_app` only until Task 11 if imports still need it. Once `rg "configure_flask_app"` returns only its definition, delete it.

- [ ] **Step 4: Add the minimal event-loop bridge used by app lifespan**

Create `src/python/websocket/event_loop_bridge.py`:

```python
import asyncio
import logging

logger = logging.getLogger(__name__)
_loop: asyncio.AbstractEventLoop | None = None


def bind_event_loop(loop: asyncio.AbstractEventLoop) -> None:
    global _loop
    _loop = loop
    logger.info("ASGI event loop bound for background notifications")


def clear_event_loop() -> None:
    global _loop
    _loop = None
    logger.info("ASGI event loop cleared for background notifications")
```

- [ ] **Step 5: Replace `src/python/app.py` with ASGI app shell**

Keep the thread-control block as the first executable block. Then build the ASGI app with this structure:

```python
import os

_THREAD_CONTROL_VARS = [
    'OMP_NUM_THREADS',
    'MKL_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    'BLIS_NUM_THREADS',
    'VECLIB_MAXIMUM_THREADS',
    'NUMEXPR_NUM_THREADS',
]

for _var in _THREAD_CONTROL_VARS:
    if _var not in os.environ:
        os.environ[_var] = '1'

import asyncio
import json
import logging
import sys
from contextlib import asynccontextmanager
from typing import Any

import socketio
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from pydantic import ValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException

from api import register_routers
from config import (
    ConfigurationError,
    configure_fastapi_app,
    determine_server_port,
    get_server_config,
)
from quantum_calc import shutdown_process_manager, shutdown_websocket_watcher
from services.exceptions import ServiceError
from websocket import register_websocket_handlers
from websocket.event_loop_bridge import bind_event_loop, clear_event_loop
```

Add the app factory and top-level exports:

```python
logger = logging.getLogger(__name__)


def _is_development_api_docs_path(path: str) -> bool:
    return path == '/api-docs' or path.startswith('/api-docs/')


def _get_state(app: FastAPI, name: str, default: Any = None) -> Any:
    return getattr(app.state, name, default)


@asynccontextmanager
async def lifespan(fastapi_app: FastAPI):
    bind_event_loop(asyncio.get_running_loop())
    try:
        yield
    finally:
        clear_event_loop()
        shutdown_websocket_watcher()
        shutdown_process_manager()


def create_socketio_server(socketio_config: dict[str, Any] | None = None) -> socketio.AsyncServer:
    socketio_config = socketio_config or {}
    cors_allowed_origins = socketio_config.get('cors_allowed_origins')
    if cors_allowed_origins not in (None, '*'):
        cors_allowed_origins = [
            origin
            for origin in cors_allowed_origins
            if not origin.endswith(':*') and origin not in {'ws://127.0.0.1:*', 'ws://localhost:*'}
        ]
        cors_allowed_origins.extend([
            'http://127.0.0.1:3000',
            'http://localhost:3000',
            'http://127.0.0.1:5173',
            'http://localhost:5173',
        ])
    return socketio.AsyncServer(
        async_mode='asgi',
        cors_allowed_origins=cors_allowed_origins or ['file://', 'http://127.0.0.1:3000', 'http://localhost:3000'],
        ping_timeout=socketio_config.get('ping_timeout', 60),
        ping_interval=socketio_config.get('ping_interval', 25),
        logger=socketio_config.get('logger', True),
        engineio_logger=socketio_config.get('engineio_logger', False),
    )


def create_fastapi_app(server_port: int | None = None, test_config: dict[str, Any] | None = None) -> FastAPI:
    server_config = get_server_config()
    if server_port is None:
        port_env = os.getenv('PYSCF_SERVER_PORT')
        server_port = determine_server_port(server_config, port_env=port_env)

    fastapi_app = FastAPI(docs_url=None, redoc_url=None, openapi_url=None, lifespan=lifespan)
    configure_fastapi_app(fastapi_app, server_config, server_port)

    if test_config:
        for key, value in test_config.items():
            setattr(fastapi_app.state, key, value)

    register_routers(fastapi_app)
    return fastapi_app


def compose_asgi_app(fastapi_instance: FastAPI, socketio_instance: socketio.AsyncServer):
    register_websocket_handlers(socketio_instance)

    from services.notification_service import bind_notification_service, get_notification_service
    from quantum_calc import initialize_process_manager_with_callback

    bind_notification_service(socketio_instance)
    try:
        notification_service = get_notification_service()
        initialize_process_manager_with_callback(
            notification_callback=notification_service.send_calculation_update
        )
    except Exception as exc:
        logger.error("Failed to initialize process manager with callback: %s", exc)

    fastapi_instance.state.socketio = socketio_instance
    return socketio.ASGIApp(socketio_instance, fastapi_instance)


def create_app(server_port: int | None = None, test_config: dict[str, Any] | None = None):
    fastapi_instance = create_fastapi_app(server_port=server_port, test_config=test_config)
    socketio_config = getattr(fastapi_instance.state, 'SOCKETIO', {})
    socketio_instance = create_socketio_server(socketio_config)
    return compose_asgi_app(fastapi_instance, socketio_instance)


fastapi_app = create_fastapi_app()
sio = create_socketio_server(getattr(fastapi_app.state, 'SOCKETIO', {}))
app = compose_asgi_app(fastapi_app, sio)
```

Keep the `if __name__ == '__main__'` path, but run Uvicorn:

```python
if __name__ == '__main__':
    import uvicorn

    server_config = get_server_config()
    port_env = os.getenv('PYSCF_SERVER_PORT')
    actual_port = determine_server_port(server_config, port_env=port_env)
    host = server_config.get('server.host', '127.0.0.1')
    uvicorn.run('app:app', host=host, port=actual_port, reload=False, log_level='info')
```

- [ ] **Step 6: Convert `health.py` and `swagger_ui.py` for startup coverage**

Replace the Flask health blueprint with:

```python
from fastapi import APIRouter, Request

router = APIRouter()


@router.get('/health')
def health_check(request: Request):
    version = getattr(request.app.state, 'VERSION', '1.0.0')
    return {
        'status': 'ok',
        'service': 'pyscf-front-api',
        'version': version,
    }
```

Replace `swagger_ui.py` with:

```python
from pathlib import Path

import yaml
from fastapi import APIRouter
from fastapi.responses import HTMLResponse, JSONResponse

router = APIRouter()


@router.get('/api-docs', response_class=HTMLResponse)
def api_docs():
    return HTMLResponse(
        """
        <!doctype html>
        <html>
          <head>
            <title>PySCF Front API Docs</title>
            <link rel="stylesheet" href="https://unpkg.com/swagger-ui-dist/swagger-ui.css">
          </head>
          <body>
            <div id="swagger-ui"></div>
            <script src="https://unpkg.com/swagger-ui-dist/swagger-ui-bundle.js"></script>
            <script>
              SwaggerUIBundle({ url: '/api-docs/spec.json', dom_id: '#swagger-ui' });
            </script>
          </body>
        </html>
        """
    )


@router.get('/api-docs/spec.json')
def api_docs_spec():
    spec_path = Path(__file__).resolve().parents[2] / 'api-spec' / 'openapi.yaml'
    with spec_path.open('r', encoding='utf-8') as handle:
        return JSONResponse(yaml.safe_load(handle))
```

- [ ] **Step 7: Convert `src/python/api/__init__.py` to minimal router registration**

Register only routers converted in this task. Later tasks extend this file as each router is migrated.

```python
import logging
import os

from fastapi import FastAPI

from .health import router as health_router
from .swagger_ui import router as swagger_router

logger = logging.getLogger(__name__)


def register_routers(app: FastAPI) -> None:
    app.include_router(health_router)
    is_packaged = os.getenv('PYSCF_RESOURCES_PATH') is not None
    if not is_packaged:
        app.include_router(swagger_router)
        logger.info("Swagger UI registered at /api-docs/ (development mode)")
    else:
        logger.info("Swagger UI not registered (packaged mode)")
    logger.info("Registered startup FastAPI routers")
```

- [ ] **Step 8: Convert test fixtures to FastAPI TestClient**

In `src/python/tests/conftest.py`, replace Flask imports and clients with:

```python
from fastapi.testclient import TestClient

from app import create_fastapi_app
```

Use this `app` fixture body:

```python
@pytest.fixture(scope='function')
def app():
    temp_dir = tempfile.mkdtemp(prefix='pyscf_test_')
    test_config = {
        'TESTING': True,
        'CALCULATIONS_DIR': temp_dir,
        'WEBSOCKET_WATCHER_ENABLED': False,
        'SOCKETIO': {
            'cors_allowed_origins': ['http://127.0.0.1:*', 'http://localhost:*', 'file://'],
            'logger': False,
            'engineio_logger': False,
        },
    }

    import services as services_module
    import quantum_calc.settings_manager as settings_manager_module
    from quantum_calc.settings_manager import SettingsManager
    from quantum_calc.process_manager import shutdown_process_manager

    test_settings_manager = SettingsManager(settings_file=os.path.join(temp_dir, "app-settings.json"))
    test_settings = test_settings_manager.get_default_settings().model_copy(
        update={'calculations_directory': temp_dir}
    )
    test_settings_manager.save_settings(test_settings)

    with (
        mock.patch.dict(os.environ, {'PYSCF_ENV': 'development'}),
        mock.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor),
        mock.patch.object(settings_manager_module, "_settings_manager", test_settings_manager),
        mock.patch.multiple(
            services_module,
            _quantum_service=None,
            _pubchem_service=None,
            _smiles_service=None,
            _settings_service=None,
            _system_service=None,
        ),
    ):
        shutdown_process_manager()
        _app = create_fastapi_app(server_port=5000, test_config=test_config)
        yield _app
        shutdown_process_manager()

    import shutil
    shutil.rmtree(temp_dir, ignore_errors=True)
```

Use this HTTP client fixture:

```python
@pytest.fixture(scope='function')
def client(app):
    with TestClient(app) as test_client:
        yield test_client
```

Remove the Flask `app.app_context()` and Flask-SocketIO `socketio_client` fixture in this task. A real Socket.IO fixture is added in Task 9.

- [ ] **Step 9: Update fixture self-tests**

In `src/python/tests/test_fixtures.py`, convert Flask config assertions to FastAPI state assertions and remove the deleted Socket.IO fixture test:

```python
def test_app_fixture(app):
    assert app.state.TESTING is True
    assert hasattr(app.state, 'CALCULATIONS_DIR')


def test_app_fixture_isolates_current_settings_directory(app):
    from quantum_calc import get_current_settings

    settings = get_current_settings()

    assert settings.calculations_directory == app.state.CALCULATIONS_DIR
```

Delete `test_socketio_client_fixture`; the real ASGI Socket.IO fixture is tested in `src/python/tests/integration/test_socketio_asgi.py`.

- [ ] **Step 10: Run startup and fixture tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_app_startup_fastapi.py -v
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/test_fixtures.py -v
```

Expected: PASS.

- [ ] **Step 11: Commit**

```bash
git add src/python/app.py src/python/config.py src/python/api/__init__.py src/python/api/health.py src/python/api/swagger_ui.py src/python/websocket/event_loop_bridge.py src/python/tests/conftest.py src/python/tests/test_fixtures.py src/python/tests/integration/test_app_startup_fastapi.py
git commit -m "feat: add FastAPI ASGI app shell"
```

---

### Task 3: HTTP Auth, CORS, And Error Envelopes

**Files:**
- Modify: `src/python/app.py`
- Create: `src/python/tests/integration/test_validation_fastapi.py`
- Modify: `src/python/tests/integration/test_auth_security.py`
- Modify: `src/python/tests/integration/test_auth_production.py`

- [ ] **Step 1: Write failing validation and CORS tests**

Create `src/python/tests/integration/test_validation_fastapi.py`:

```python
from fastapi import Query
from fastapi.testclient import TestClient
from pydantic import BaseModel


class SampleBody(BaseModel):
    required_name: str


def test_validation_error_returns_400_envelope(app):
    @app.post('/_test/validation-body')
    def validation_body(body: SampleBody):
        return {'ok': True}

    with TestClient(app) as client:
        response = client.post('/_test/validation-body', json={})

    assert response.status_code == 400
    assert response.json()['success'] is False
    assert response.json()['error'].startswith('Validation failed:')


def test_malformed_json_returns_400_envelope(app):
    @app.post('/_test/validation-json')
    def validation_json(body: SampleBody):
        return {'ok': True}

    with TestClient(app) as client:
        response = client.post(
            '/_test/validation-json',
            content=b'{"required_name":',
            headers={'Content-Type': 'application/json'},
        )

    assert response.status_code == 400
    assert response.json()['success'] is False
    assert response.json()['error'].startswith('Validation failed:')


def test_invalid_query_parameter_returns_400_envelope(app):
    @app.get('/_test/validation-query')
    def validation_query(limit: int = Query(10)):
        return {'limit': limit}

    with TestClient(app) as client:
        response = client.get('/_test/validation-query?limit=not-an-int')

    assert response.status_code == 400
    assert response.json()['success'] is False
    assert response.json()['error'].startswith('Validation failed:')


def test_options_preflight_succeeds_without_auth(monkeypatch, app):
    monkeypatch.setenv('PYSCF_AUTH_TOKEN', 'secret-token')

    with TestClient(app) as client:
        response = client.options(
            '/health',
            headers={
                'Origin': 'http://127.0.0.1:3000',
                'Access-Control-Request-Method': 'PATCH',
                'Access-Control-Request-Headers': 'X-Auth-Token, Content-Type',
            },
        )

    assert response.status_code in {200, 204}
    assert 'X-Auth-Token' in response.headers['access-control-allow-headers']


def test_auth_error_includes_cors_headers(monkeypatch, app):
    monkeypatch.setenv('PYSCF_AUTH_TOKEN', 'secret-token')

    with TestClient(app) as client:
        response = client.get(
            '/health',
            headers={
                'Origin': 'http://127.0.0.1:3000',
                'X-Auth-Token': 'wrong-token',
            },
        )

    assert response.status_code == 401
    assert response.json()['success'] is False
    assert response.headers['access-control-allow-origin'] == 'http://127.0.0.1:3000'


def test_unhandled_exception_returns_500_envelope(app):
    @app.get('/_test/unhandled-error')
    def unhandled_error():
        raise RuntimeError('boom')

    with TestClient(app, raise_server_exceptions=False) as client:
        response = client.get('/_test/unhandled-error')

    assert response.status_code == 500
    assert response.json() == {
        'success': False,
        'error': 'An internal server error occurred.',
    }
```

- [ ] **Step 2: Run validation tests to verify they fail**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_validation_fastapi.py -v
```

Expected: FAIL because exception handlers, auth middleware, and CORS middleware are not registered yet.

- [ ] **Step 3: Add FastAPI auth middleware**

In `src/python/app.py`, add:

```python
from fastapi import FastAPI, Request
```

Then add:

```python
def register_auth_middleware(fastapi_app: FastAPI) -> None:
    @fastapi_app.middleware('http')
    async def verify_auth_token(request: Request, call_next):
        if _is_development_api_docs_path(request.url.path) and os.getenv('PYSCF_ENV') == 'development':
            return await call_next(request)

        if request.method == 'OPTIONS':
            return await call_next(request)

        auth_token = os.getenv('PYSCF_AUTH_TOKEN')
        if auth_token:
            client_token = request.headers.get('X-Auth-Token')
            if client_token != auth_token:
                client_host = request.client.host if request.client else 'unknown'
                logger.warning("Unauthorized access attempt from %s", client_host)
                return JSONResponse({'success': False, 'error': 'Unauthorized'}, status_code=401)
        else:
            is_testing = bool(getattr(request.app.state, 'TESTING', False))
            env = os.getenv('PYSCF_ENV')
            if is_testing and env != 'production':
                return await call_next(request)
            if env not in {'development', 'test'}:
                logger.warning("Unauthorized access attempt: Missing authentication token in production mode")
                return JSONResponse(
                    {'success': False, 'error': 'Unauthorized: Missing authentication token'},
                    status_code=401,
                )
            logger.warning("Running without authentication token in debug/development mode!")

        return await call_next(request)
```

Call `register_auth_middleware(fastapi_app)` inside `create_fastapi_app` after `configure_fastapi_app` and before `register_routers`.

- [ ] **Step 4: Add HTTP CORS middleware**

In `src/python/app.py`, add:

```python
def register_cors_middleware(fastapi_app: FastAPI) -> None:
    fastapi_app.add_middleware(
        CORSMiddleware,
        allow_origins=[
            'http://127.0.0.1',
            'http://localhost',
            'file://',
            'null',
        ],
        allow_origin_regex=r'^(https?://(127\.0\.0\.1|localhost)(:\d+)?)$',
        allow_credentials=True,
        allow_methods=['GET', 'POST', 'PUT', 'PATCH', 'DELETE', 'OPTIONS'],
        allow_headers=['Content-Type', 'X-Auth-Token'],
    )
```

Call `register_auth_middleware(fastapi_app)` first and `register_cors_middleware(fastapi_app)` after it. Starlette wraps the last added middleware outermost, so CORS must be added after auth to attach CORS headers to auth-error responses.

- [ ] **Step 5: Add error envelope handlers**

In `src/python/app.py`, add:

```python
def _format_validation_errors(errors: list[dict[str, Any]]) -> str:
    messages: list[str] = []
    for err in errors:
        loc = '.'.join(str(part) for part in err.get('loc', []))
        msg = err.get('msg', 'Invalid value')
        messages.append(f"{loc}: {msg}" if loc else msg)
    return 'Validation failed: ' + '; '.join(messages)


def register_exception_handlers(fastapi_app: FastAPI) -> None:
    @fastapi_app.exception_handler(json.JSONDecodeError)
    async def json_decode_error_handler(request: Request, error: json.JSONDecodeError):
        message = f"Validation failed: malformed JSON body: {error.msg}"
        logger.warning("Malformed JSON on %s: %s", request.url.path, message)
        return JSONResponse({'success': False, 'error': message}, status_code=400)

    @fastapi_app.exception_handler(RequestValidationError)
    async def request_validation_error_handler(request: Request, error: RequestValidationError):
        message = _format_validation_errors(error.errors())
        logger.warning("Validation error on %s: %s", request.url.path, message)
        return JSONResponse({'success': False, 'error': message}, status_code=400)

    @fastapi_app.exception_handler(ValidationError)
    async def pydantic_validation_error_handler(request: Request, error: ValidationError):
        message = _format_validation_errors(error.errors())
        logger.warning("Pydantic validation error on %s: %s", request.url.path, message)
        return JSONResponse({'success': False, 'error': message}, status_code=400)

    @fastapi_app.exception_handler(StarletteHTTPException)
    async def http_exception_handler(request: Request, error: StarletteHTTPException):
        if error.status_code == 404:
            return JSONResponse({'success': False, 'error': 'Not Found'}, status_code=404)
        if error.status_code == 400:
            return JSONResponse(
                {'success': False, 'error': f'Validation failed: {error.detail}'},
                status_code=400,
            )
        return JSONResponse(
            {'success': False, 'error': str(error.detail)},
            status_code=error.status_code,
        )

    @fastapi_app.exception_handler(ServiceError)
    async def service_error_handler(request: Request, error: ServiceError):
        return JSONResponse(
            {'success': False, 'error': error.message},
            status_code=error.status_code,
        )

    @fastapi_app.exception_handler(Exception)
    async def unhandled_exception_handler(request: Request, error: Exception):
        logger.error("Unhandled exception on %s: %s", request.url.path, error, exc_info=True)
        return JSONResponse(
            {'success': False, 'error': 'An internal server error occurred.'},
            status_code=500,
        )
```

Call `register_exception_handlers(fastapi_app)` inside `create_fastapi_app` before `register_routers`.

- [ ] **Step 6: Update auth tests to use TestClient**

In `src/python/tests/integration/test_auth_security.py` and `src/python/tests/integration/test_auth_production.py`, replace Flask-specific app construction with:

```python
from fastapi.testclient import TestClient

from app import create_fastapi_app


def make_auth_client(monkeypatch, token='test-token', env='development'):
    monkeypatch.setenv('PYSCF_AUTH_TOKEN', token)
    monkeypatch.setenv('PYSCF_ENV', env)
    app = create_fastapi_app(server_port=5000, test_config={'TESTING': env != 'production'})
    return TestClient(app)
```

Keep existing assertions for:

```python
assert response.status_code == 401
assert response.json() == {'success': False, 'error': 'Unauthorized'}
```

For missing token in production, assert:

```python
assert response.json() == {
    'success': False,
    'error': 'Unauthorized: Missing authentication token',
}
```

- [ ] **Step 7: Run auth and validation tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_validation_fastapi.py tests/integration/test_auth_security.py tests/integration/test_auth_production.py -v
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add src/python/app.py src/python/tests/integration/test_validation_fastapi.py src/python/tests/integration/test_auth_security.py src/python/tests/integration/test_auth_production.py
git commit -m "feat: add FastAPI auth and error handling"
```

---

### Task 4: Convert Simple HTTP Routers

**Files:**
- Modify: `src/python/api/__init__.py`
- Modify: `src/python/api/settings.py`
- Modify: `src/python/api/pubchem.py`
- Modify: `src/python/api/smiles.py`
- Modify: `src/python/api/chat_history.py`
- Create: `src/python/tests/integration/test_api_endpoints/test_settings_api.py`
- Modify: tests under `src/python/tests/integration/test_api_endpoints/`

- [ ] **Step 1: Run simple endpoint tests to capture failing baseline**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_health_api.py tests/integration/test_api_endpoints/test_pubchem_api.py tests/integration/test_api_endpoints/test_smiles_api.py -v
```

Expected: FAIL on Flask blueprint/test client assumptions.

- [ ] **Step 2: Add settings endpoint coverage**

Create `src/python/tests/integration/test_api_endpoints/test_settings_api.py`:

```python
class TestSettingsAPI:
    def _settings_payload(self, calculations_directory='/tmp/pyscf-test'):
        return {
            'max_parallel_instances': 4,
            'max_cpu_utilization_percent': 95.0,
            'max_memory_utilization_percent': 95.0,
            'gpu_acceleration_enabled': False,
            'system_total_cores': 8,
            'system_total_memory_mb': 16384,
            'calculations_directory': calculations_directory,
            'timezone': 'UTC',
            'gemini_api_key': None,
            'research_email': None,
        }

    def test_get_settings_success(self, client, mocker):
        mock_settings = self._settings_payload()
        mock_service = mocker.patch('api.settings.get_settings_service')
        mock_service.return_value.get_settings.return_value = mock_settings

        response = client.get('/api/settings')

        assert response.status_code == 200
        assert response.json()['success'] is True
        assert response.json()['data']['settings'] == mock_settings

    def test_update_settings_success(self, client, mocker):
        updated = self._settings_payload('/tmp/pyscf-updated')
        mock_service = mocker.patch('api.settings.get_settings_service')
        mock_service.return_value.update_settings.return_value = updated

        response = client.put('/api/settings', json=self._settings_payload('/tmp/pyscf-updated'))

        assert response.status_code == 200
        assert response.json()['success'] is True
        assert response.json()['data']['settings'] == updated
        mock_service.return_value.update_settings.assert_called_once()
```

- [ ] **Step 3: Convert service-backed routers with explicit body models**

For `src/python/api/settings.py`, use:

```python
from fastapi import APIRouter

from generated_models import SettingsUpdateRequest
from services import get_settings_service

router = APIRouter(prefix='/api/settings')


@router.get('')
def get_settings():
    settings = get_settings_service().get_settings()
    return {'success': True, 'data': {'settings': settings}}


@router.put('')
def update_settings(body: SettingsUpdateRequest):
    new_settings = body.root if hasattr(body, 'root') else body
    updated_settings = get_settings_service().update_settings(new_settings.model_dump())
    return {'success': True, 'data': {'settings': updated_settings}}
```

For `src/python/api/pubchem.py`, use:

```python
from fastapi import APIRouter

from generated_models import PubChemSearchRequest, XYZValidateRequest
from services import get_pubchem_service

router = APIRouter(prefix='/api/pubchem')


@router.post('/search')
def search_pubchem(body: PubChemSearchRequest):
    search_type_value = body.searchType or "name"
    search_type = search_type_value.value if hasattr(search_type_value, "value") else str(search_type_value)
    result = get_pubchem_service().search_compound(body.query, search_type)
    return {'success': True, 'data': result}


@router.post('/validate')
def validate_xyz_endpoint(body: XYZValidateRequest):
    validation_result = get_pubchem_service().validate_xyz(body.xyz)
    return {'success': True, 'data': validation_result}
```

For `src/python/api/smiles.py`, use:

```python
from fastapi import APIRouter

from generated_models import SMILESConvertRequest
from services import get_smiles_service

router = APIRouter(prefix='/api/smiles')


@router.post('/convert')
def convert_smiles(body: SMILESConvertRequest):
    result = get_smiles_service().convert_smiles(body.smiles)
    return {'success': True, 'data': result}
```

For `src/python/api/chat_history.py`, use:

```python
from fastapi import APIRouter
from fastapi.responses import JSONResponse

from generated_models import CreateChatSessionRequest, UpdateChatSessionRequest
from services.chat_history_service import get_chat_history_service

router = APIRouter(prefix='/api/chat-history')


@router.get('/sessions')
def get_chat_sessions():
    data = get_chat_history_service().list_sessions()
    return {'success': True, 'data': data}


@router.post('/sessions', status_code=201)
def create_chat_session(body: CreateChatSessionRequest):
    session = get_chat_history_service().create_session(name=body.name)
    return {'success': True, 'data': {'session': session}}


@router.get('/sessions/{session_id}')
def get_chat_session(session_id: str):
    session_data = get_chat_history_service().get_session_with_messages(session_id)
    if session_data is None:
        return JSONResponse(
            {'success': False, 'error': f'Chat session not found: {session_id}'},
            status_code=404,
        )
    return {'success': True, 'data': session_data}


@router.patch('/sessions/{session_id}')
def update_chat_session(session_id: str, body: UpdateChatSessionRequest):
    session = get_chat_history_service().update_session(session_id, body.name)
    if session is None:
        return JSONResponse(
            {'success': False, 'error': f'Chat session not found: {session_id}'},
            status_code=404,
        )
    return {'success': True, 'data': {'session': session}}


@router.delete('/sessions/{session_id}')
def delete_chat_session(session_id: str):
    deleted = get_chat_history_service().delete_session(session_id)
    if not deleted:
        return JSONResponse(
            {'success': False, 'error': f'Chat session not found: {session_id}'},
            status_code=404,
        )
    return {
        'success': True,
        'data': {
            'message': 'Chat session deleted successfully',
            'deleted_id': session_id,
        },
    }
```

For service errors, let `ServiceError` propagate to the handler from Task 3.

- [ ] **Step 4: Extend `api/__init__.py` with the converted routers**

Add imports:

```python
from .chat_history import router as chat_history_router
from .pubchem import router as pubchem_router
from .settings import router as settings_router
from .smiles import router as smiles_router
```

Add includes after `health_router` and before `swagger_router`:

```python
app.include_router(pubchem_router)
app.include_router(smiles_router)
app.include_router(settings_router)
app.include_router(chat_history_router)
```

- [ ] **Step 5: Convert endpoint tests from `get_json()` to `.json()`**

Replace patterns like:

```python
data = response.get_json()
```

with:

```python
data = response.json()
```

Replace Flask test client calls that pass `data=` JSON strings with FastAPI TestClient calls:

```python
response = client.post('/api/pubchem/validate', json={'xyz': 'H 0 0 0'})
```

Also convert the remaining Flask response/client attributes in endpoint and auth tests:

```python
response.content_type
response.data
response.json['success']
client.post('/api/pubchem/search', data='invalid json', content_type='application/json')
```

becomes:

```python
response.headers['content-type']
response.content
response.json()['success']
client.post(
    '/api/pubchem/search',
    content=b'invalid json',
    headers={'Content-Type': 'application/json'},
)
```

Keep all endpoint path and envelope assertions unchanged.

- [ ] **Step 6: Run simple endpoint tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_health_api.py tests/integration/test_api_endpoints/test_settings_api.py tests/integration/test_api_endpoints/test_pubchem_api.py tests/integration/test_api_endpoints/test_smiles_api.py -v
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/python/api/__init__.py src/python/api/settings.py src/python/api/pubchem.py src/python/api/smiles.py src/python/api/chat_history.py src/python/tests/integration/test_api_endpoints
git commit -m "feat: convert simple API routers to FastAPI"
```

---

### Task 5: Convert Quantum And System Routers

**Files:**
- Modify: `src/python/api/__init__.py`
- Modify: `src/python/api/quantum.py`
- Modify: `src/python/api/system.py`
- Modify: `src/python/tests/integration/test_api_endpoints/test_quantum_api.py`
- Modify: `src/python/tests/integration/test_validation_fastapi.py`
- Modify: `src/python/tests/integration/test_api_endpoints/test_system_api.py`

- [ ] **Step 1: Add failing regression tests for quantum raw-body order and GPU loopback guard**

Extend `src/python/tests/integration/test_validation_fastapi.py`:

```python
def test_quantum_calculate_requires_calculation_method_before_model_validation(client):
    response = client.post('/api/quantum/calculate', json={'molecule': {'atoms': []}})

    assert response.status_code == 400
    body = response.json()
    assert body['success'] is False
    assert 'calculation_method' in body['error']


def test_quantum_calculate_rejects_inapplicable_method_parameter(client):
    response = client.post(
        '/api/quantum/calculate',
        json={
            'calculation_method': 'HF',
            'xyz': 'H 0 0 0',
            'basis_function': 'sto-3g',
            'exchange_correlation': 'b3lyp',
        },
    )

    assert response.status_code == 400
    body = response.json()
    assert body['success'] is False
    assert 'exchange_correlation' in body['error']


def test_gpu4pyscf_install_rejects_non_loopback_client(client, monkeypatch):
    import api.system as system_api

    monkeypatch.setattr(system_api, '_get_client_host', lambda request: '203.0.113.10')

    response = client.post('/api/system/gpu4pyscf-install')

    assert response.status_code == 403
    assert response.json() == {
        'success': False,
        'error': 'GPU4PySCF installation is only available from the local machine.',
    }
```

- [ ] **Step 2: Run quantum/system tests to verify they fail**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_validation_fastapi.py tests/integration/test_api_endpoints/test_quantum_api.py tests/integration/test_api_endpoints/test_system_api.py -v
```

Expected: FAIL because routers still import Flask and `flask_pydantic`.

- [ ] **Step 3: Convert `quantum.py` to FastAPI while preserving validation order**

Use this shape for `/api/quantum/calculate`:

```python
import json
from datetime import datetime

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import JSONResponse

from generated_models import CalculationUpdateRequest, QuantumCalculationRequest
from quantum_calc.method_defaults import validate_parameters_for_method
from services import get_quantum_service

router = APIRouter(prefix='/api/quantum')


@router.post('/calculate')
async def quantum_calculate(request: Request):
    raw_body = await request.body()
    if not raw_body:
        raise HTTPException(status_code=400, detail='Request body is required')

    raw_data = json.loads(raw_body)
    if not raw_data:
        raise HTTPException(status_code=400, detail='Request body is required')

    calculation_method = raw_data.get('calculation_method')
    if not calculation_method:
        raise HTTPException(status_code=400, detail='calculation_method is required')

    is_valid, applicability_error = validate_parameters_for_method(calculation_method, raw_data)
    if not is_valid:
        return JSONResponse({'success': False, 'error': applicability_error}, status_code=400)

    body = QuantumCalculationRequest.model_validate(raw_data)
    validated_model = body.root if hasattr(body, 'root') else body
    parameters = validated_model.model_dump(exclude_none=False, mode='python')
    for key, value in list(parameters.items()):
        parameters[key] = value.value if hasattr(value, 'value') else value
    parameters['created_at'] = datetime.now().isoformat()

    result = get_quantum_service().start_calculation(parameters)
    return JSONResponse({'success': True, 'data': {'calculation': result}}, status_code=202)
```

Convert the remaining quantum routes with explicit FastAPI path/query/body parameters:

```python
from typing import Annotated

from fastapi import Query


@router.get('/supported-parameters')
def get_supported_parameters():
    parameters = get_quantum_service().get_supported_parameters()
    return {'success': True, 'data': parameters}


@router.get('/calculations')
def list_calculations(
    name_query: str | None = None,
    status: str | None = None,
    calculation_method: str | None = None,
    basis_function: str | None = None,
    date_from: str | None = None,
    date_to: str | None = None,
):
    result = get_quantum_service().list_calculations(
        name_query=name_query,
        status=status,
        calculation_method=calculation_method,
        basis_function=basis_function,
        date_from=date_from,
        date_to=date_to,
    )
    return {'success': True, 'data': result}


@router.get('/status')
def get_calculation_status():
    result = get_quantum_service().get_calculation_status()
    return {'success': True, 'data': result}


@router.get('/calculations/{calculation_id}')
def get_calculation_details(calculation_id: str):
    result = get_quantum_service().get_calculation_details(calculation_id)
    return {'success': True, 'data': result}


@router.put('/calculations/{calculation_id}')
def update_calculation(calculation_id: str, body: CalculationUpdateRequest):
    result = get_quantum_service().update_calculation(calculation_id, body.name)
    return {'success': True, 'data': result}


@router.post('/calculations/{calculation_id}/pause')
def pause_calculation(calculation_id: str):
    result = get_quantum_service().pause_calculation(calculation_id)
    return JSONResponse({'success': True, 'data': result}, status_code=202)


@router.post('/calculations/{calculation_id}/resume')
def resume_calculation(calculation_id: str):
    result = get_quantum_service().resume_calculation(calculation_id)
    return JSONResponse({'success': True, 'data': result}, status_code=202)


@router.delete('/calculations/{calculation_id}')
def delete_calculation(calculation_id: str):
    result = get_quantum_service().delete_calculation(calculation_id)
    return {'success': True, 'data': result}


@router.get('/calculations/{calculation_id}/orbitals')
def get_orbitals(calculation_id: str):
    orbital_summary = get_quantum_service().get_molecular_orbitals(calculation_id)
    return {'success': True, 'data': orbital_summary}


@router.get('/calculations/{calculation_id}/orbitals/{orbital_index}/cube')
def get_orbital_cube(
    calculation_id: str,
    orbital_index: int,
    gridSize: Annotated[int, Query(alias='gridSize')] = 80,
    isovaluePos: Annotated[float | None, Query(alias='isovaluePos')] = None,
    isovalueNeg: Annotated[float | None, Query(alias='isovalueNeg')] = None,
):
    cube_data = get_quantum_service().generate_orbital_cube(
        calculation_id,
        orbital_index,
        grid_size=gridSize,
        isovalue_pos=isovaluePos,
        isovalue_neg=isovalueNeg,
    )
    return {'success': True, 'data': cube_data}


@router.get('/calculations/{calculation_id}/orbitals/cube-files')
def list_cube_files(calculation_id: str):
    result = get_quantum_service().list_cube_files(calculation_id)
    return {'success': True, 'data': result}


@router.delete('/calculations/{calculation_id}/orbitals/cube-files')
def delete_cube_files(calculation_id: str, orbital_index: int | None = None):
    result = get_quantum_service().delete_cube_files(calculation_id, orbital_index)
    return {'success': True, 'data': result}


@router.get('/calculations/{calculation_id}/ir-spectrum')
def get_ir_spectrum(
    calculation_id: str,
    broadening_fwhm: float = 100.0,
    x_min: float = 400.0,
    x_max: float = 4000.0,
    show_peaks: bool = True,
):
    result = get_quantum_service().generate_ir_spectrum(
        calculation_id,
        broadening_fwhm=broadening_fwhm,
        x_min=x_min,
        x_max=x_max,
        show_peaks=show_peaks,
    )
    return {'success': True, 'data': result}
```

Check these paths against `src/api-spec/openapi.yaml` before committing.

- [ ] **Step 4: Convert `system.py` to FastAPI and keep loopback helper logic**

Add an injectable client host helper:

```python
import json
from datetime import datetime

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from generated_models import AllocatedResources, ResourceConstraints, SystemResourceInfo, SystemResourceSummary
from services import get_system_service

router = APIRouter()


def _get_client_host(request: Request) -> str | None:
    return request.client.host if request.client else None
```

Keep `_is_loopback_address` behavior for IPv4, IPv6, and IPv4-mapped IPv6. Use it in the GPU installer route:

```python
@router.post('/api/system/gpu4pyscf-install')
async def install_gpu4pyscf(request: Request):
    client_host = _get_client_host(request)
    if not _is_loopback_address(client_host):
        return JSONResponse(
            {
                'success': False,
                'error': 'GPU4PySCF installation is only available from the local machine.',
            },
            status_code=403,
        )

    raw_body = await request.body()
    payload = json.loads(raw_body) if raw_body else {}
    include_cutensor = bool(payload.get('include_cutensor', True))
    force_reinstall = bool(payload.get('force_reinstall', False))
    result = get_system_service().install_gpu4pyscf(
        include_cutensor=include_cutensor,
        force_reinstall=force_reinstall,
    )
    return {'success': True, 'data': result}
```

Convert diagnostics routes:

```python
@router.get('/api/system/resource-status')
def get_system_resource_status():
    resource_summary = get_system_service().get_resource_status()
    system_info = SystemResourceInfo(
        total_cpu_cores=resource_summary['system_info']['total_cpu_cores'],
        total_memory_mb=resource_summary['system_info']['total_memory_mb'],
        available_memory_mb=resource_summary['system_info']['available_memory_mb'],
        cpu_usage_percent=resource_summary['system_info']['cpu_usage_percent'],
        memory_usage_percent=resource_summary['system_info']['memory_usage_percent'],
        timestamp=datetime.fromisoformat(resource_summary['system_info']['timestamp'].replace('Z', '+00:00')),
    )
    constraints = ResourceConstraints(**resource_summary['resource_constraints'])
    allocated = AllocatedResources(**resource_summary['allocated_resources'])
    summary = SystemResourceSummary(
        system_info=system_info,
        resource_constraints=constraints,
        allocated_resources=allocated,
    )
    return {'success': True, 'data': summary.model_dump(mode='json')}


@router.get('/api/system/gpu4pyscf-status')
def get_gpu4pyscf_status():
    status = get_system_service().get_gpu4pyscf_status()
    return {'success': True, 'data': status}

@router.get('/api/debug/system-diagnostics')
def get_system_diagnostics():
    diagnostics = get_system_service().get_system_diagnostics()
    return {'success': True, 'data': diagnostics}

@router.get('/api/debug/process-manager-diagnostics')
def get_process_manager_diagnostics():
    diagnostics = get_system_service().get_process_manager_diagnostics()
    return {'success': True, 'data': diagnostics}

@router.get('/api/debug/resource-manager-diagnostics')
def get_resource_manager_diagnostics():
    diagnostics = get_system_service().get_resource_manager_diagnostics()
    return {'success': True, 'data': diagnostics}
```

- [ ] **Step 5: Extend `api/__init__.py` with quantum and system routers**

Add imports:

```python
from .quantum import router as quantum_router
from .system import router as system_router
```

Add includes after `settings_router`:

```python
app.include_router(system_router)
app.include_router(quantum_router)
```

- [ ] **Step 6: Run quantum/system tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_validation_fastapi.py tests/integration/test_api_endpoints/test_quantum_api.py tests/integration/test_api_endpoints/test_system_api.py -v
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/python/api/__init__.py src/python/api/quantum.py src/python/api/system.py src/python/tests/integration/test_validation_fastapi.py src/python/tests/integration/test_api_endpoints/test_quantum_api.py src/python/tests/integration/test_api_endpoints/test_system_api.py
git commit -m "feat: convert quantum and system APIs to FastAPI"
```

---

### Task 6: Convert Agent SSE Route

**Files:**
- Modify: `src/python/api/__init__.py`
- Modify: `src/python/api/agent.py`
- Create or modify: `src/python/tests/integration/test_api_endpoints/test_agent_api.py`

- [ ] **Step 1: Write failing SSE route test**

Create or extend `src/python/tests/integration/test_api_endpoints/test_agent_api.py`:

```python
import json


def test_agent_chat_stream_preserves_sse_format(client, mocker):
    def fake_stream(*args, **kwargs):
        yield {'type': 'agent_status', 'payload': {'status': 'started'}}
        yield {'type': 'chunk', 'payload': {'text': 'hello'}}
        yield {'type': 'done', 'payload': {'message_id': 'msg-1'}}

    mocker.patch('api.agent.stream_chat_response', side_effect=fake_stream)

    with client.stream(
        'POST',
        '/api/agent/chat',
        json={'message': 'Hello', 'history': [], 'session_id': 'session-1'},
    ) as response:
        chunks = ''.join(response.iter_text())

    assert response.status_code == 200
    assert response.headers['content-type'].startswith('text/event-stream')
    assert 'data: {"type": "agent_status", "payload": {"status": "started"}}\n\n' in chunks
    assert 'data: {"type": "chunk", "payload": {"text": "hello"}}\n\n' in chunks
    assert 'data: {"type": "done", "payload": {"message_id": "msg-1"}}\n\n' in chunks
```

- [ ] **Step 2: Run SSE test to verify it fails**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_agent_api.py -v
```

Expected: FAIL because the route still uses Flask `Response` and `stream_with_context`.

- [ ] **Step 3: Convert `agent.py` to `StreamingResponse`**

Replace Flask imports with:

```python
import json

from fastapi import APIRouter
from fastapi.responses import StreamingResponse

from generated_models import AgentChatRequest
```

Keep `_format_sse_event`:

```python
def _format_sse_event(event: dict) -> str:
    return f"data: {json.dumps(event)}\n\n"
```

Use a synchronous generator wrapped by `StreamingResponse`:

```python
router = APIRouter(prefix='/api/agent')


@router.post('/chat')
def chat(request: AgentChatRequest):
    def event_generator():
        for event in stream_chat_response(request):
            yield _format_sse_event(event)

    return StreamingResponse(event_generator(), media_type='text/event-stream')
```

Preserve existing chat history persistence and Gemini call behavior by moving the current Flask route body into `stream_chat_response(request)`.

- [ ] **Step 4: Extend `api/__init__.py` with the agent router**

Add import:

```python
from .agent import router as agent_router
```

Add include after `quantum_router`:

```python
app.include_router(agent_router)
```

- [ ] **Step 5: Run SSE test**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_agent_api.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/python/api/__init__.py src/python/api/agent.py src/python/tests/integration/test_api_endpoints/test_agent_api.py
git commit -m "feat: serve agent chat SSE with FastAPI"
```

---

### Task 7: Remove Flask Config Coupling

**Files:**
- Modify: `src/python/quantum_calc/config_manager.py`
- Modify: `src/python/api/health.py`
- Modify: `src/python/tests/integration/test_api_endpoints/test_health_api.py`
- Create: `src/python/tests/unit/test_quantum_config_manager_fastapi.py`

- [ ] **Step 1: Write failing config-manager tests**

Create `src/python/tests/unit/test_quantum_config_manager_fastapi.py`:

```python
from config import get_server_config
from quantum_calc.config_manager import QuantumCalculationConfigManager


def test_quantum_config_manager_does_not_import_flask():
    import quantum_calc.config_manager as module

    assert not hasattr(module, 'current_app')


def test_quantum_config_manager_reads_server_config():
    manager = QuantumCalculationConfigManager()
    config = manager._get_config()

    server_config = get_server_config()
    expected = server_config.get('quantum_calculation_defaults', {})
    assert config['quantum_calculation_defaults'] == expected
    assert manager.get_memory_setting('DFT') == int(expected['memory_settings']['DFT'])
    assert manager.get_config_info() == {
        'config_source': 'ServerConfig',
        'config_available': True,
        'using_fallback': False,
    }
```

- [ ] **Step 2: Run config tests to verify they fail**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_quantum_config_manager_fastapi.py -v
```

Expected: FAIL because `quantum_calc.config_manager` imports Flask `current_app`.

- [ ] **Step 3: Refactor `config_manager.py`**

Remove:

```python
from flask import current_app
```

Add:

```python
from config import get_server_config
```

Replace Flask config reads with:

```python
def _get_config(self) -> Dict[str, Any]:
    server_config = get_server_config()
    config = server_config.get('quantum_calculation_defaults', {})
    if config:
        return {'quantum_calculation_defaults': config}
    logger.warning("Quantum calculation defaults not found in ServerConfig, using fallback")
    return self._get_fallback_config()
```

Update `get_config_info` so it reports:

```python
"config_source": "ServerConfig"
```

Do not retain a Flask fallback branch.

- [ ] **Step 4: Run config and health tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_quantum_config_manager_fastapi.py tests/integration/test_api_endpoints/test_health_api.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/python/quantum_calc/config_manager.py src/python/api/health.py src/python/tests/unit/test_quantum_config_manager_fastapi.py src/python/tests/integration/test_api_endpoints/test_health_api.py
git commit -m "refactor: remove Flask config coupling"
```

---

### Task 8: Socket.IO Async Handlers And Event Loop Bridge

**Files:**
- Modify: `src/python/websocket/event_loop_bridge.py`
- Modify: `src/python/websocket/handlers.py`
- Modify: `src/python/services/notification_service.py`
- Create: `src/python/tests/unit/test_event_loop_bridge.py`
- Modify: `src/python/tests/integration/test_websocket_handlers.py`
- Modify: `src/python/tests/integration/test_calculation_workflow.py`

- [ ] **Step 1: Write failing event-loop bridge tests**

Create `src/python/tests/unit/test_event_loop_bridge.py`:

```python
import asyncio
import threading

from websocket.event_loop_bridge import bind_event_loop, clear_event_loop, schedule_coroutine


def test_schedule_coroutine_from_background_thread():
    async def main():
        loop = asyncio.get_running_loop()
        bind_event_loop(loop)
        results = []

        async def append_value(value):
            results.append(value)
            return value

        def worker():
            future = schedule_coroutine(append_value('sent'))
            assert future.result(timeout=2) == 'sent'

        thread = threading.Thread(target=worker)
        thread.start()
        await asyncio.to_thread(thread.join, 2)

        assert not thread.is_alive()
        assert results == ['sent']
        clear_event_loop()

    asyncio.run(main())


def test_schedule_coroutine_without_bound_loop_returns_none(caplog):
    clear_event_loop()

    async def noop():
        return None

    future = schedule_coroutine(noop())

    assert future is None
    assert 'ASGI event loop is not bound' in caplog.text
```

- [ ] **Step 2: Run bridge tests to verify they fail**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_event_loop_bridge.py -v
```

Expected: FAIL because `schedule_coroutine` is not implemented yet.

- [ ] **Step 3: Extend `event_loop_bridge.py` with callback scheduling**

Replace the minimal Task 2 bridge with:

```python
import asyncio
import logging
from concurrent.futures import Future
from typing import Coroutine, TypeVar

logger = logging.getLogger(__name__)
T = TypeVar('T')

_loop: asyncio.AbstractEventLoop | None = None


def bind_event_loop(loop: asyncio.AbstractEventLoop) -> None:
    global _loop
    _loop = loop
    logger.info("ASGI event loop bound for background notifications")


def clear_event_loop() -> None:
    global _loop
    _loop = None
    logger.info("ASGI event loop cleared for background notifications")


def schedule_coroutine(coro: Coroutine[object, object, T]) -> Future[T] | None:
    if _loop is None or _loop.is_closed():
        logger.warning("ASGI event loop is not bound; dropping scheduled coroutine")
        coro.close()
        return None

    future = asyncio.run_coroutine_threadsafe(coro, _loop)

    def log_failure(done_future: Future[T]) -> None:
        try:
            done_future.result()
        except Exception:
            logger.exception("Scheduled ASGI coroutine failed")

    future.add_done_callback(log_failure)
    return future
```

- [ ] **Step 4: Convert websocket handlers to async server APIs**

In `src/python/websocket/handlers.py`, remove Flask imports:

```python
from flask import session
from flask_socketio import emit, join_room, leave_room
```

Use per-SID state and keep the existing file-watcher semantics:

```python
import os
from typing import Any

from quantum_calc import CalculationRepository, get_websocket_watcher
from websocket.event_loop_bridge import schedule_coroutine

_sid_state: dict[str, dict[str, Any]] = {}


def _state_for(sid: str) -> dict[str, Any]:
    return _sid_state.setdefault(sid, {'callbacks': {}})
```

Register async handlers:

```python
def register_websocket_handlers(sio):
    if getattr(sio, '_pyscf_handlers_registered', False):
        return
    setattr(sio, '_pyscf_handlers_registered', True)

    @sio.event
    async def connect(sid, environ, auth):
        expected_token = os.getenv('PYSCF_AUTH_TOKEN')
        if expected_token and (not auth or auth.get('token') != expected_token):
            return False
        _state_for(sid)
        return True

    @sio.event
    async def join_global_updates(sid):
        await sio.enter_room(sid, 'global_updates')
        logger.info("Client joined global_updates room for real-time monitoring of all calculations")

    @sio.event
    async def leave_global_updates(sid):
        await sio.leave_room(sid, 'global_updates')
        logger.info("Client left global_updates room")

    @sio.event
    async def join_calculation(sid, data):
        calculation_id = (data or {}).get('calculation_id')
        if not calculation_id:
            await sio.emit('error', {'error': 'calculation_id is required'}, to=sid)
            return

        if calculation_id.startswith('new-calculation-'):
            logger.info("Socket.IO connection attempt for temporary calculation ID: %s", calculation_id)
        else:
            logger.info("Socket.IO connection established for calculation %s", calculation_id)

        from quantum_calc import get_current_settings

        settings = get_current_settings()
        file_manager = CalculationRepository(base_dir=settings.calculations_directory)
        try:
            calc_path = str(file_manager.resolve_calculation_path(calculation_id))
        except ValueError:
            await sio.emit(
                'error',
                {
                    'error': 'Invalid calculation ID.',
                    'id': calculation_id,
                    'is_temporary': calculation_id.startswith('new-calculation-'),
                },
                to=sid,
            )
            return

        if not os.path.isdir(calc_path):
            if calculation_id.startswith('new-calculation-'):
                error_message = f'Temporary calculation ID "{calculation_id}" does not exist on server.'
            else:
                error_message = f'Calculation "{calculation_id}" not found.'
            await sio.emit(
                'error',
                {
                    'error': error_message,
                    'id': calculation_id,
                    'is_temporary': calculation_id.startswith('new-calculation-'),
                },
                to=sid,
            )
            return

        room = f'calculation_{calculation_id}'
        await sio.enter_room(sid, room)
        state = _state_for(sid)

        def on_file_change(file_data: dict) -> None:
            async def emit_update() -> None:
                try:
                    calculation_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
                    await sio.emit('calculation_update', calculation_instance, room=room)
                    if calculation_instance['status'] in ['completed', 'error']:
                        logger.info(
                            "Calculation %s finished with status '%s'.",
                            calculation_id,
                            calculation_instance['status'],
                        )
                except Exception:
                    logger.exception("Error in file change callback for %s", calculation_id)
                    await sio.emit(
                        'error',
                        {'error': 'Failed to read calculation data', 'id': calculation_id},
                        room=room,
                    )

            schedule_coroutine(emit_update())

        try:
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.add_connection(calculation_id, on_file_change)
            state['callbacks'][calculation_id] = on_file_change
            initial_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
            await sio.emit('calculation_update', initial_instance, to=sid)
        except Exception:
            logger.exception("Error setting up Socket.IO monitoring for %s", calculation_id)
            await sio.emit(
                'error',
                {'error': 'Failed to set up calculation monitoring', 'id': calculation_id},
                to=sid,
            )

    @sio.event
    async def leave_calculation(sid, data):
        calculation_id = (data or {}).get('calculation_id')
        if not calculation_id:
            await sio.emit('error', {'error': 'calculation_id is required'}, to=sid)
            return
        room = f'calculation_{calculation_id}'
        await sio.leave_room(sid, room)
        state = _state_for(sid)
        callback = state['callbacks'].pop(calculation_id, None)
        if callback is not None:
            from quantum_calc import get_current_settings

            settings = get_current_settings()
            file_manager = CalculationRepository(base_dir=settings.calculations_directory)
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.remove_connection(calculation_id, callback)
        logger.info("Client left calculation %s", calculation_id)

    @sio.event
    async def disconnect(sid):
        state = _sid_state.pop(sid, {'callbacks': {}})
        if not state['callbacks']:
            return

        from quantum_calc import get_current_settings

        settings = get_current_settings()
        file_manager = CalculationRepository(base_dir=settings.calculations_directory)
        watcher = get_websocket_watcher(file_manager.get_base_directory())
        for calculation_id, callback in state['callbacks'].items():
            watcher.remove_connection(calculation_id, callback)
            logger.info("Cleaned up file watcher for disconnected client (calculation %s)", calculation_id)
```

The callback must build and emit the same `CalculationInstance` payload returned by `build_calculation_instance`; do not replace it with a reduced status-only payload.

- [ ] **Step 5: Update notification service**

In `src/python/services/notification_service.py`, keep the public synchronous method and schedule async emits:

```python
from websocket.event_loop_bridge import schedule_coroutine


class NotificationService:
    def __init__(self):
        self._socketio = None
        logger.info("NotificationService initialized (socketio not yet bound)")

    def bind_socketio(self, socketio) -> None:
        self._socketio = socketio
        logger.info("Socket.IO AsyncServer bound to NotificationService")

    def send_calculation_update(
        self,
        calculation_id: str,
        status: str,
        error_message: str | None = None,
    ) -> None:
        if self._socketio is None:
            logger.warning("Cannot send notification for %s: Socket.IO not bound", calculation_id)
            return

        try:
            from quantum_calc import CalculationRepository, get_current_settings
            import os
            from datetime import datetime

            settings = get_current_settings()
            file_manager = CalculationRepository(base_dir=settings.calculations_directory)
            calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
            if not os.path.exists(calc_dir):
                logger.warning("Calculation directory not found: %s", calc_dir)
                return

            parameters = file_manager.read_calculation_parameters(calc_dir) or {}
            results = file_manager.read_calculation_results(calc_dir)
            display_name = file_manager.get_display_name(calculation_id, parameters)
            calculation_instance = {
                'id': calculation_id,
                'name': display_name,
                'status': status,
                'createdAt': parameters.get('created_at', datetime.now().isoformat()),
                'updatedAt': datetime.now().isoformat(),
                'parameters': parameters,
                'results': results,
                'workingDirectory': calc_dir,
            }
            if error_message:
                calculation_instance['error'] = error_message
                calculation_instance['errorMessage'] = error_message
        except Exception:
            logger.exception("Error building WebSocket notification payload")
            return

        schedule_coroutine(
            self._socketio.emit(
                'calculation_update',
                calculation_instance,
                room='global_updates',
            )
        )
```

- [ ] **Step 6: Convert workflow tests away from Flask-SocketIO fixture**

In `src/python/tests/integration/test_calculation_workflow.py` and `src/python/tests/integration/test_pause_resume_workflow.py`, replace `.get_json()` with `.json()` for every FastAPI response:

```python
submit_data = response_submit.json()
details_data = response_details.json()
calc_id = response.get_json()['data']['calculation']['id']
```

becomes:

```python
submit_data = response_submit.json()
details_data = response_details.json()
calc_id = response.json()['data']['calculation']['id']
```

Replace the Flask-SocketIO test-client workflow test with an async handler or ASGI smoke assertion. A minimal conversion is to remove the `socketio_client` fixture argument and verify notification delivery through the real ASGI smoke test added in Task 9:

```python
def test_workflow_with_websocket_integration(client, mocker, valid_hf_params):
    mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
    mock_mol = mocker.MagicMock()
    mock_scf = mocker.MagicMock()
    mock_scf.kernel.return_value = -1.06
    mock_scf.mo_energy = [-0.5, 0.3]
    mock_scf.mo_occ = [2.0, 0.0]
    mocker.patch('quantum_calc.hf_calculator.gto.M', return_value=mock_mol)
    mocker.patch('quantum_calc.hf_calculator.scf.RHF', return_value=mock_scf)

    response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)

    assert response_submit.status_code == 202
    calc_id = response_submit.json()['data']['calculation']['id']
    assert calc_id
```

Task 9 is responsible for proving actual Socket.IO ASGI `join_calculation` and `calculation_update` delivery.

In `src/python/tests/integration/test_api_endpoints/test_system_api.py`, convert Flask-only response/client APIs:

```python
response_text = response.text
data = response.json()
```

For loopback and remote-address tests, patch the helper instead of passing Flask `environ_base`:

```python
mock_get_host = mocker.patch('api.system._get_client_host', return_value='127.0.0.1')
response = client.post('/api/system/gpu4pyscf-install')
mock_get_host.assert_called_once()
```

and:

```python
mocker.patch('api.system._get_client_host', return_value='10.10.10.10')
response = client.post('/api/system/gpu4pyscf-install')
```

- [ ] **Step 7: Run bridge tests, handler-level tests, and workflow tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/unit/test_event_loop_bridge.py tests/integration/test_websocket_handlers.py tests/integration/test_calculation_workflow.py tests/integration/test_pause_resume_workflow.py tests/integration/test_api_endpoints/test_system_api.py -v
```

Expected: PASS after handler tests are converted to async helpers or moved to the ASGI smoke test in Task 9.

- [ ] **Step 8: Commit**

```bash
git add src/python/websocket/event_loop_bridge.py src/python/websocket/handlers.py src/python/services/notification_service.py src/python/tests/unit/test_event_loop_bridge.py src/python/tests/integration/test_websocket_handlers.py src/python/tests/integration/test_calculation_workflow.py src/python/tests/integration/test_pause_resume_workflow.py src/python/tests/integration/test_api_endpoints/test_system_api.py
git commit -m "feat: migrate Socket.IO handlers to ASGI"
```

---

### Task 9: Real Socket.IO ASGI Smoke Test

**Files:**
- Create: `src/python/tests/integration/test_socketio_asgi.py`
- Modify: `src/python/tests/conftest.py`

- [ ] **Step 1: Add ASGI server fixture**

In `src/python/tests/conftest.py`, add:

```python
import socket
import threading
import time
from pathlib import Path
from unittest import mock

import services as services_module
import uvicorn
import quantum_calc.settings_manager as settings_manager_module
from quantum_calc.settings_manager import SettingsManager

from app import create_app


def _get_free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(('127.0.0.1', 0))
        return sock.getsockname()[1]


@pytest.fixture
def asgi_server(monkeypatch, tmp_path):
    port = _get_free_port()
    monkeypatch.setenv('PYSCF_AUTH_TOKEN', 'socket-token')
    monkeypatch.setenv('PYSCF_ENV', 'development')
    calculations_dir = tmp_path / 'socket-calculations'
    settings_file = tmp_path / 'app-settings.json'
    test_settings_manager = SettingsManager(settings_file=str(settings_file))
    test_settings = test_settings_manager.get_default_settings().model_copy(
        update={'calculations_directory': str(calculations_dir)}
    )
    test_settings_manager.save_settings(test_settings)

    with (
        mock.patch.object(settings_manager_module, '_settings_manager', test_settings_manager),
        mock.patch.multiple(
            services_module,
            _quantum_service=None,
            _pubchem_service=None,
            _smiles_service=None,
            _settings_service=None,
            _system_service=None,
        ),
    ):
        asgi_app = create_app(
            server_port=port,
            test_config={
                'TESTING': True,
                'CALCULATIONS_DIR': str(calculations_dir),
                'WEBSOCKET_WATCHER_ENABLED': False,
                'SOCKETIO': {
                    'cors_allowed_origins': ['http://127.0.0.1:3000', 'http://localhost:3000', 'file://'],
                    'logger': False,
                    'engineio_logger': False,
                },
            },
        )
        config = uvicorn.Config(asgi_app, host='127.0.0.1', port=port, log_level='warning')
        server = uvicorn.Server(config)
        thread = threading.Thread(target=server.run, daemon=True)
        thread.start()

        deadline = time.time() + 5
        while not server.started and time.time() < deadline:
            time.sleep(0.05)

        assert server.started
        yield f'http://127.0.0.1:{port}'

        server.should_exit = True
        thread.join(timeout=5)
```

- [ ] **Step 2: Write real Socket.IO smoke tests**

Create `src/python/tests/integration/test_socketio_asgi.py`:

```python
import asyncio
from pathlib import Path

import pytest
import socketio

from services.notification_service import get_notification_service


def test_socketio_rejects_missing_token(asgi_server):
    async def scenario():
        client = socketio.AsyncClient()
        await client.connect(
            asgi_server,
            transports=['websocket'],
            headers={'Origin': 'http://127.0.0.1:3000'},
        )

    with pytest.raises(socketio.exceptions.ConnectionError):
        asyncio.run(scenario())


def test_socketio_rejects_wrong_token(asgi_server):
    async def scenario():
        client = socketio.AsyncClient()
        await client.connect(
            asgi_server,
            transports=['websocket'],
            auth={'token': 'wrong-token'},
            headers={'Origin': 'http://127.0.0.1:3000'},
        )

    with pytest.raises(socketio.exceptions.ConnectionError):
        asyncio.run(scenario())


def test_socketio_accepts_token_and_receives_calculation_update(asgi_server):
    async def scenario():
        client = socketio.AsyncClient()
        received = []

        @client.on('calculation_update')
        async def on_calculation_update(data):
            received.append(data)

        from quantum_calc import CalculationRepository, get_current_settings

        settings = get_current_settings()
        repository = CalculationRepository(base_dir=settings.calculations_directory)
        calc_dir = Path(repository.get_base_directory()) / 'calc-123'
        calc_dir.mkdir(parents=True, exist_ok=True)
        (calc_dir / 'parameters.json').write_text(
            '{"name":"Socket Smoke","calculation_method":"HF","basis_function":"sto-3g","created_at":"2024-01-01T00:00:00"}',
            encoding='utf-8',
        )
        (calc_dir / 'status.json').write_text('{"status":"completed"}', encoding='utf-8')

        await client.connect(
            asgi_server,
            transports=['websocket'],
            auth={'token': 'socket-token'},
            headers={'Origin': 'http://127.0.0.1:3000'},
        )
        await client.emit('join_global_updates')
        await client.emit('join_calculation', {'calculation_id': 'calc-123'})

        deadline = asyncio.get_running_loop().time() + 2
        while not received and asyncio.get_running_loop().time() < deadline:
            await asyncio.sleep(0.05)
        assert received
        received.clear()

        service = get_notification_service()
        service.send_calculation_update('calc-123', 'completed')

        deadline = asyncio.get_running_loop().time() + 2
        while not received and asyncio.get_running_loop().time() < deadline:
            await asyncio.sleep(0.05)

        await client.disconnect()
        return received

    received = asyncio.run(scenario())
    assert received
    assert received[0]['id'] == 'calc-123'
    assert received[0]['status'] == 'completed'
    assert 'updatedAt' in received[0]
```

- [ ] **Step 3: Run ASGI smoke tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_socketio_asgi.py -v
```

Expected: PASS.

- [ ] **Step 4: Run combined websocket regression set**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_websocket_handlers.py tests/integration/test_socketio_asgi.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/python/tests/conftest.py src/python/tests/integration/test_socketio_asgi.py
git commit -m "test: add Socket.IO ASGI smoke coverage"
```

---

### Task 10: Electron And Gunicorn ASGI Startup

**Files:**
- Modify: `config/server-config.json`
- Modify: `src/main/config.ts`
- Modify: `src/main/python-server.ts`
- Modify: `src/main/python-env.ts`
- Modify: `src/main/menu.ts`
- Modify: `package.json`

- [ ] **Step 1: Write failing config assertion in TypeScript tests or script smoke**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
node - <<'NODE'
const fs = require('fs');
const config = JSON.parse(fs.readFileSync('config/server-config.json', 'utf8'));
if (config.gunicorn.worker_class !== 'uvicorn.workers.UvicornWorker') {
  throw new Error(`worker_class must be uvicorn.workers.UvicornWorker, got ${config.gunicorn.worker_class}`);
}
if (config.gunicorn.workers !== 1) {
  throw new Error(`Socket.IO ASGI runtime requires workers=1, got ${config.gunicorn.workers}`);
}
if (config.gunicorn.preload_app !== false) {
  throw new Error('ASGI backend must not preload process-manager/socket bindings');
}
console.log('Gunicorn ASGI config verified');
NODE
```

Expected: FAIL while `worker_class` remains the WSGI worker.

- [ ] **Step 2: Update `config/server-config.json`**

Set:

```json
"gunicorn": {
  "workers": 1,
  "threads": 4,
  "worker_class": "uvicorn.workers.UvicornWorker",
  "timeout": 0,
  "keep_alive": 2,
  "access_logfile": "-",
  "log_level": "info",
  "preload_app": false
}
```

Keep `workers` at `1` because in-memory Socket.IO rooms, notification callbacks, and file-watcher state are single-process. Keep `preload_app` false so process-manager and Socket.IO callback binding happen in the worker process that owns the ASGI event loop.

- [ ] **Step 3: Update `src/main/config.ts`**

Keep the existing `gunicorn` settings shape and ensure `worker_class` accepts the ASGI worker string:

```typescript
gunicorn: {
  workers: number;
  threads: number;
  worker_class: string;
  timeout: number;
  keep_alive: number;
  access_logfile: string | null;
  log_level: string;
  preload_app: boolean;
};
```

- [ ] **Step 4: Update `src/main/python-server.ts` Gunicorn args**

Keep `app:app` and build args from config:

```typescript
const serverSettings = serverConfig.server;
const args = [
  '-m',
  'gunicorn',
  '--workers',
  String(gunicornSettings.workers),
  '--worker-class',
  gunicornSettings.worker_class,
  '--bind',
  `${serverSettings.host}:${port}`,
  '--timeout',
  String(gunicornSettings.timeout),
  '--keep-alive',
  String(gunicornSettings.keep_alive),
  '--log-level',
  gunicornSettings.log_level,
  ...(gunicornSettings.access_logfile !== null && gunicornSettings.access_logfile !== undefined
    ? ['--access-logfile', gunicornSettings.access_logfile]
    : []),
  ...(gunicornSettings.preload_app === true ? ['--preload'] : []),
  'app:app',
];
```

Remove Flask-specific diagnostics and replace text with Python/FastAPI backend wording:

```typescript
const backendName = 'Python/FastAPI backend';
```

Update troubleshooting text:

```text
Test the backend manually: cd src/python && python app.py
```

Do not rename `flaskPort` in this task unless TypeScript references can be changed in one contained commit.

- [ ] **Step 5: Update menu/env text**

In `src/main/menu.ts`, change:

```text
Backend: Python, Flask, PySCF, RDKit
```

to:

```text
Backend: Python, FastAPI, PySCF, RDKit
```

In `src/main/python-env.ts`, update comments from Flask server to Python backend server while keeping env var names stable.

- [ ] **Step 6: Verify TypeScript**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
npm run typecheck
```

Expected: PASS.

- [ ] **Step 7: Verify Gunicorn smoke**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
npm run test:gunicorn-local
```

Expected: Gunicorn starts with `uvicorn.workers.UvicornWorker` and imports `app:app`.

- [ ] **Step 8: Commit**

```bash
git add config/server-config.json src/main/config.ts src/main/python-server.ts src/main/python-env.ts src/main/menu.ts package.json
git commit -m "feat: start backend with Gunicorn Uvicorn worker"
```

---

### Task 11: OpenAPI Contract, Codegen, And Flask Removal

**Files:**
- Modify: `src/python/tests/integration/test_api_endpoints/test_openapi_contract.py`
- Modify: `src/api-spec/openapi.yaml` only if implementation exposes a deliberate API contract change
- Regenerate: `src/python/generated_models.py` only through `npm run codegen`
- Regenerate: `src/web/types/generated-api.ts` only through `npm run codegen`
- Modify: docs/comments with Flask wording under `src/python/tests/README.md` and related test docs if they are retained
- Modify: `.github/environment.yml` if any Flask dependency remains from Task 1 until now

- [ ] **Step 1: Convert the implementation contract extractor to FastAPI**

Keep the current OpenAPI extraction, route equality test, query-parameter equality test, and method-default tests. Replace only the Flask implementation extraction helpers with FastAPI-aware parsing:

```python
_FASTAPI_METHODS = {"get", "post", "put", "patch", "delete"}
_FASTAPI_PATH_PARAM_PATTERN = re.compile(r"\{[^}]+\}")


def _normalize_impl_path(path: str) -> str:
    path = _IMPL_PATH_PARAM_PATTERN.sub("{}", path)
    return _FASTAPI_PATH_PARAM_PATTERN.sub("{}", path)


def _literal_path(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _router_prefix(tree: ast.Module) -> str:
    for item in ast.walk(tree):
        if not isinstance(item, ast.Assign):
            continue
        if not isinstance(item.value, ast.Call):
            continue
        func = item.value.func
        if isinstance(func, ast.Name) and func.id == "APIRouter":
            for keyword in item.value.keywords:
                if keyword.arg == "prefix":
                    value = _literal_path(keyword.value)
                    if value is not None:
                        return value
    return ""


def _parse_route_decorator(decorator: ast.expr) -> tuple[str, set[str]] | None:
    if not isinstance(decorator, ast.Call):
        return None
    if not isinstance(decorator.func, ast.Attribute):
        return None
    method = decorator.func.attr.lower()
    if method not in _FASTAPI_METHODS:
        return None
    if not decorator.args:
        return None
    path = _literal_path(decorator.args[0])
    if path is None:
        return None
    return path, {method.upper()}
```

Replace `_RequestArgsVisitor` with a function that reads FastAPI function signatures and excludes path params:

```python
def _path_param_names(path: str) -> set[str]:
    return {
        match.strip("{}").split(":")[-1]
        for match in _FASTAPI_PATH_PARAM_PATTERN.findall(path)
    }


def _query_params_for_function(node: ast.FunctionDef | ast.AsyncFunctionDef, route_path: str) -> set[str]:
    path_params = _path_param_names(route_path)
    args = node.args.args
    defaults = [None] * (len(args) - len(node.args.defaults)) + list(node.args.defaults)
    query_params: set[str] = set()

    def query_alias_from_call(call: ast.Call) -> str | None:
        if not isinstance(call.func, ast.Name) or call.func.id != "Query":
            return None
        for keyword in call.keywords:
            if keyword.arg == "alias":
                alias = _literal_path(keyword.value)
                if alias is not None:
                    return alias
        return ""

    def query_alias_from_annotation(annotation: ast.AST | None) -> str | None:
        if not isinstance(annotation, ast.Subscript):
            return None
        elts = annotation.slice.elts if isinstance(annotation.slice, ast.Tuple) else [annotation.slice]
        for elt in elts:
            if isinstance(elt, ast.Call):
                alias = query_alias_from_call(elt)
                if alias is not None:
                    return alias
        return None

    for arg, default in zip(args, defaults):
        if arg.arg in {"request", "body"} or arg.arg in path_params:
            continue
        default_alias = query_alias_from_call(default) if isinstance(default, ast.Call) else None
        annotation_alias = query_alias_from_annotation(arg.annotation)
        if default_alias is not None:
            query_params.add(default_alias or arg.arg)
        elif annotation_alias is not None:
            query_params.add(annotation_alias or arg.arg)
        elif default is not None or arg.annotation is not None:
            query_params.add(arg.arg)

    return query_params
```

Update `_extract_implementation_contract()` so it preserves both route-direction checks and query checks:

```python
@lru_cache(maxsize=1)
def _extract_implementation_contract() -> tuple[set[tuple[str, str]], dict[tuple[str, str], set[str]]]:
    routes: set[tuple[str, str]] = set()
    query_params_by_route: dict[tuple[str, str], set[str]] = {}

    for api_file in sorted(API_DIR.glob("*.py")):
        if api_file.name == "__init__.py":
            continue

        tree = ast.parse(api_file.read_text(encoding="utf-8"), filename=str(api_file))
        prefix = _router_prefix(tree)

        for node in tree.body:
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue

            route_specs = []
            for decorator in node.decorator_list:
                route_spec = _parse_route_decorator(decorator)
                if route_spec:
                    route_specs.append(route_spec)

            for raw_path, methods in route_specs:
                full_path = f"{prefix}{raw_path}"
                if not _is_public_contract_path(full_path):
                    continue
                normalized_path = _normalize_impl_path(full_path)
                query_names = _query_params_for_function(node, raw_path)
                for method in methods:
                    key = (method, normalized_path)
                    routes.add(key)
                    query_params_by_route.setdefault(key, set()).update(query_names)

    return routes, query_params_by_route
```

- [ ] **Step 2: Run contract tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_openapi_contract.py -v
```

Expected: PASS when route methods, paths, and query parameters match `src/api-spec/openapi.yaml`. If it fails, update implementation or the OpenAPI source in the same task; do not weaken the contract test.

- [ ] **Step 3: Run codegen**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
npm run codegen
```

Expected: generated Python and TypeScript API files are unchanged unless `src/api-spec/openapi.yaml` changed.

- [ ] **Step 4: Remove remaining production Flask imports**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
rg -n "from flask|import flask|flask(=|-cors|-pydantic|-socketio|-sock)|pytest-flask|flask_socketio|flask_pydantic|werkzeug|Flask-SocketIO|configure_flask_app" src/python scripts src/main package.json .github/environment.yml
```

Expected after cleanup: no production imports or dependency entries. `flaskPort` / `window.flaskPort` may remain in frontend/main process compatibility code because the design explicitly keeps that global name for now; do not count those hits as runtime Flask remnants. Historical test documentation may contain old explanatory text; update retained docs so current instructions say FastAPI.

Delete `configure_flask_app` from `src/python/config.py` when:

```bash
rg -n "configure_flask_app" src/python
```

returns only the function definition.

- [ ] **Step 5: Run contract and codegen checks**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_openapi_contract.py -v
```

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
npm run codegen
git diff -- src/python/generated_models.py src/web/types/generated-api.ts
```

Expected: contract test PASS. Generated-file diff is empty unless the OpenAPI contract was intentionally updated.

- [ ] **Step 6: Commit**

```bash
git add src/python/tests/integration/test_api_endpoints/test_openapi_contract.py src/api-spec/openapi.yaml src/python/generated_models.py src/web/types/generated-api.ts src/python/config.py src/python/tests/README.md .github/environment.yml
git commit -m "chore: remove Flask contract assumptions"
```

---

### Task 12: Full Verification And Final Review Loop

**Files:**
- Review only; edit the file named by any verification failure.

- [ ] **Step 1: Run full targeted backend API tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints -v
```

Expected: PASS.

- [ ] **Step 2: Run auth tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_auth_security.py tests/integration/test_auth_production.py -v
```

Expected: PASS.

- [ ] **Step 3: Run websocket tests**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_websocket_handlers.py tests/integration/test_socketio_asgi.py -v
```

Expected: PASS.

- [ ] **Step 4: Run workflow tests affected by websocket and process manager callbacks**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi/src/python
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_calculation_workflow.py tests/integration/test_pause_resume_workflow.py -v
```

Expected: PASS.

- [ ] **Step 5: Run frontend and generated-type checks**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
npm run codegen
npm run typecheck
```

Expected: both PASS.

- [ ] **Step 6: Run startup and packaging verification**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
npm run test:gunicorn-local
npm run test:python-standalone
npm run verify-env
npm run verify-build-env
```

Expected: all PASS. For `npm run test:gunicorn-local`, stop the process after the Gunicorn/Uvicorn worker imports `app:app` and reaches startup if the command is intentionally long-running.

- [ ] **Step 7: Search for Flask runtime remnants**

Run:

```bash
cd /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi
rg -n "from flask|import flask|flask(=|-cors|-pydantic|-socketio|-sock)|pytest-flask|flask_socketio|flask_pydantic|Flask-SocketIO|Werkzeug|werkzeug|configure_flask_app" src/python scripts src/main package.json .github/environment.yml config
```

Expected: no runtime references. `flaskPort` / `window.flaskPort` compatibility names may remain; if documentation references Flask as historical context, update it to describe the current FastAPI backend.

- [ ] **Step 8: Request subagent code review and loop until Important+ is zero**

Use `superpowers:requesting-code-review` and ask the reviewer to classify findings as `Critical`, `Important`, `Minor`, or `Nit`. The review prompt must include:

```text
Review the FastAPI migration branch in /Users/goodapple/workspace/PySCF_native_app/.worktrees/migration-fastapi.
Focus on regressions relative to docs/superpowers/specs/2026-05-21-fastapi-migration-design.md:
- no Flask runtime remains
- FastAPI validation errors return HTTP 400 envelopes, never 422
- /api/quantum/calculate preserves raw-body validation order
- Socket.IO ASGI CORS/auth/room/update behavior is correct
- sync callbacks schedule async emits through the shared event-loop bridge
- Electron Gunicorn startup uses uvicorn.workers.UvicornWorker with workers=1
- API paths and response envelopes match src/api-spec/openapi.yaml
Return findings ordered by severity with file/line references.
```

For each `Critical` or `Important` finding, use `superpowers:receiving-code-review`, reproduce the issue with a targeted test or command, implement the smallest fix, rerun the targeted test, and request review again with the latest diff. Stop the review loop only when the reviewer reports no `Critical` or `Important` findings.

- [ ] **Step 9: Final commit**

After the review loop is clean and verification passes:

```bash
git status --short
git add src/python .github/environment.yml package.json scripts config src/main src/api-spec docs/superpowers
git commit -m "feat: migrate backend to FastAPI"
```

Before this commit, confirm `git status --short` contains only migration files from this plan. If all prior tasks were committed separately and the worktree is clean, do not create an empty final commit.

---

## Final Handoff Checklist

- [ ] `src/python/app.py` exposes `fastapi_app`, `sio`, and `app`.
- [ ] Default FastAPI docs routes `/docs`, `/redoc`, and `/openapi.json` are disabled.
- [ ] Existing `/api-docs` and `/api-docs/spec.json` work in development mode.
- [ ] HTTP auth behavior matches development, test, and production expectations.
- [ ] Validation failures return `400` envelopes instead of FastAPI `422`.
- [ ] `/api/quantum/calculate` rejects missing `calculation_method` before Pydantic model validation.
- [ ] GPU4PySCF install rejects non-loopback clients with the existing `403` envelope.
- [ ] `/api/agent/chat` preserves `data: {"type": "chunk"}\n\n` style SSE formatting.
- [ ] Socket.IO accepts token auth, Electron-like origins, room joins, and `calculation_update` delivery.
- [ ] `NotificationService` schedules async emits through `websocket.event_loop_bridge`.
- [ ] Electron starts Gunicorn with `uvicorn.workers.UvicornWorker`.
- [ ] `workers` remains `1` in `config/server-config.json`.
- [ ] `rg "from flask|import flask|flask(=|-cors|-pydantic|-socketio|-sock)|pytest-flask|flask_socketio|flask_pydantic|configure_flask_app"` has no runtime hits.
- [ ] `npm run codegen` and `npm run typecheck` pass.
- [ ] Targeted pytest suites pass.
- [ ] `npm run test:gunicorn-local`, `npm run test:python-standalone`, `npm run verify-env`, and `npm run verify-build-env` pass.
