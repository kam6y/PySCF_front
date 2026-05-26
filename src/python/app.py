import os

# =============================================================================
# Thread Control: Set default thread count BEFORE importing any numerical libraries
# =============================================================================
# OpenBLAS, MKL, and other BLAS/LAPACK libraries read environment variables
# only at import time. On Linux, they may ignore later changes.
# Setting these to '1' ensures that when PySCF/NumPy are imported (via quantum_calc),
# the libraries initialize with a single thread by default.
# The actual thread count is set per-calculation in the worker process.
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
# =============================================================================

import asyncio
import json
import logging
from contextlib import asynccontextmanager
from typing import Any

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import ValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException

from api import register_routers
from config import (
    ConfigurationError,
    configure_fastapi_app,
    determine_server_port,
    get_server_config,
)
from quantum_calc import shutdown_process_manager
from services.exceptions import ServiceError
from websocket.event_loop_bridge import bind_event_loop, clear_event_loop

try:
    _server_config = get_server_config()
    log_level = getattr(logging, _server_config.get_logging_level().upper())
    log_format = _server_config.get_logging_format()
except ConfigurationError as exc:
    print(f"WARNING: Configuration error: {exc}. Using default logging settings.")
    log_level = logging.INFO
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

logging.basicConfig(level=log_level, format=log_format)
logger = logging.getLogger(__name__)


def _is_development_api_docs_path(path: str) -> bool:
    return path == '/api-docs' or path.startswith('/api-docs/')


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
        allow_headers=['Cache-Control', 'Content-Type', 'X-Auth-Token'],
    )


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


@asynccontextmanager
async def lifespan(fastapi_app: FastAPI):
    bind_event_loop(asyncio.get_running_loop())
    initialize_process_manager_notifications()
    try:
        yield
    finally:
        clear_event_loop()
        shutdown_process_manager(wait=False, force=True)


def initialize_process_manager_notifications() -> None:
    from quantum_calc import initialize_process_manager_with_callback
    from services.notification_service import get_notification_service

    try:
        notification_service = get_notification_service()
        initialize_process_manager_with_callback(
            notification_callback=notification_service.send_calculation_update
        )
    except Exception as exc:
        logger.error("Failed to initialize process manager with callback: %s", exc)


def create_fastapi_app(server_port: int | None = None, test_config: dict[str, Any] | None = None) -> FastAPI:
    server_config = get_server_config()
    if server_port is None:
        port_env = os.getenv('PYSCF_SERVER_PORT')
        server_port = determine_server_port(server_config, port_env=port_env)

    fastapi_app = FastAPI(docs_url=None, redoc_url=None, openapi_url=None, lifespan=lifespan)
    configure_fastapi_app(fastapi_app, server_config, server_port)
    register_auth_middleware(fastapi_app)
    register_cors_middleware(fastapi_app)
    register_exception_handlers(fastapi_app)

    if test_config:
        for key, value in test_config.items():
            setattr(fastapi_app.state, key, value)

    register_routers(fastapi_app)
    return fastapi_app


def create_app(
    server_port: int | None = None,
    test_config: dict[str, Any] | None = None,
) -> FastAPI:
    return create_fastapi_app(server_port=server_port, test_config=test_config)


app = create_app()
fastapi_app = app


if __name__ == '__main__':
    import uvicorn

    server_config = get_server_config()
    port_env = os.getenv('PYSCF_SERVER_PORT')
    actual_port = determine_server_port(server_config, port_env=port_env)
    host = server_config.get('server.host', '127.0.0.1')
    uvicorn.run('app:app', host=host, port=actual_port, reload=False, log_level='info')
