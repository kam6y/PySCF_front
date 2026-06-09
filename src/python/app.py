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
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]

for _var in _THREAD_CONTROL_VARS:
    if _var not in os.environ:
        os.environ[_var] = "1"
# =============================================================================

import asyncio
import json
import logging
from collections.abc import Awaitable, Callable
from contextlib import asynccontextmanager
from typing import Any

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import ValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.responses import Response

# Maximum allowed request body size in bytes (5 MB).
# This must be larger than the largest legitimate field (e.g. xyz data for big
# molecules can be ~1 MB, ketcher JSON similarly large).  5 MB provides ample
# headroom while still preventing abuse.
MAX_REQUEST_BODY_BYTES: int = 5 * 1024 * 1024

# ServiceError status codes whose messages are curated, developer-authored,
# and safe to return to the client (e.g. for user-facing toast notifications).
# True 500 / unknown 5xx may contain interpolated exception text and are
# redacted to a generic message.
CURATED_5XX_CODES: frozenset[int] = frozenset({503, 507})

# Loopback addresses recognised as safe for TrustedHostMiddleware.
# IPv6 loopback (::1) is intentionally excluded: Starlette's
# TrustedHostMiddleware splits the Host header on the first colon
# (``host.split(":")[0]``), so a bracketed IPv6 Host like ``[::1]:5000``
# becomes ``"["`` and bare ``::1`` becomes ``""`` — neither can match the
# literal string ``"::1"``.  This app binds to 127.0.0.1; IPv6 loopback
# is not used.
_LOOPBACK_HOSTS: frozenset[str] = frozenset({"127.0.0.1", "localhost"})

# HTTP methods that may carry a request body.  For these methods the
# middleware requires a valid Content-Length header so that body-size
# enforcement cannot be bypassed via chunked Transfer-Encoding (which
# omits Content-Length).
_BODY_METHODS: frozenset[str] = frozenset({"POST", "PUT", "PATCH"})

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
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

logging.basicConfig(level=log_level, format=log_format)
logger = logging.getLogger(__name__)


def _is_development_api_docs_path(path: str) -> bool:
    return path == "/api-docs" or path.startswith("/api-docs/")


def register_request_size_middleware(fastapi_app: FastAPI) -> None:
    """Reject requests whose body size is not properly declared or exceeds the limit.

    For body-bearing methods (POST/PUT/PATCH):
      - Missing Content-Length header -> HTTP 411 (Length Required).
        This closes the chunked Transfer-Encoding bypass where a client
        omits Content-Length to evade the size check.
      - Invalid (non-integer) Content-Length -> HTTP 400.
      - Content-Length exceeds MAX_REQUEST_BODY_BYTES -> HTTP 413.

    Other methods (GET/HEAD/DELETE/OPTIONS) pass through unconditionally
    because they typically carry no body.

    The check is header-only so it never buffers or consumes the request
    body, which keeps SSE/streaming response endpoints and WebSocket
    upgrades working.  (WebSocket upgrades use scope type "websocket" and
    are not intercepted by Starlette's HTTP middleware at all.)
    """

    @fastapi_app.middleware("http")
    async def enforce_request_body_size(
        request: Request,
        call_next: Callable[[Request], Awaitable[Response]],
    ) -> Response:
        content_length_header = request.headers.get("content-length")
        if content_length_header is not None:
            try:
                content_length = int(content_length_header)
            except (ValueError, TypeError):
                content_length = -1  # force the invalid-header path below
            if content_length < 0:
                return JSONResponse(
                    {"success": False, "error": "Invalid Content-Length header"},
                    status_code=400,
                )
            if content_length > MAX_REQUEST_BODY_BYTES:
                logger.warning(
                    "Request body too large on %s: %d bytes (limit %d)",
                    request.url.path,
                    content_length,
                    MAX_REQUEST_BODY_BYTES,
                )
                return JSONResponse(
                    {
                        "success": False,
                        "error": (
                            f"Request body too large. "
                            f"Maximum size is {MAX_REQUEST_BODY_BYTES} bytes."
                        ),
                    },
                    status_code=413,
                )
        elif request.method in _BODY_METHODS:
            # Body-bearing method without Content-Length header.  Reject to
            # prevent chunked-transfer bypass of the size limit.
            logger.warning(
                "Missing Content-Length on %s %s",
                request.method,
                request.url.path,
            )
            return JSONResponse(
                {
                    "success": False,
                    "error": (
                        "Content-Length header is required for "
                        f"{request.method} requests."
                    ),
                },
                status_code=411,
            )
        return await call_next(request)


def register_auth_middleware(fastapi_app: FastAPI) -> None:
    @fastapi_app.middleware("http")
    async def verify_auth_token(
        request: Request,
        call_next: Callable[[Request], Awaitable[Response]],
    ) -> Response:
        env = os.getenv("PYSCF_ENV", "").lower()
        if _is_development_api_docs_path(request.url.path) and env == "development":
            return await call_next(request)

        if request.method == "OPTIONS":
            return await call_next(request)

        auth_token = os.getenv("PYSCF_AUTH_TOKEN")
        if auth_token:
            client_token = request.headers.get("X-Auth-Token")
            if client_token != auth_token:
                client_host = request.client.host if request.client else "unknown"
                logger.warning("Unauthorized access attempt from %s", client_host)
                return JSONResponse(
                    {"success": False, "error": "Unauthorized"}, status_code=401
                )
        else:
            # Narrowly-scoped TESTING bypass: allow requests only when the
            # application has been explicitly constructed with TESTING=True
            # (pytest fixtures) and we are NOT in production.
            is_testing = bool(getattr(request.app.state, "TESTING", False))
            if is_testing and env != "production":
                return await call_next(request)
            logger.warning("Unauthorized access attempt: Missing authentication token")
            return JSONResponse(
                {
                    "success": False,
                    "error": "Unauthorized: Missing authentication token",
                },
                status_code=401,
            )

        return await call_next(request)


def register_cors_middleware(fastapi_app: FastAPI) -> None:
    env = os.getenv("PYSCF_ENV", "").lower()

    if env in {"development", "test"}:
        # Development / test: only loopback HTTP origins with an explicit
        # port are accepted, matching the Electron-side
        # ``isAllowedDevRendererUrl`` policy which requires a non-default
        # port (portless / :80 URLs are rejected because dev servers
        # always bind to a non-default port).  Bare ``http://localhost``
        # and ``http://127.0.0.1`` (portless) are intentionally omitted.
        origins: list[str] = []
        # Port group limited to valid TCP range 1-65535.
        _VALID_PORT = (
            r"[1-9]|[1-9]\d|[1-9]\d{2}|[1-9]\d{3}"
            r"|[1-5]\d{4}|6[0-4]\d{3}|65[0-4]\d{2}|655[0-2]\d|6553[0-5]"
        )
        origin_regex: str | None = rf"^http://(127\.0\.0\.1|localhost):({_VALID_PORT})$"
    else:
        # Production / unknown: Packaged Electron builds load the renderer
        # via ``file://``, which causes the browser to send ``Origin: null``
        # (or ``file://``).  Only those origins are permitted.
        # CORS is NOT the authentication boundary — the mandatory
        # ``X-Auth-Token`` custom header is.  A cross-origin attacker cannot
        # read or forge this header (browsers block cross-origin custom
        # headers unless the preflight succeeds).  ``allow_credentials`` is
        # False so no cookies are ever reflected.
        origins = ["file://", "null"]
        origin_regex = None

    fastapi_app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_origin_regex=origin_regex,
        allow_credentials=False,
        allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"],
        allow_headers=["Cache-Control", "Content-Type", "X-Auth-Token"],
    )


def _format_validation_errors(errors: list[dict[str, Any]]) -> str:
    messages: list[str] = []
    for err in errors:
        loc = ".".join(str(part) for part in err.get("loc", []))
        msg = err.get("msg", "Invalid value")
        messages.append(f"{loc}: {msg}" if loc else msg)
    return "Validation failed: " + "; ".join(messages)


def register_exception_handlers(fastapi_app: FastAPI) -> None:
    @fastapi_app.exception_handler(json.JSONDecodeError)
    async def json_decode_error_handler(
        request: Request, error: json.JSONDecodeError
    ) -> JSONResponse:
        message = f"Validation failed: malformed JSON body: {error.msg}"
        logger.warning("Malformed JSON on %s: %s", request.url.path, message)
        return JSONResponse({"success": False, "error": message}, status_code=400)

    @fastapi_app.exception_handler(RequestValidationError)
    async def request_validation_error_handler(
        request: Request, error: RequestValidationError
    ) -> JSONResponse:
        message = _format_validation_errors(error.errors())
        logger.warning("Validation error on %s: %s", request.url.path, message)
        return JSONResponse({"success": False, "error": message}, status_code=400)

    @fastapi_app.exception_handler(ValidationError)
    async def pydantic_validation_error_handler(
        request: Request, error: ValidationError
    ) -> JSONResponse:
        message = _format_validation_errors(error.errors())
        logger.warning("Pydantic validation error on %s: %s", request.url.path, message)
        return JSONResponse({"success": False, "error": message}, status_code=400)

    @fastapi_app.exception_handler(StarletteHTTPException)
    async def http_exception_handler(
        request: Request, error: StarletteHTTPException
    ) -> JSONResponse:
        if error.status_code == 404:
            return JSONResponse(
                {"success": False, "error": "Not Found"}, status_code=404
            )
        if error.status_code == 400:
            return JSONResponse(
                {"success": False, "error": f"Validation failed: {error.detail}"},
                status_code=400,
            )
        if error.status_code >= 500:
            # Application code uses ServiceError / ResourceUnavailableError for
            # curated 5xx responses; StarletteHTTPException 5xx are
            # framework-internal and always redacted to prevent leaking
            # raw detail strings.
            logger.error(
                "HTTP %d on %s: %s",
                error.status_code,
                request.url.path,
                error.detail,
                exc_info=error,
            )
            return JSONResponse(
                {"success": False, "error": "An internal server error occurred."},
                status_code=error.status_code,
            )
        return JSONResponse(
            {"success": False, "error": str(error.detail)},
            status_code=error.status_code,
        )

    @fastapi_app.exception_handler(ServiceError)
    async def service_error_handler(
        request: Request, error: ServiceError
    ) -> JSONResponse:
        if error.status_code >= 500:
            logger.error(
                "Service error on %s: %s",
                request.url.path,
                error.message,
                exc_info=error,
            )
            # Curated operational messages (503/507) are safe to surface;
            # all other 5xx are redacted to prevent leaking raw exception
            # text or internal paths.
            client_message = (
                error.message
                if error.status_code in CURATED_5XX_CODES
                else "An internal server error occurred."
            )
            return JSONResponse(
                {"success": False, "error": client_message},
                status_code=error.status_code,
            )
        return JSONResponse(
            {"success": False, "error": error.message},
            status_code=error.status_code,
        )

    @fastapi_app.exception_handler(Exception)
    async def unhandled_exception_handler(
        request: Request, error: Exception
    ) -> JSONResponse:
        logger.error(
            "Unhandled exception on %s: %s", request.url.path, error, exc_info=error
        )
        return JSONResponse(
            {"success": False, "error": "An internal server error occurred."},
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


def create_fastapi_app(
    server_port: int | None = None, test_config: dict[str, Any] | None = None
) -> FastAPI:
    server_config = get_server_config()
    if server_port is None:
        port_env = os.getenv("PYSCF_SERVER_PORT")
        server_port = determine_server_port(server_config, port_env=port_env)

    fastapi_app = FastAPI(
        docs_url=None, redoc_url=None, openapi_url=None, lifespan=lifespan
    )
    configure_fastapi_app(fastapi_app, server_config, server_port)
    # Middleware registration order matters: Starlette executes the *last*
    # registered middleware first (outermost).  Desired execution order:
    #   0. Host validation  (outermost - reject invalid Host headers first)
    #   1. CORS  (must wrap 413/401 responses with CORS headers)
    #   2. Body-size check  (reject oversized bodies before auth/routing)
    #   3. Auth  (innermost - only reached if body size is OK)
    # So register in reverse: auth -> body-size -> CORS -> host validation.
    register_auth_middleware(fastapi_app)
    register_request_size_middleware(fastapi_app)
    register_cors_middleware(fastapi_app)
    register_exception_handlers(fastapi_app)

    # Host validation: restrict to loopback addresses appropriate for this
    # local desktop app.  ``testserver`` is included only in development/test
    # because Starlette's TestClient sends ``Host: testserver`` by default;
    # all other environments (production, unset, typos) must NOT trust that
    # synthetic hostname (fail-closed allowlist).
    configured_host = server_config.get("server.host", "127.0.0.1")
    if not isinstance(configured_host, str) or not configured_host:
        logger.warning(
            "Invalid server.host configuration value %r; falling back to 127.0.0.1",
            configured_host,
        )
        configured_host = "127.0.0.1"
    elif configured_host not in _LOOPBACK_HOSTS:
        logger.warning(
            "server.host %r is not a recognised loopback address; "
            "falling back to 127.0.0.1",
            configured_host,
        )
        configured_host = "127.0.0.1"

    env = os.getenv("PYSCF_ENV", "").lower()
    allowed_hosts: list[str] = list({"127.0.0.1", "localhost", configured_host})
    if env in {"development", "test"}:
        allowed_hosts.append("testserver")
    fastapi_app.add_middleware(TrustedHostMiddleware, allowed_hosts=allowed_hosts)

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


if __name__ == "__main__":
    import uvicorn

    server_config = get_server_config()
    port_env = os.getenv("PYSCF_SERVER_PORT")
    actual_port = determine_server_port(server_config, port_env=port_env)
    host = server_config.get("server.host", "127.0.0.1")
    uvicorn.run("app:app", host=host, port=actual_port, reload=False, log_level="info")
