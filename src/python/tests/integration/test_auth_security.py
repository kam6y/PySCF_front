import pytest
from fastapi.testclient import TestClient

from app import create_fastapi_app


TEST_TOKEN = "pytest-test-token-12345"


def make_auth_client(
    monkeypatch: pytest.MonkeyPatch,
    token: str | None = TEST_TOKEN,
    env: str = "development",
) -> TestClient:
    if token is None:
        monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    else:
        monkeypatch.setenv("PYSCF_AUTH_TOKEN", token)
    monkeypatch.setenv("PYSCF_ENV", env)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(
        server_port=5000, test_config={"TESTING": env != "production"}
    )
    # Production builds exclude ``testserver`` from TrustedHostMiddleware
    # allowed_hosts, so TestClient must send a recognised loopback Host.
    base_url = "http://127.0.0.1" if env == "production" else "http://testserver"
    return TestClient(app, base_url=base_url)


def make_production_client(
    monkeypatch: pytest.MonkeyPatch,
) -> TestClient:
    """Create a TestClient configured for production mode."""
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    return TestClient(app, base_url="http://127.0.0.1")


def test_request_without_token(monkeypatch):
    """Test that requests without a token are rejected."""
    with make_auth_client(monkeypatch) as client:
        response = client.get("/health")

    assert response.status_code == 401
    assert response.json() == {"success": False, "error": "Unauthorized"}


def test_request_with_wrong_token(monkeypatch):
    """Test that requests with an incorrect token are rejected."""
    with make_auth_client(monkeypatch) as client:
        response = client.get("/health", headers={"X-Auth-Token": "wrong-token"})

    assert response.status_code == 401
    assert response.json() == {"success": False, "error": "Unauthorized"}


def test_request_with_correct_token(monkeypatch):
    """Test that requests with the correct token are accepted."""
    with make_auth_client(monkeypatch) as client:
        response = client.get("/health", headers={"X-Auth-Token": TEST_TOKEN})

    assert response.status_code == 200


def test_options_request_cors(monkeypatch):
    """Test that OPTIONS requests are allowed without token (for CORS)."""
    with make_auth_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:3000",
                "Access-Control-Request-Method": "PATCH",
                "Access-Control-Request-Headers": "X-Auth-Token, Content-Type",
            },
        )

    assert response.status_code in {200, 204}
    assert "X-Auth-Token" in response.headers["access-control-allow-headers"]


def test_sse_options_request_accepts_event_stream_headers(monkeypatch):
    """Test that SSE preflight requests accept browser event-stream headers."""
    with make_auth_client(monkeypatch) as client:
        response = client.options(
            "/api/quantum/calculations/updates/stream",
            headers={
                "Origin": "http://localhost:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "cache-control,x-auth-token",
            },
        )

    assert response.status_code in {200, 204}
    allowed_headers = response.headers["access-control-allow-headers"]
    assert "Cache-Control" in allowed_headers
    assert "X-Auth-Token" in allowed_headers


def test_api_docs_html_without_token_in_development(monkeypatch):
    """Test that development API docs HTML is accessible without a token."""
    with make_auth_client(monkeypatch) as client:
        response = client.get("/api-docs")

    assert response.status_code == 200
    assert b"PySCF Front API Docs" in response.content


def test_api_docs_spec_without_token_in_development(monkeypatch):
    """Test that development API docs spec is accessible without a token."""
    with make_auth_client(monkeypatch) as client:
        response = client.get("/api-docs/spec.json")

    assert response.status_code == 200
    assert response.json()["openapi"]


def test_api_docs_requires_token_outside_development(monkeypatch):
    """Test that API docs auth bypass is limited to development mode."""
    with make_auth_client(monkeypatch, env="production") as client:
        response = client.get("/api-docs")

    assert response.status_code == 401
    assert response.json() == {"success": False, "error": "Unauthorized"}


def test_missing_token_in_production(monkeypatch):
    """Test that production requests fail when no auth token is configured."""
    with make_auth_client(monkeypatch, token=None, env="production") as client:
        response = client.get("/health")

    assert response.status_code == 401
    assert response.json() == {
        "success": False,
        "error": "Unauthorized: Missing authentication token",
    }


def test_missing_token_in_development(monkeypatch):
    """Test that development requests are rejected when no auth token is configured.

    Even in development mode, PYSCF_AUTH_TOKEN must be set for non-test
    runtime requests.  The TESTING bypass is not active here (TESTING=False)
    to simulate a real runtime scenario.
    """
    monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    monkeypatch.setenv("PYSCF_ENV", "development")
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": False})
    with TestClient(app) as client:
        response = client.get("/health")

    assert response.status_code == 401
    assert response.json() == {
        "success": False,
        "error": "Unauthorized: Missing authentication token",
    }


# ---------------------------------------------------------------------------
# PYSCF_ENABLE_DEV_DOCS gate tests (SEC-005 hardening)
# ---------------------------------------------------------------------------

_DOCS_DISABLED_BODY = {"success": False, "error": "Developer docs are disabled."}


@pytest.mark.parametrize(
    "path",
    ["/api-docs", "/api-docs/init.js", "/api-docs/spec.json"],
    ids=["html", "init-js", "spec-json"],
)
def test_dev_docs_blocked_when_disabled(monkeypatch, path):
    """GET dev-docs paths return 403 when PYSCF_ENABLE_DEV_DOCS=0."""
    monkeypatch.setenv("PYSCF_ENABLE_DEV_DOCS", "0")
    with make_auth_client(monkeypatch) as client:
        response = client.get(path, headers={"X-Auth-Token": TEST_TOKEN})

    assert response.status_code == 403
    assert response.json() == _DOCS_DISABLED_BODY


def _check_bytes(expected: bytes):
    """Return a content checker that asserts expected bytes in response.content."""

    def _check(response):
        assert expected in response.content

    return _check


def _check_openapi_truthy(response):
    """Assert that the OpenAPI spec field is truthy."""
    assert response.json()["openapi"]


@pytest.mark.parametrize(
    ("env_action", "path", "content_check"),
    [
        ("delenv", "/api-docs", _check_bytes(b"PySCF Front API Docs")),
        ("delenv", "/api-docs/init.js", _check_bytes(b"SwaggerUIBundle")),
        ("delenv", "/api-docs/spec.json", _check_openapi_truthy),
        ("setenv", "/api-docs", _check_bytes(b"PySCF Front API Docs")),
        ("setenv", "/api-docs/init.js", _check_bytes(b"SwaggerUIBundle")),
        ("setenv", "/api-docs/spec.json", _check_openapi_truthy),
    ],
    ids=[
        "absent-html",
        "absent-init-js",
        "absent-spec-json",
        "enabled-html",
        "enabled-init-js",
        "enabled-spec-json",
    ],
)
def test_dev_docs_served_when_allowed(monkeypatch, env_action, path, content_check):
    """GET dev-docs paths return 200 with correct content when docs are allowed."""
    if env_action == "delenv":
        monkeypatch.delenv("PYSCF_ENABLE_DEV_DOCS", raising=False)
    else:
        monkeypatch.setenv("PYSCF_ENABLE_DEV_DOCS", "1")
    with make_auth_client(monkeypatch) as client:
        response = client.get(path, headers={"X-Auth-Token": TEST_TOKEN})

    assert response.status_code == 200
    content_check(response)


# ---------------------------------------------------------------------------
# SEC-005: CORS origin conditioning on packaged vs development mode
# ---------------------------------------------------------------------------


def test_cors_development_allows_loopback_origin(monkeypatch):
    """In development mode, loopback HTTP origins are allowed."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    assert (
        response.headers.get("access-control-allow-origin") == "http://127.0.0.1:5173"
    )


def test_cors_development_allows_localhost_origin(monkeypatch):
    """In development mode, localhost HTTP origins are allowed."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://localhost:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    assert (
        response.headers.get("access-control-allow-origin") == "http://localhost:5173"
    )


def test_cors_development_rejects_null_origin(monkeypatch):
    """In development mode, null origin is not allowed."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "null",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    # CORSMiddleware will not set access-control-allow-origin for disallowed origins
    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "null"


def test_cors_production_allows_null_origin_transitional(monkeypatch):
    """In production mode, the literal Origin 'null' is ACCEPTED (transitional).

    Opaque origins (data: URIs, sandboxed iframes, local file:// pages)
    send ``Origin: null``.  This is a deliberate transition fallback while
    the app:// custom protocol migration is validated.

    REMOVAL CONDITION (see ``register_cors_middleware`` in app.py):
      Remove 'null' from ``allow_origins`` once:
        1. app:// is GUI-smoke-tested on all target platforms.
        2. Electron serializes the origin as 'app://renderer' in every
           IPC and fetch path.
        3. No production crash/error reports reference a CORS rejection
           for 'null' or 'file://' origins.
      When 'null' is removed, this test SHOULD fail -- update or delete
      it to confirm the intentional policy change.
    """
    with make_production_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "null",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    assert response.headers.get("access-control-allow-origin") == "null"


def test_cors_production_rejects_loopback_origin(monkeypatch):
    """In production (packaged) mode, loopback HTTP origins are not allowed."""
    with make_production_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "http://127.0.0.1:5173"


def test_cors_production_allows_file_origin_transitional(monkeypatch):
    """In production mode, 'file://' origin is ACCEPTED (transitional).

    Electron's file:// renderer sends ``Origin: file://``.  This is a
    deliberate transition fallback alongside 'null' while the app://
    custom protocol migration is validated.

    REMOVAL CONDITION (see ``register_cors_middleware`` in app.py):
      Same conditions as ``test_cors_production_allows_null_origin_transitional``.
      When 'file://' is removed from the origin list, this test SHOULD
      fail -- update or delete it to confirm the intentional change.
    """
    with make_production_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "file://",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    assert response.headers.get("access-control-allow-origin") == "file://"


def test_cors_production_allows_app_renderer_origin(monkeypatch):
    """In production mode, Origin 'app://renderer' is ACCEPTED.

    This is the PRIMARY origin for the packaged Electron app (NOT
    transitional).  Electron's custom ``app://`` protocol scheme causes
    the renderer to send ``Origin: app://renderer`` on every fetch to
    the local backend.  If this origin were accidentally dropped from
    ``allow_origins``, the packaged app would silently break with no
    test signal.
    """
    with make_production_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "app://renderer",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    assert response.headers.get("access-control-allow-origin") == "app://renderer"


def test_cors_production_rejects_external_origin(monkeypatch):
    """In production mode, an external origin is rejected."""
    with make_production_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://evil.com",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "http://evil.com"


def test_cors_development_rejects_external_origin(monkeypatch):
    """In development mode, an external origin is rejected."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://evil.com",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "http://evil.com"


@pytest.mark.parametrize(
    "invalid_port",
    ["0", "65536", "99999"],
    ids=["port-0", "port-65536", "port-99999"],
)
def test_cors_development_rejects_invalid_port(monkeypatch, invalid_port):
    """In development mode, loopback origins with invalid ports are rejected."""
    origin = f"http://localhost:{invalid_port}"
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": origin,
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != origin
    assert allow_origin != "*"


def test_cors_development_rejects_https_loopback(monkeypatch):
    """In development mode, https loopback origins are rejected (http only)."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "https://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "https://127.0.0.1:5173"


@pytest.mark.parametrize(
    "origin",
    ["http://localhost", "http://127.0.0.1"],
    ids=["localhost-portless", "127-portless"],
)
def test_cors_development_rejects_portless_loopback(monkeypatch, origin):
    """In development mode, portless loopback origins are rejected (parity with Electron)."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": origin,
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != origin
    assert allow_origin != "*"


# ---------------------------------------------------------------------------
# SEC-005: CORS fail-closed for unknown / unset PYSCF_ENV
# ---------------------------------------------------------------------------


def test_cors_unset_env_uses_production_origins(monkeypatch):
    """When PYSCF_ENV is unset, CORS falls through to restrictive (production) origins."""
    monkeypatch.delenv("PYSCF_ENV", raising=False)
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    # Use loopback base_url so TrustedHostMiddleware does not short-circuit
    # the request before CORS runs (testserver is not trusted when env is
    # unset under the fail-closed allowlist).
    with TestClient(app, base_url="http://127.0.0.1") as client:
        # Loopback should be rejected (production only allows file:// and null)
        response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "http://127.0.0.1:5173"


def test_cors_typo_env_uses_production_origins(monkeypatch):
    """When PYSCF_ENV is a typo like 'prod', CORS falls through to restrictive origins."""
    monkeypatch.setenv("PYSCF_ENV", "prod")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    # Use loopback base_url so TrustedHostMiddleware does not short-circuit
    # the request before CORS runs (testserver is not trusted when env is a
    # typo under the fail-closed allowlist).
    with TestClient(app, base_url="http://127.0.0.1") as client:
        # Loopback should be rejected
        response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    allow_origin = response.headers.get("access-control-allow-origin")
    assert allow_origin != "http://127.0.0.1:5173"


def test_mixed_case_env_treated_as_development_by_auth(monkeypatch):
    """Mixed-case PYSCF_ENV without a token is rejected regardless of env value.

    Even though .lower() normalizes 'Development' to 'development', auth
    now requires PYSCF_AUTH_TOKEN for all non-test runtime requests.
    """
    monkeypatch.setenv("PYSCF_ENV", "Development")
    monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": False})
    with TestClient(app) as client:
        response = client.get("/health")

    assert response.status_code == 401
    assert response.json() == {
        "success": False,
        "error": "Unauthorized: Missing authentication token",
    }


def test_mixed_case_env_consistent_across_auth_and_cors(monkeypatch):
    """Mixed-case PYSCF_ENV: auth rejects without token; CORS still uses dev origins.

    Auth now requires PYSCF_AUTH_TOKEN regardless of env, so a tokenless
    request is rejected.  CORS still normalises the env correctly and
    allows loopback origins in development mode (OPTIONS bypass auth).
    """
    monkeypatch.setenv("PYSCF_ENV", "Development")
    monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": False})
    with TestClient(app) as client:
        # Auth: should reject without token (fail-closed)
        auth_response = client.get("/health")
        # CORS: should still allow loopback origin (development mode, OPTIONS
        # requests bypass auth)
        cors_response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert auth_response.status_code == 401
    assert (
        cors_response.headers.get("access-control-allow-origin")
        == "http://127.0.0.1:5173"
    )


def test_cors_production_does_not_send_credentials_header(monkeypatch):
    """In production mode, Access-Control-Allow-Credentials is not sent."""
    with make_production_client(monkeypatch) as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "null",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    # allow_credentials=False means the header should not be 'true'
    creds_header = response.headers.get("access-control-allow-credentials")
    assert creds_header != "true"


def test_cors_development_does_not_send_credentials_header(monkeypatch):
    """In development mode, Access-Control-Allow-Credentials is not sent."""
    with make_auth_client(monkeypatch, env="development") as client:
        response = client.options(
            "/health",
            headers={
                "Origin": "http://localhost:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert response.status_code in {200, 204}
    # allow_credentials=False means the header should not be 'true'
    creds_header = response.headers.get("access-control-allow-credentials")
    assert creds_header != "true"


# ---------------------------------------------------------------------------
# TrustedHostMiddleware: Host header validation
# ---------------------------------------------------------------------------


def test_trusted_host_rejects_evil_host(monkeypatch):
    """A request with Host: evil.com is rejected by TrustedHostMiddleware."""
    with make_auth_client(monkeypatch) as client:
        response = client.get("/health", headers={"host": "evil.com"})

    assert response.status_code == 400
    assert "Invalid host header" in response.text


def test_trusted_host_accepts_loopback(monkeypatch):
    """A request with Host: 127.0.0.1 passes TrustedHostMiddleware."""
    with make_auth_client(monkeypatch) as client:
        response = client.get(
            "/health",
            headers={"host": "127.0.0.1", "X-Auth-Token": TEST_TOKEN},
        )

    assert response.status_code == 200


# ---------------------------------------------------------------------------
# TESTING bypass production boundary
# ---------------------------------------------------------------------------


def test_testing_bypass_blocked_in_production(monkeypatch):
    """TESTING=True does NOT bypass auth when PYSCF_ENV=production.

    Even when the app is constructed with test_config={"TESTING": True},
    the TESTING bypass must not activate in production mode.
    """
    monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    with TestClient(app, base_url="http://127.0.0.1") as client:
        response = client.get("/health")

    assert response.status_code == 401
    assert response.json() == {
        "success": False,
        "error": "Unauthorized: Missing authentication token",
    }


# ---------------------------------------------------------------------------
# IMP-1: testserver host rejected in production
# ---------------------------------------------------------------------------


def test_testserver_host_rejected_in_production(monkeypatch):
    """Host: testserver is rejected by TrustedHostMiddleware in production.

    Production builds exclude ``testserver`` from allowed_hosts so the
    synthetic hostname that Starlette's TestClient sends by default is
    not trusted.
    """
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    # Default TestClient base_url sends Host: testserver
    with TestClient(app) as client:
        response = client.get("/health")

    assert response.status_code == 400
    assert "Invalid host header" in response.text


# ---------------------------------------------------------------------------
# IMP-2: non-loopback configured host falls back to 127.0.0.1
# ---------------------------------------------------------------------------


def test_non_loopback_configured_host_falls_back(monkeypatch, caplog):
    """A non-loopback server.host config value is rejected with a warning.

    TrustedHostMiddleware should still accept 127.0.0.1 and the
    non-loopback value should NOT be in allowed_hosts.
    """
    import logging
    from unittest.mock import MagicMock

    monkeypatch.setenv("PYSCF_ENV", "development")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)

    # Stub server config to return a non-loopback host
    fake_config = MagicMock()
    fake_config.get.side_effect = lambda key, default=None: {
        "server.host": "0.0.0.0",
        "server.port": 5000,
    }.get(key, default)
    fake_config.get_logging_level.return_value = "INFO"
    fake_config.get_logging_format.return_value = (
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    monkeypatch.setattr("app.get_server_config", lambda: fake_config)

    with caplog.at_level(logging.WARNING, logger="app"):
        app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})

    # Warning should have been emitted about the non-loopback host
    assert any(
        "not a recognised loopback address" in record.getMessage()
        for record in caplog.records
    ), "Expected a warning about non-loopback host fallback"

    with TestClient(app) as client:
        # 127.0.0.1 should be accepted (loopback fallback)
        ok_response = client.get(
            "/health",
            headers={"host": "127.0.0.1", "X-Auth-Token": TEST_TOKEN},
        )
        assert ok_response.status_code == 200

        # The original non-loopback host should be rejected
        bad_response = client.get(
            "/health",
            headers={"host": "0.0.0.0", "X-Auth-Token": TEST_TOKEN},
        )
        assert bad_response.status_code == 400
        assert "Invalid host header" in bad_response.text


# ---------------------------------------------------------------------------
# IMP-T3: StarletteHTTPException 5xx redaction
# ---------------------------------------------------------------------------


def test_starlette_http_exception_5xx_redacts_detail():
    """A StarletteHTTPException with status 500 must redact the raw detail.

    The response body should be the generic envelope, not the raw detail.
    """
    from fastapi import FastAPI
    from starlette.exceptions import HTTPException as StarletteHTTPException

    from app import register_exception_handlers

    test_app = FastAPI()
    register_exception_handlers(test_app)

    @test_app.get("/boom")
    def boom():
        raise StarletteHTTPException(status_code=500, detail="secret db info")

    with TestClient(test_app, raise_server_exceptions=False) as client:
        response = client.get("/boom")

    assert response.status_code == 500
    body = response.json()
    assert body == {"success": False, "error": "An internal server error occurred."}
    assert "secret db info" not in response.text


# ---------------------------------------------------------------------------
# IMP-T4: ServiceError 507 curated message preserved
# ---------------------------------------------------------------------------


def test_service_error_507_preserves_curated_message():
    """A ServiceError with status 507 should surface its curated message.

    507 is in CURATED_5XX_CODES, so the developer-authored message must
    be returned to the client rather than the generic redacted text.
    """
    from fastapi import FastAPI

    from app import register_exception_handlers
    from services.exceptions import InsufficientResourcesError

    test_app = FastAPI()
    register_exception_handlers(test_app)

    @test_app.get("/disk-full")
    def disk_full():
        raise InsufficientResourcesError("Disk space exhausted")

    with TestClient(test_app, raise_server_exceptions=False) as client:
        response = client.get("/disk-full")

    assert response.status_code == 507
    body = response.json()
    assert body == {"success": False, "error": "Disk space exhausted"}


# ---------------------------------------------------------------------------
# F2: testserver fail-closed when PYSCF_ENV is unset/empty
# ---------------------------------------------------------------------------


def test_testserver_host_rejected_when_env_unset(monkeypatch):
    """Host: testserver is rejected when PYSCF_ENV is unset (fail-closed).

    The fail-closed allowlist only includes ``testserver`` when the
    lowercased env is in {"development", "test"}.  When PYSCF_ENV is
    unset (empty string after .lower()), ``testserver`` must NOT be
    trusted.
    """
    monkeypatch.delenv("PYSCF_ENV", raising=False)
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    # Default TestClient base_url sends Host: testserver
    with TestClient(app) as client:
        response = client.get("/health", headers={"X-Auth-Token": TEST_TOKEN})

    assert response.status_code == 400
    assert "Invalid host header" in response.text


# ---------------------------------------------------------------------------
# F3: Host header with port (production-realistic)
# ---------------------------------------------------------------------------


def test_trusted_host_accepts_loopback_with_port(monkeypatch):
    """A request with Host: 127.0.0.1:5000 passes TrustedHostMiddleware.

    Production runtime sends ``Host: 127.0.0.1:<port>``.  Starlette
    strips the port via ``host.split(":")[0]`` before matching, so
    ``127.0.0.1`` must be in allowed_hosts.  This test pins that
    port-stripping behaviour so a future Starlette upgrade that changes
    it is caught.
    """
    with make_auth_client(monkeypatch) as client:
        response = client.get(
            "/health",
            headers={"host": "127.0.0.1:5000", "X-Auth-Token": TEST_TOKEN},
        )

    assert response.status_code == 200


# ---------------------------------------------------------------------------
# L4: ServiceError 500 redacts raw detail
# ---------------------------------------------------------------------------


def test_service_error_500_redacts_raw_detail():
    """A ServiceError with status 500 must redact the raw detail message.

    500 is NOT in CURATED_5XX_CODES, so the response body should be the
    generic envelope and the raw detail must not leak to the client.
    """
    from fastapi import FastAPI

    from app import register_exception_handlers
    from services.exceptions import ServiceError

    test_app = FastAPI()
    register_exception_handlers(test_app)

    @test_app.get("/internal-boom")
    def internal_boom():
        raise ServiceError("secret traceback /var/db/creds.py", status_code=500)

    with TestClient(test_app, raise_server_exceptions=False) as client:
        response = client.get("/internal-boom")

    assert response.status_code == 500
    body = response.json()
    assert body == {"success": False, "error": "An internal server error occurred."}
    assert "secret traceback" not in response.text
    assert "creds.py" not in response.text
