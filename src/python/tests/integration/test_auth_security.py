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
    return TestClient(app)


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
    """Test that development requests are allowed when no auth token is configured."""
    with make_auth_client(monkeypatch, token=None, env="development") as client:
        response = client.get("/health")

    assert response.status_code != 401


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


def test_cors_production_allows_null_origin(monkeypatch):
    """In production (packaged) mode, null origin is allowed for file:// renderer."""
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    with TestClient(app) as client:
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
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    with TestClient(app) as client:
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


def test_cors_production_allows_file_origin(monkeypatch):
    """In production (packaged) mode, file:// origin is allowed for Electron renderer."""
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    with TestClient(app) as client:
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


def test_cors_production_rejects_external_origin(monkeypatch):
    """In production mode, an external origin is rejected."""
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    with TestClient(app) as client:
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
    with TestClient(app) as client:
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
    with TestClient(app) as client:
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
    """Mixed-case PYSCF_ENV like 'Development' is normalized to 'development' by auth."""
    monkeypatch.setenv("PYSCF_ENV", "Development")
    monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": False})
    with TestClient(app) as client:
        response = client.get("/health")

    # With .lower() normalization, "Development" is treated as "development"
    # so the request is allowed (not 401).
    assert response.status_code != 401


def test_mixed_case_env_consistent_across_auth_and_cors(monkeypatch):
    """Mixed-case PYSCF_ENV is handled consistently by auth and CORS middleware."""
    monkeypatch.setenv("PYSCF_ENV", "Development")
    monkeypatch.delenv("PYSCF_AUTH_TOKEN", raising=False)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": False})
    with TestClient(app) as client:
        # Auth: should allow without token (development mode)
        auth_response = client.get("/health")
        # CORS: should allow loopback origin (development mode)
        cors_response = client.options(
            "/health",
            headers={
                "Origin": "http://127.0.0.1:5173",
                "Access-Control-Request-Method": "GET",
                "Access-Control-Request-Headers": "X-Auth-Token",
            },
        )

    assert auth_response.status_code != 401
    assert (
        cors_response.headers.get("access-control-allow-origin")
        == "http://127.0.0.1:5173"
    )


def test_cors_production_does_not_send_credentials_header(monkeypatch):
    """In production mode, Access-Control-Allow-Credentials is not sent."""
    monkeypatch.setenv("PYSCF_ENV", "production")
    monkeypatch.setenv("PYSCF_AUTH_TOKEN", TEST_TOKEN)
    monkeypatch.delenv("PYSCF_RESOURCES_PATH", raising=False)
    app = create_fastapi_app(server_port=5000, test_config={"TESTING": True})
    with TestClient(app) as client:
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
