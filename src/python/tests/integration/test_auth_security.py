import pytest
from fastapi.testclient import TestClient

from app import create_fastapi_app


TEST_TOKEN = "pytest-test-token-12345"


def make_auth_client(monkeypatch, token=TEST_TOKEN, env="development"):
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
