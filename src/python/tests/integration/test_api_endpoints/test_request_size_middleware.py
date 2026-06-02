"""
Integration and unit tests for the request-body-size middleware.

Tests the ``register_request_size_middleware`` / ``enforce_request_body_size``
behaviour defined in ``app.py``, covering:
  - Normal requests that pass through the middleware.
  - Oversized Content-Length (413).
  - Missing Content-Length on body-bearing methods (411).
  - Invalid (non-integer) Content-Length (400).
  - Non-body methods (GET) passing through unconditionally.

Cases 3 (missing Content-Length) and 4 (invalid Content-Length) are exercised
as focused unit tests against the middleware dispatch function because
httpx's TestClient automatically computes and sets the Content-Length header,
making it impossible to omit or corrupt via the normal request API.
"""

import asyncio
import json

import pytest
from fastapi import FastAPI
from starlette.requests import Request
from starlette.responses import Response

from app import MAX_REQUEST_BODY_BYTES, register_request_size_middleware


# ============================================================================
# Helpers for unit-testing the middleware dispatch function directly
# ============================================================================


def _extract_middleware_dispatch() -> callable:
    """Build a minimal FastAPI app and extract the middleware dispatch callable.

    ``enforce_request_body_size`` is a closure created inside
    ``register_request_size_middleware``.  We register the middleware on a
    throw-away app and pull the dispatch function from the middleware stack
    so that we can invoke it with hand-crafted Request objects.
    """
    probe_app = FastAPI()
    register_request_size_middleware(probe_app)
    # Starlette stores user middleware as Middleware objects; after
    # registration via @app.middleware("http") the dispatch is stored
    # in the kwargs under "dispatch".
    for middleware in probe_app.user_middleware:
        kwargs = getattr(middleware, "kwargs", {})
        if "dispatch" in kwargs:
            return kwargs["dispatch"]
    raise RuntimeError("Could not locate middleware dispatch function")


def _make_scope(
    method: str = "POST",
    path: str = "/api/quantum/calculate",
    headers: list[tuple[bytes, bytes]] | None = None,
) -> dict:
    """Return a minimal ASGI scope dict suitable for constructing a Request."""
    if headers is None:
        headers = []
    return {
        "type": "http",
        "method": method,
        "path": path,
        "query_string": b"",
        "root_path": "",
        "scheme": "http",
        "server": ("127.0.0.1", 8000),
        "headers": headers,
    }


async def _passthrough_call_next(request: Request) -> Response:
    """Dummy call_next that signals the middleware let the request through."""
    return Response(
        content=json.dumps({"passed": True}),
        status_code=200,
        media_type="application/json",
    )


# ============================================================================
# Integration tests (via the full TestClient stack)
# ============================================================================


class TestRequestSizeMiddlewareIntegration:
    """Tests that exercise the middleware through the real ASGI stack."""

    def test_post_within_limit_passes_middleware(
        self, client, mocker, valid_dft_params
    ):
        """POST with Content-Length within the limit is forwarded past the middleware.

        We mock the downstream service so that the test isolates middleware
        behaviour rather than quantum calculation logic.
        """
        # Arrange
        mock_calc_instance = {
            "id": "calc-mw-test",
            "name": "Test H2 DFT",
            "status": "pending",
            "createdAt": "2024-01-01T00:00:00",
            "parameters": valid_dft_params,
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.start_calculation.return_value = mock_calc_instance

        # Act
        response = client.post("/api/quantum/calculate", json=valid_dft_params)

        # Assert -- the request reaches the handler (202) and is NOT blocked
        # by the middleware (which would return 411/413/400).
        assert response.status_code == 202
        data = response.json()
        assert data["success"] is True
        mock_service.return_value.start_calculation.assert_called_once()

    def test_post_exceeding_limit_returns_413(self, client):
        """POST with Content-Length > MAX_REQUEST_BODY_BYTES returns 413
        with the correct JSON error envelope.
        """
        # Arrange -- create a body that exceeds the limit.
        oversized_body = b"x" * (MAX_REQUEST_BODY_BYTES + 1)

        # Act
        response = client.post(
            "/api/quantum/calculate",
            content=oversized_body,
            headers={"Content-Type": "application/json"},
        )

        # Assert
        assert response.status_code == 413
        data = response.json()
        assert data["success"] is False
        expected_error = (
            f"Request body too large. "
            f"Maximum size is {MAX_REQUEST_BODY_BYTES} bytes."
        )
        assert data["error"] == expected_error

    def test_get_without_content_length_passes_through(self, client):
        """GET requests (not in _BODY_METHODS) pass through the middleware
        even without a Content-Length header.
        """
        # Act -- /health is a simple GET endpoint that requires no auth.
        response = client.get("/health")

        # Assert -- middleware does not interfere.
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "ok"


# ============================================================================
# Unit tests (calling the middleware dispatch function directly)
# ============================================================================


class TestRequestSizeMiddlewareUnit:
    """Tests that call the middleware dispatch function directly.

    This approach is necessary for cases where httpx's TestClient
    automatically manages the Content-Length header, making it impossible
    to omit or set an invalid value through the normal request API.
    """

    @pytest.fixture(scope="class")
    def dispatch(self):
        """Extract the middleware dispatch function once per test class."""
        return _extract_middleware_dispatch()

    # ---- Helper to dispatch and parse response --------------------------------

    def _dispatch_and_parse(
        self,
        dispatch,
        method: str = "POST",
        headers: list[tuple[bytes, bytes]] | None = None,
    ) -> tuple[int, dict]:
        """Build a Request from a scope, dispatch it, and return (status, body)."""
        scope = _make_scope(method=method, headers=headers or [])
        request = Request(scope)
        response = asyncio.run(dispatch(request, _passthrough_call_next))
        return response.status_code, json.loads(response.body)

    # ---- Case 3: Missing Content-Length on body-bearing method → 411 -------

    @pytest.mark.parametrize(
        "method",
        ["POST", "PUT", "PATCH"],
        ids=["post", "put", "patch"],
    )
    def test_missing_content_length_returns_411(self, dispatch, method):
        """Body-bearing method without Content-Length returns 411."""
        status, body = self._dispatch_and_parse(dispatch, method=method, headers=[])

        assert status == 411
        assert body["success"] is False
        assert (
            body["error"] == f"Content-Length header is required for {method} requests."
        )

    # ---- Case 4: Invalid (non-integer) Content-Length → 400 ----------------

    @pytest.mark.parametrize(
        ("method", "content_length_value"),
        [
            ("POST", b"not-a-number"),
            ("POST", b""),
            # GET with garbage Content-Length is also rejected because the
            # middleware checks validity before checking the method.
            ("GET", b"garbage"),
        ],
        ids=["post-not-a-number", "post-empty", "get-garbage"],
    )
    def test_invalid_content_length_returns_400(
        self, dispatch, method, content_length_value
    ):
        """Request with non-integer Content-Length returns 400."""
        status, body = self._dispatch_and_parse(
            dispatch,
            method=method,
            headers=[(b"content-length", content_length_value)],
        )

        assert status == 400
        assert body["success"] is False
        assert body["error"] == "Invalid Content-Length header"

    # ---- Case 5 (unit): GET with no Content-Length passes through ----------

    def test_get_no_content_length_passes_through(self, dispatch):
        """GET without Content-Length passes through the middleware."""
        status, body = self._dispatch_and_parse(dispatch, method="GET", headers=[])

        assert status == 200
        assert body["passed"] is True

    # ---- Additional edge cases ---------------------------------------------

    def test_post_valid_content_length_passes_through(self, dispatch):
        """POST with a valid, within-limit Content-Length passes through."""
        status, body = self._dispatch_and_parse(
            dispatch,
            method="POST",
            headers=[(b"content-length", b"1024")],
        )

        assert status == 200
        assert body["passed"] is True

    def test_post_exactly_at_limit_passes_through(self, dispatch):
        """POST with Content-Length exactly equal to the limit passes."""
        status, body = self._dispatch_and_parse(
            dispatch,
            method="POST",
            headers=[(b"content-length", str(MAX_REQUEST_BODY_BYTES).encode())],
        )

        assert status == 200
        assert body["passed"] is True

    def test_post_one_byte_over_limit_returns_413(self, dispatch):
        """POST with Content-Length one byte over the limit returns 413."""
        over_limit = MAX_REQUEST_BODY_BYTES + 1
        status, body = self._dispatch_and_parse(
            dispatch,
            method="POST",
            headers=[(b"content-length", str(over_limit).encode())],
        )

        assert status == 413
        assert body["success"] is False
        expected_error = (
            f"Request body too large. "
            f"Maximum size is {MAX_REQUEST_BODY_BYTES} bytes."
        )
        assert body["error"] == expected_error
