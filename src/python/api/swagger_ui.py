"""
Development-only Swagger UI endpoint.

Security hardening (SEC-005):
- Pinned exact swagger-ui-dist version with SRI integrity hashes
- Content-Security-Policy header restricts resource loading
- Inline init script moved to a served endpoint (/api-docs/init.js)
- Only registered when not in packaged mode (see api/__init__.py)
"""

import os
from pathlib import Path

import yaml
from fastapi import APIRouter
from fastapi.responses import HTMLResponse, JSONResponse, Response

router = APIRouter()

# Pinned swagger-ui-dist version and SRI hashes (sha384).
# Update these when upgrading swagger-ui-dist.
#
# To regenerate SRI hashes for a given version, run:
#   VER=5.18.2  # replace with the target version
#   curl -sSL "https://unpkg.com/swagger-ui-dist@${VER}/swagger-ui.css" \
#     | openssl dgst -sha384 -binary | openssl base64 -A
#   curl -sSL "https://unpkg.com/swagger-ui-dist@${VER}/swagger-ui-bundle.js" \
#     | openssl dgst -sha384 -binary | openssl base64 -A
# Prefix each output with "sha384-" to form the final integrity value.
SWAGGER_UI_VERSION = "5.18.2"
SWAGGER_UI_CSS_SRI = (
    "sha384-rcbEi6xgdPk0iWkAQzT2F3FeBJXdG+ydrawGlfHAFIZG7wU6aKbQaRewysYpmrlW"
)
SWAGGER_UI_BUNDLE_JS_SRI = (
    "sha384-NXtFPpN61oWCuN4D42K6Zd5Rt2+uxeIT36R7kpXBuY9tLnZorzrJ4ykpqwJfgjpZ"
)

_SWAGGER_UI_CSS_URL = (
    f"https://unpkg.com/swagger-ui-dist@{SWAGGER_UI_VERSION}/swagger-ui.css"
)
_SWAGGER_UI_BUNDLE_JS_URL = (
    f"https://unpkg.com/swagger-ui-dist@{SWAGGER_UI_VERSION}/swagger-ui-bundle.js"
)

# Content-Security-Policy for the docs page.
# Allows styles/scripts only from unpkg with the pinned version, plus the
# local init.js endpoint.
_CSP_HEADER = (
    "default-src 'none'; "
    f"style-src https://unpkg.com/swagger-ui-dist@{SWAGGER_UI_VERSION}/; "
    f"script-src https://unpkg.com/swagger-ui-dist@{SWAGGER_UI_VERSION}/ 'self'; "
    "connect-src 'self'; "
    "img-src 'self' data:; "
    f"font-src https://unpkg.com/swagger-ui-dist@{SWAGGER_UI_VERSION}/; "
    "base-uri 'none'; "
    "form-action 'none'"
)


def _dev_docs_allowed() -> bool:
    """Check whether serving dev docs is permitted.

    The env var ``PYSCF_ENABLE_DEV_DOCS`` is evaluated at **request time** so
    that it can be toggled without restarting the process.

    Semantics:
    - Env var **absent** → docs are served (router is only registered in
      non-packaged / development mode, so this is safe).
    - Env var **present and == "1"** → docs are served.
    - Env var **present and != "1"** → docs are blocked (explicit opt-out).
    """
    return os.getenv("PYSCF_ENABLE_DEV_DOCS", "1") == "1"


@router.get("/api-docs", response_class=HTMLResponse)
def api_docs() -> Response:
    if not _dev_docs_allowed():
        return JSONResponse(
            {"success": False, "error": "Developer docs are disabled."},
            status_code=403,
        )

    html = f"""\
<!doctype html>
<html>
  <head>
    <title>PySCF Front API Docs</title>
    <link
      rel="stylesheet"
      href="{_SWAGGER_UI_CSS_URL}"
      integrity="{SWAGGER_UI_CSS_SRI}"
      crossorigin="anonymous"
    />
  </head>
  <body>
    <div id="swagger-ui"></div>
    <script
      src="{_SWAGGER_UI_BUNDLE_JS_URL}"
      integrity="{SWAGGER_UI_BUNDLE_JS_SRI}"
      crossorigin="anonymous"
    ></script>
    <script src="/api-docs/init.js"></script>
  </body>
</html>"""
    return HTMLResponse(
        content=html,
        headers={"Content-Security-Policy": _CSP_HEADER},
    )


@router.get("/api-docs/init.js")
def api_docs_init_js() -> Response:
    """Serve the Swagger UI initialization script as a static file."""
    if not _dev_docs_allowed():
        return JSONResponse(
            {"success": False, "error": "Developer docs are disabled."},
            status_code=403,
        )

    js_content = (
        'SwaggerUIBundle({ url: "/api-docs/spec.json", dom_id: "#swagger-ui" });'
    )
    return Response(
        content=js_content,
        media_type="application/javascript",
        headers={"Cache-Control": "public, max-age=3600"},
    )


@router.get("/api-docs/spec.json")
def api_docs_spec() -> JSONResponse:
    """Serve the OpenAPI spec as JSON.

    Handles missing/malformed YAML gracefully so that errors do not leak
    absolute filesystem paths to the client.
    """
    if not _dev_docs_allowed():
        return JSONResponse(
            {"success": False, "error": "Developer docs are disabled."},
            status_code=403,
        )

    spec_path = Path(__file__).resolve().parents[2] / "api-spec" / "openapi.yaml"
    try:
        with spec_path.open("r", encoding="utf-8") as handle:
            spec = yaml.safe_load(handle)
    except OSError:
        return JSONResponse(
            {"success": False, "error": "OpenAPI spec file not found or unreadable."},
            status_code=404,
        )
    except yaml.YAMLError:
        return JSONResponse(
            {"success": False, "error": "OpenAPI spec file is malformed."},
            status_code=500,
        )
    return JSONResponse(spec)
