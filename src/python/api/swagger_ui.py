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
