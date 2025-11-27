"""
Swagger UI blueprint for API documentation (development mode only).

This module provides interactive API documentation via Swagger UI.
It is only registered when running in development mode (npm run dev).
"""

from flask import Blueprint, jsonify, current_app
from pathlib import Path
import yaml

swagger_bp = Blueprint(
    'swagger',
    __name__,
    url_prefix='/api-docs'
)


@swagger_bp.route('/', methods=['GET'])
def swagger_ui():
    """Serve Swagger UI HTML."""
    return '''
    <!DOCTYPE html>
    <html>
      <head>
        <title>PySCF Native App API Documentation</title>
        <meta charset="utf-8"/>
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <link rel="stylesheet" href="https://unpkg.com/swagger-ui-dist@5/swagger-ui.css">
      </head>
      <body>
        <div id="swagger-ui"></div>
        <script src="https://unpkg.com/swagger-ui-dist@5/swagger-ui-bundle.js"></script>
        <script>
          window.onload = function() {
            SwaggerUIBundle({
              url: '/api-docs/spec.json',
              dom_id: '#swagger-ui',
              presets: [
                SwaggerUIBundle.presets.apis,
                SwaggerUIBundle.SwaggerUIStandalonePreset
              ],
              layout: "BaseLayout",
              deepLinking: true
            })
          }
        </script>
      </body>
    </html>
    ''', 200, {'Content-Type': 'text/html; charset=utf-8'}


@swagger_bp.route('/spec.json', methods=['GET'])
def swagger_spec():
    """Serve the OpenAPI specification as JSON."""
    try:
        # Get path to OpenAPI spec file
        spec_path = Path(__file__).parent.parent.parent / 'api-spec' / 'openapi.yaml'

        if not spec_path.exists():
            current_app.logger.error(f"OpenAPI spec not found at: {spec_path}")
            return jsonify({'error': 'OpenAPI spec not found'}), 404

        # Parse YAML to JSON
        with open(spec_path, 'r', encoding='utf-8') as f:
            spec = yaml.safe_load(f)

        return jsonify(spec)
    except Exception as e:
        current_app.logger.error(f"Error serving OpenAPI spec: {e}")
        return jsonify({'error': 'Failed to load OpenAPI spec'}), 500
