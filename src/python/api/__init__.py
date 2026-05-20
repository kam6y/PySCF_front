# API Blueprints for PySCF Front Backend
# This module contains all API endpoint definitions organized by responsibility

import os

# Import all blueprint modules
from .health import health_bp
from .pubchem import pubchem_bp
from .smiles import smiles_bp
from .settings import settings_bp
from .system import system_bp
from .quantum import quantum_bp
from .agent import agent_bp
from .chat_history import chat_history_bp

# List of all blueprints to register with the main app
all_blueprints = [
    health_bp,
    pubchem_bp,
    smiles_bp,
    settings_bp,
    system_bp,
    quantum_bp,
    agent_bp,
    chat_history_bp
]


def register_blueprints(app):
    """
    Register all API blueprints with the Flask app.

    Swagger UI is conditionally registered in development mode only.
    Development mode is detected by checking if PYSCF_RESOURCES_PATH
    environment variable is not set (set only by Electron in packaged mode).
    """
    # Register all standard blueprints
    for blueprint in all_blueprints:
        app.register_blueprint(blueprint)

    # Conditionally register Swagger UI in development mode only
    is_packaged = os.getenv('PYSCF_RESOURCES_PATH') is not None
    if not is_packaged:
        from .swagger_ui import swagger_bp
        app.register_blueprint(swagger_bp)
        app.logger.info("✓ Swagger UI registered at /api-docs/ (development mode)")
    else:
        app.logger.info("✗ Swagger UI not registered (packaged mode)")
