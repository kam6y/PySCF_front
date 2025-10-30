# API Blueprints for PySCF Front Backend
# This module contains all API endpoint definitions organized by responsibility

from flask import Blueprint

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
    """Register all API blueprints with the Flask app."""
    for blueprint in all_blueprints:
        app.register_blueprint(blueprint)