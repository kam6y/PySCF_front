import logging
import os

from fastapi import FastAPI

from .agent import router as agent_router
from .calculation_updates import router as calculation_updates_router
from .chat_history import router as chat_history_router
from .health import router as health_router
from .pubchem import router as pubchem_router
from .quantum import router as quantum_router
from .settings import router as settings_router
from .smiles import router as smiles_router
from .swagger_ui import router as swagger_router
from .system import router as system_router

logger = logging.getLogger(__name__)


def register_routers(app: FastAPI) -> None:
    app.include_router(health_router)
    app.include_router(pubchem_router)
    app.include_router(smiles_router)
    app.include_router(settings_router)
    app.include_router(system_router)
    app.include_router(quantum_router)
    app.include_router(calculation_updates_router)
    app.include_router(agent_router)
    app.include_router(chat_history_router)
    is_packaged = os.getenv('PYSCF_RESOURCES_PATH') is not None
    if not is_packaged:
        app.include_router(swagger_router)
        logger.info("Swagger UI registered at /api-docs/ (development mode)")
    else:
        logger.info("Swagger UI not registered (packaged mode)")
    logger.info("Registered startup FastAPI routers")
