import logging
import os

from fastapi import FastAPI

from .health import router as health_router
from .swagger_ui import router as swagger_router

logger = logging.getLogger(__name__)


def register_routers(app: FastAPI) -> None:
    app.include_router(health_router)
    is_packaged = os.getenv('PYSCF_RESOURCES_PATH') is not None
    if not is_packaged:
        app.include_router(swagger_router)
        logger.info("Swagger UI registered at /api-docs/ (development mode)")
    else:
        logger.info("Swagger UI not registered (packaged mode)")
    logger.info("Registered startup FastAPI routers")
