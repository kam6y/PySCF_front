"""PySCF_Front Backend Main Entry Point"""

import asyncio
import logging
import os
from concurrent import futures
from contextlib import asynccontextmanager

import grpc
import uvicorn
from fastapi import FastAPI

from api.grpc_server import CalculationServicer
from core.config import settings
from core.database import DatabaseManager
from core.pyscf_engine import PySCFEngine


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager"""
    # Startup
    logger.info("Starting PySCF_Front Backend...")
    
    # Initialize database
    await DatabaseManager.initialize()
    
    # Start gRPC server
    grpc_server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
    servicer = CalculationServicer()
    
    # Import the generated gRPC service
    from generated import calculation_pb2_grpc
    calculation_pb2_grpc.add_CalculationServiceServicer_to_server(servicer, grpc_server)
    
    listen_addr = f'[::]:{settings.GRPC_PORT}'
    grpc_server.add_insecure_port(listen_addr)
    grpc_server.start()
    
    logger.info(f"gRPC server started on {listen_addr}")
    
    yield
    
    # Shutdown
    logger.info("Shutting down PySCF_Front Backend...")
    grpc_server.stop(grace=5)
    await DatabaseManager.close()


# FastAPI app for REST API and debugging
app = FastAPI(
    title="PySCF_Front Backend",
    description="Quantum chemistry calculation backend using PySCF",
    version="1.0.0",
    lifespan=lifespan
)


@app.get("/")
async def root():
    """Root endpoint"""
    return {"message": "PySCF_Front Backend is running", "version": "1.0.0"}


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    try:
        # Check database connection
        db_status = await DatabaseManager.health_check()
        
        # Check PySCF availability
        pyscf_status = PySCFEngine.health_check()
        
        return {
            "status": "healthy",
            "database": db_status,
            "pyscf": pyscf_status,
            "grpc_port": settings.GRPC_PORT
        }
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return {"status": "unhealthy", "error": str(e)}


@app.get("/system-info")
async def system_info():
    """System information endpoint"""
    try:
        info = PySCFEngine.get_system_info()
        return {
            "status": "success",
            "system_info": info
        }
    except Exception as e:
        logger.error(f"System info failed: {e}")
        return {"status": "error", "error": str(e)}


if __name__ == "__main__":
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=settings.API_PORT,
        reload=settings.DEBUG,
        log_level="info"
    )