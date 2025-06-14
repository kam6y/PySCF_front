"""
PySCF_Front Backend Main Application
"""
import os
import threading
import logging
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Import gRPC service
from src.grpc_service import serve as serve_grpc
from src.database.connection import init_database
from src.core.pyscf_config import pyscf_config

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="PySCF_Front Backend",
    description="Quantum chemistry calculation backend",
    version="0.1.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify allowed origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "message": "PySCF_Front Backend API",
        "version": "0.1.0",
        "status": "running"
    }

@app.get("/health")
async def health():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "database": "connected",  # TODO: Add actual database check
        "pyscf": "available"      # TODO: Add actual PySCF check
    }

@app.get("/info")
async def info():
    """System information endpoint"""
    import sys
    import platform
    
    return {
        "python_version": sys.version,
        "platform": platform.platform(),
        "environment": os.getenv("ENVIRONMENT", "development"),
        "database_url": os.getenv("DATABASE_URL", "not_configured").split("@")[-1] if "@" in os.getenv("DATABASE_URL", "") else "not_configured"
    }

@app.on_event("startup")
async def startup_event():
    """Initialize application on startup"""
    logger.info("Starting PySCF_Front Backend...")
    
    # Initialize database
    db_success = init_database()
    if not db_success:
        logger.error("Failed to initialize database")
        raise RuntimeError("Database initialization failed")
    
    logger.info("Database initialized successfully")
    
    # Initialize PySCF
    logger.info(f"PySCF initialized with {len(pyscf_config.get_available_methods())} methods and {len(pyscf_config.get_available_basis_sets())} basis sets")
    
    # Start gRPC server
    grpc_thread = threading.Thread(target=serve_grpc, daemon=True)
    grpc_thread.start()
    logger.info("gRPC server started in background thread")
    
    logger.info("PySCF_Front Backend startup complete")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)