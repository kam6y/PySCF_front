"""
PySCF_Front Backend Main Application
"""
import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

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

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)