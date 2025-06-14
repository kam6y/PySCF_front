"""Configuration management for PySCF_Front Backend"""

import os
from typing import Optional

class Settings:
    """Application settings"""
    
    # Application
    DEBUG: bool = os.getenv("DEBUG", "false").lower() == "true"
    API_PORT: int = int(os.getenv("API_PORT", "8000"))
    GRPC_PORT: int = int(os.getenv("GRPC_PORT", "50051"))
    
    # Database
    DATABASE_URL: str = os.getenv(
        "DATABASE_URL", 
        "mysql+pymysql://pyscf_user:pyscf_password@localhost:3306/pyscf_dev"
    )
    DATABASE_ECHO: bool = DEBUG
    
    # PySCF Configuration
    PYSCF_TMPDIR: str = os.getenv("PYSCF_TMPDIR", "/tmp/pyscf")
    PYSCF_MAX_MEMORY: int = int(os.getenv("PYSCF_MAX_MEMORY", "4000"))  # MB
    OMP_NUM_THREADS: int = int(os.getenv("OMP_NUM_THREADS", "4"))
    
    # GPU Configuration
    USE_GPU: bool = os.getenv("PYSCF_USE_GPU", "false").lower() == "true"
    CUDA_VISIBLE_DEVICES: Optional[str] = os.getenv("CUDA_VISIBLE_DEVICES")
    
    # Security
    SECRET_KEY: str = os.getenv("SECRET_KEY", "dev-secret-key")
    ACCESS_TOKEN_EXPIRE_MINUTES: int = int(os.getenv("ACCESS_TOKEN_EXPIRE_MINUTES", "30"))
    
    # MCP Server (Optional)
    MCP_SERVER_ENABLED: bool = os.getenv("MCP_SERVER_ENABLED", "false").lower() == "true"
    MCP_SERVER_PORT: int = int(os.getenv("MCP_SERVER_PORT", "50053"))
    
    # Logging
    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")
    LOG_FILE: Optional[str] = os.getenv("LOG_FILE")
    


# Global settings instance
settings = Settings()


def get_database_url() -> str:
    """Get database URL with proper formatting"""
    return settings.DATABASE_URL


def get_pyscf_config() -> dict:
    """Get PySCF configuration dictionary"""
    return {
        "tmpdir": settings.PYSCF_TMPDIR,
        "max_memory": settings.PYSCF_MAX_MEMORY,
        "num_threads": settings.OMP_NUM_THREADS,
        "use_gpu": settings.USE_GPU,
        "cuda_devices": settings.CUDA_VISIBLE_DEVICES,
    }


def is_development() -> bool:
    """Check if running in development mode"""
    return settings.DEBUG


def is_gpu_enabled() -> bool:
    """Check if GPU acceleration is enabled"""
    return settings.USE_GPU


def is_mcp_enabled() -> bool:
    """Check if MCP server is enabled"""
    return settings.MCP_SERVER_ENABLED