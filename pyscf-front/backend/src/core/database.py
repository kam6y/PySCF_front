"""Database management for PySCF_Front"""

import logging
from typing import Optional, Dict, Any
from contextlib import asynccontextmanager

from sqlalchemy import create_engine, text
from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession
from sqlalchemy.orm import sessionmaker, declarative_base
from sqlalchemy.pool import QueuePool

from .config import get_database_url, settings

logger = logging.getLogger(__name__)

# SQLAlchemy Base
Base = declarative_base()

class DatabaseManager:
    """Database connection and session management"""
    
    _engine = None
    _async_engine = None
    _SessionLocal = None
    _AsyncSessionLocal = None
    
    @classmethod
    async def initialize(cls):
        """Initialize database connections"""
        try:
            # Synchronous engine for simple operations
            cls._engine = create_engine(
                get_database_url(),
                poolclass=QueuePool,
                pool_size=20,
                max_overflow=0,
                pool_pre_ping=True,
                pool_recycle=3600,
                echo=settings.DATABASE_ECHO
            )
            
            cls._SessionLocal = sessionmaker(
                autocommit=False,
                autoflush=False,
                bind=cls._engine
            )
            
            # Test connection
            with cls._engine.connect() as conn:
                conn.execute(text("SELECT 1"))
            
            logger.info("Database initialized successfully")
            
        except Exception as e:
            logger.error(f"Database initialization failed: {e}")
            # For development, we'll continue without database
            logger.warning("Continuing without database connection (development mode)")
    
    @classmethod
    async def close(cls):
        """Close database connections"""
        if cls._engine:
            cls._engine.dispose()
        if cls._async_engine:
            await cls._async_engine.dispose()
        logger.info("Database connections closed")
    
    @classmethod
    async def health_check(cls) -> Dict[str, Any]:
        """Check database health"""
        if not cls._engine:
            return {
                "status": "disconnected",
                "error": "Database not initialized"
            }
        
        try:
            with cls._engine.connect() as conn:
                result = conn.execute(text("SELECT 1 as health_check"))
                row = result.fetchone()
                
                if row and row[0] == 1:
                    return {
                        "status": "healthy",
                        "connection": "ok"
                    }
                else:
                    return {
                        "status": "unhealthy",
                        "error": "Unexpected result from health check"
                    }
                    
        except Exception as e:
            return {
                "status": "unhealthy",
                "error": str(e)
            }
    
    @classmethod
    @asynccontextmanager
    async def get_session(cls):
        """Get database session context manager"""
        if not cls._SessionLocal:
            raise RuntimeError("Database not initialized")
        
        session = cls._SessionLocal()
        try:
            yield session
            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()
    
    @classmethod
    def get_sync_session(cls):
        """Get synchronous database session"""
        if not cls._SessionLocal:
            raise RuntimeError("Database not initialized")
        return cls._SessionLocal()


# Database models
from sqlalchemy import Column, String, Integer, Float, Text, DateTime, Boolean, JSON, ForeignKey
from sqlalchemy.sql import func
from sqlalchemy.orm import relationship
import uuid


def generate_uuid():
    """Generate UUID string"""
    return str(uuid.uuid4())


class Instance(Base):
    """Calculation instance model"""
    __tablename__ = "instances"
    
    id = Column(String(36), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    description = Column(Text)
    status = Column(String(20), default="draft")  # draft, ready, running, completed, error
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now())
    user_id = Column(String(36), nullable=True)
    project_id = Column(String(36), nullable=True)
    
    # Relationships
    molecules = relationship("Molecule", back_populates="instance", cascade="all, delete-orphan")
    calculations = relationship("Calculation", back_populates="instance", cascade="all, delete-orphan")


class Molecule(Base):
    """Molecule model"""
    __tablename__ = "molecules"
    
    id = Column(String(36), primary_key=True, default=generate_uuid)
    instance_id = Column(String(36), ForeignKey("instances.id"), nullable=False)
    name = Column(String(255))
    formula = Column(String(100))
    molecular_weight = Column(Float)
    geometry_type = Column(String(20), default="xyz")  # xyz, zmatrix
    geometry_data = Column(JSON)
    charge = Column(Integer, default=0)
    multiplicity = Column(Integer, default=1)
    symmetry = Column(String(10))
    
    # Relationships
    instance = relationship("Instance", back_populates="molecules")


class Calculation(Base):
    """Calculation model"""
    __tablename__ = "calculations"
    
    id = Column(String(36), primary_key=True, default=generate_uuid)
    instance_id = Column(String(36), ForeignKey("instances.id"), nullable=False)
    method = Column(String(50), nullable=False)
    basis_set = Column(String(50), nullable=False)
    parameters = Column(JSON)
    convergence_criteria = Column(JSON)
    max_iterations = Column(Integer, default=100)
    start_time = Column(DateTime)
    end_time = Column(DateTime)
    status = Column(String(20), default="pending")  # pending, running, completed, failed, cancelled
    error_message = Column(Text)
    
    # Relationships
    instance = relationship("Instance", back_populates="calculations")
    results = relationship("Result", back_populates="calculation", cascade="all, delete-orphan")


class Result(Base):
    """Calculation result model"""
    __tablename__ = "results"
    
    id = Column(String(36), primary_key=True, default=generate_uuid)
    calculation_id = Column(String(36), ForeignKey("calculations.id"), nullable=False)
    result_type = Column(String(50), nullable=False)  # energy, orbitals, properties, etc.
    result_data = Column(JSON)
    file_path = Column(String(500))
    file_size = Column(Integer)
    checksum = Column(String(64))
    created_at = Column(DateTime, default=func.now())
    
    # Relationships
    calculation = relationship("Calculation", back_populates="results")


class JobQueue(Base):
    """Job queue model"""
    __tablename__ = "job_queue"
    
    id = Column(String(36), primary_key=True, default=generate_uuid)
    calculation_id = Column(String(36), ForeignKey("calculations.id"), nullable=False)
    priority = Column(Integer, default=5)
    status = Column(String(20), default="waiting")  # waiting, running, completed, failed
    assigned_worker = Column(String(100))
    created_at = Column(DateTime, default=func.now())
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    
    # Relationships
    calculation = relationship("Calculation")


# Initialize database tables
def create_tables():
    """Create all database tables"""
    if DatabaseManager._engine:
        Base.metadata.create_all(bind=DatabaseManager._engine)
        logger.info("Database tables created successfully")
    else:
        logger.warning("Cannot create tables: database not initialized")