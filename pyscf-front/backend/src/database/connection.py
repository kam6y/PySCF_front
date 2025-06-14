"""
Database connection and configuration
"""
import os
import logging
from typing import Optional
from sqlalchemy import create_engine, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.pool import StaticPool
import mysql.connector
from mysql.connector import Error

logger = logging.getLogger(__name__)

# Database URL from environment
DATABASE_URL = os.getenv(
    'DATABASE_URL', 
    'mysql+pymysql://pyscf_user:pyscf_password@localhost:3307/pyscf_dev'
)

# SQLAlchemy configuration
engine = create_engine(
    DATABASE_URL,
    echo=os.getenv('DATABASE_ECHO', 'false').lower() == 'true',
    pool_pre_ping=True,
    pool_recycle=300,
    connect_args={
        "check_same_thread": False,
        "charset": "utf8mb4"
    } if "sqlite" in DATABASE_URL else {"charset": "utf8mb4"}
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

class DatabaseManager:
    """Database connection manager"""
    
    def __init__(self):
        self.engine = engine
        self.SessionLocal = SessionLocal
        self.Base = Base
    
    def get_session(self) -> Session:
        """Get database session"""
        return self.SessionLocal()
    
    def close_session(self, session: Session):
        """Close database session"""
        session.close()
    
    def test_connection(self) -> tuple[bool, str]:
        """Test database connection"""
        try:
            from sqlalchemy import text
            with self.engine.connect() as connection:
                result = connection.execute(text("SELECT 1"))
                row = result.fetchone()
                if row and row[0] == 1:
                    return True, "Database connection successful"
                else:
                    return False, "Unexpected result from test query"
        except Exception as e:
            logger.error(f"Database connection test failed: {e}")
            return False, f"Connection failed: {str(e)}"
    
    def create_tables(self):
        """Create all tables"""
        try:
            self.Base.metadata.create_all(bind=self.engine)
            logger.info("Database tables created successfully")
        except Exception as e:
            logger.error(f"Failed to create tables: {e}")
            raise
    
    def get_engine_info(self) -> dict:
        """Get engine information"""
        return {
            'url': str(self.engine.url).replace(self.engine.url.password or '', '***'),
            'driver': self.engine.dialect.name,
            'pool_size': self.engine.pool.size(),
            'checked_in': self.engine.pool.checkedin(),
            'checked_out': self.engine.pool.checkedout(),
        }

# Global database manager instance
db_manager = DatabaseManager()

def get_db() -> Session:
    """Dependency to get database session"""
    db = db_manager.get_session()
    try:
        yield db
    finally:
        db.close()

def init_database():
    """Initialize database"""
    try:
        # Test connection
        success, message = db_manager.test_connection()
        if not success:
            raise RuntimeError(f"Database connection failed: {message}")
        
        # Create tables
        db_manager.create_tables()
        
        logger.info("Database initialized successfully")
        return True
    except Exception as e:
        logger.error(f"Database initialization failed: {e}")
        return False