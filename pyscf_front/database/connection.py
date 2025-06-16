"""
PySCF_Front データベース接続管理
"""
import os
from typing import Optional
from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.pool import QueuePool
from loguru import logger
from dotenv import load_dotenv

from .models import Base

# 環境変数を読み込み
load_dotenv()


class DatabaseConfig:
    """データベース設定クラス"""
    
    def __init__(self):
        self.db_type = os.getenv('DB_TYPE', 'sqlite').lower()
        
        # SQLite設定
        self.db_path = os.getenv('DB_PATH', 'data/pyscf_front.db')
        
        # MySQL設定
        self.host = os.getenv('DB_HOST', 'localhost')
        self.port = int(os.getenv('DB_PORT', '3306'))
        self.database = os.getenv('DB_NAME', 'pyscf_front')
        self.username = os.getenv('DB_USER', 'root')
        self.password = os.getenv('DB_PASSWORD', '')
        self.charset = os.getenv('DB_CHARSET', 'utf8mb4')
        self.pool_size = int(os.getenv('DB_POOL_SIZE', '10'))
        self.max_overflow = int(os.getenv('DB_MAX_OVERFLOW', '20'))
        self.pool_timeout = int(os.getenv('DB_POOL_TIMEOUT', '30'))
        self.pool_recycle = int(os.getenv('DB_POOL_RECYCLE', '3600'))
    
    @property
    def connection_url(self) -> str:
        """データベース接続URLを生成"""
        if self.db_type == 'sqlite':
            # SQLiteの場合は絶対パスに変換
            import os.path
            if not os.path.isabs(self.db_path):
                # 相対パスの場合はプロジェクトルートからのパスとして扱う
                project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                self.db_path = os.path.join(project_root, self.db_path)
            
            # ディレクトリが存在しない場合は作成
            db_dir = os.path.dirname(self.db_path)
            if not os.path.exists(db_dir):
                os.makedirs(db_dir, exist_ok=True)
                
            return f"sqlite:///{self.db_path}"
        else:
            # MySQL
            return (
                f"mysql+mysqlconnector://{self.username}:{self.password}"
                f"@{self.host}:{self.port}/{self.database}"
                f"?charset={self.charset}"
            )


class DatabaseManager:
    """データベース接続管理クラス"""
    
    def __init__(self, config: Optional[DatabaseConfig] = None):
        self.config = config or DatabaseConfig()
        self._engine = None
        self._session_factory = None
    
    @property
    def engine(self):
        """SQLAlchemyエンジンを取得"""
        if self._engine is None:
            if self.config.db_type == 'sqlite':
                # SQLiteの場合はプーリング設定を無効化
                self._engine = create_engine(
                    self.config.connection_url,
                    echo=os.getenv('DB_ECHO', 'false').lower() == 'true'
                )
                logger.info(f"Database engine created: SQLite at {self.config.db_path}")
            else:
                # MySQLの場合は従来の設定
                self._engine = create_engine(
                    self.config.connection_url,
                    poolclass=QueuePool,
                    pool_size=self.config.pool_size,
                    max_overflow=self.config.max_overflow,
                    pool_timeout=self.config.pool_timeout,
                    pool_recycle=self.config.pool_recycle,
                    echo=os.getenv('DB_ECHO', 'false').lower() == 'true'
                )
                logger.info(f"Database engine created: MySQL at {self.config.host}:{self.config.port}")
        return self._engine
    
    @property
    def session_factory(self):
        """セッションファクトリーを取得"""
        if self._session_factory is None:
            self._session_factory = sessionmaker(bind=self.engine)
        return self._session_factory
    
    def get_session(self) -> Session:
        """新しいセッションを作成"""
        return self.session_factory()
    
    def create_tables(self):
        """テーブルを作成"""
        try:
            Base.metadata.create_all(self.engine)
            logger.info("Database tables created successfully")
        except Exception as e:
            logger.error(f"Failed to create database tables: {e}")
            raise
    
    def drop_tables(self):
        """テーブルを削除（開発用）"""
        try:
            Base.metadata.drop_all(self.engine)
            logger.info("Database tables dropped successfully")
        except Exception as e:
            logger.error(f"Failed to drop database tables: {e}")
            raise
    
    def test_connection(self) -> bool:
        """データベース接続をテスト"""
        try:
            with self.engine.connect() as conn:
                conn.execute(text("SELECT 1"))
            logger.info("Database connection test successful")
            return True
        except Exception as e:
            logger.error(f"Database connection test failed: {e}")
            return False


# グローバルインスタンス
db_manager = DatabaseManager()


def get_db_session() -> Session:
    """データベースセッションを取得するヘルパー関数"""
    return db_manager.get_session()


def init_database():
    """データベースを初期化"""
    logger.info("Initializing database...")
    
    # 接続テスト
    if not db_manager.test_connection():
        raise ConnectionError("Cannot connect to database")
    
    # テーブル作成
    db_manager.create_tables()
    
    logger.info("Database initialization completed")


def reset_database():
    """データベースをリセット（開発用）"""
    logger.warning("Resetting database (development only)...")
    db_manager.drop_tables()
    db_manager.create_tables()
    logger.info("Database reset completed")