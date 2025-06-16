"""
PySCF_Front データベースモデル
SQLAlchemyを使用したORM定義
"""
import uuid
from datetime import datetime
from typing import Optional, Dict, Any
from sqlalchemy import (
    create_engine, Column, String, Integer, Text, JSON, TIMESTAMP,
    Enum, Boolean, BigInteger, ForeignKey, Index
)
from sqlalchemy.types import Numeric
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, Session
from sqlalchemy.dialects.mysql import VARCHAR
import enum

Base = declarative_base()


class InstanceStatus(enum.Enum):
    """インスタンスステータス"""
    DRAFT = "draft"
    READY = "ready"
    RUNNING = "running"
    COMPLETED = "completed"
    ERROR = "error"


class CalculationStatus(enum.Enum):
    """計算ステータス"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobStatus(enum.Enum):
    """ジョブステータス"""
    WAITING = "waiting"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class GeometryType(enum.Enum):
    """分子形状データタイプ"""
    XYZ = "xyz"
    ZMATRIX = "zmatrix"


class InteractionType(enum.Enum):
    """MCPインタラクションタイプ"""
    MOLECULE_GENERATION = "molecule_generation"
    PARAMETER_RECOMMENDATION = "parameter_recommendation"
    RESULT_INTERPRETATION = "result_interpretation"


class Instance(Base):
    """インスタンス管理テーブル"""
    __tablename__ = 'instances'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    name = Column(String(255), nullable=False)
    description = Column(Text)
    created_at = Column(TIMESTAMP, default=datetime.utcnow)
    updated_at = Column(TIMESTAMP, default=datetime.utcnow, onupdate=datetime.utcnow)
    status = Column(Enum(InstanceStatus), default=InstanceStatus.DRAFT)
    user_id = Column(VARCHAR(36))
    project_id = Column(VARCHAR(36))
    
    # リレーションシップ
    molecules = relationship("Molecule", back_populates="instance", cascade="all, delete-orphan")
    calculations = relationship("Calculation", back_populates="instance", cascade="all, delete-orphan")
    mcp_interactions = relationship("MCPInteraction", back_populates="instance", cascade="all, delete-orphan")
    
    # インデックス
    __table_args__ = (
        Index('idx_status', 'status'),
        Index('idx_user_id', 'user_id'),
    )

    def __repr__(self):
        return f"<Instance(id='{self.id}', name='{self.name}', status='{self.status.value}')>"


class Molecule(Base):
    """分子情報テーブル"""
    __tablename__ = 'molecules'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    instance_id = Column(VARCHAR(36), ForeignKey('instances.id', ondelete="CASCADE"), nullable=False)
    name = Column(String(255))
    formula = Column(String(100))
    molecular_weight = Column(Numeric(10, 4))
    geometry_type = Column(Enum(GeometryType))
    geometry_data = Column(JSON)
    charge = Column(Integer, default=0)
    multiplicity = Column(Integer, default=1)
    symmetry = Column(String(10))
    
    # リレーションシップ
    instance = relationship("Instance", back_populates="molecules")
    
    # インデックス
    __table_args__ = (
        Index('idx_instance_id', 'instance_id'),
    )

    def __repr__(self):
        return f"<Molecule(id='{self.id}', name='{self.name}', formula='{self.formula}')>"


class Calculation(Base):
    """計算設定テーブル"""
    __tablename__ = 'calculations'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    instance_id = Column(VARCHAR(36), ForeignKey('instances.id', ondelete="CASCADE"), nullable=False)
    method = Column(String(50), nullable=False)
    basis_set = Column(String(50), nullable=False)
    parameters = Column(JSON)
    convergence_criteria = Column(JSON)
    max_iterations = Column(Integer, default=100)
    start_time = Column(TIMESTAMP, nullable=True)
    end_time = Column(TIMESTAMP, nullable=True)
    status = Column(Enum(CalculationStatus), default=CalculationStatus.PENDING)
    error_message = Column(Text)
    
    # リレーションシップ
    instance = relationship("Instance", back_populates="calculations")
    results = relationship("Result", back_populates="calculation", cascade="all, delete-orphan")
    job_queue = relationship("JobQueue", back_populates="calculation", cascade="all, delete-orphan")
    mcp_recommendations = relationship("MCPRecommendation", back_populates="calculation", cascade="all, delete-orphan")
    
    # インデックス
    __table_args__ = (
        Index('idx_instance_status', 'instance_id', 'status'),
    )

    def __repr__(self):
        return f"<Calculation(id='{self.id}', method='{self.method}', basis_set='{self.basis_set}', status='{self.status.value}')>"


class Result(Base):
    """計算結果テーブル"""
    __tablename__ = 'results'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    calculation_id = Column(VARCHAR(36), ForeignKey('calculations.id', ondelete="CASCADE"), nullable=False)
    result_type = Column(String(50), nullable=False)
    result_data = Column(JSON)
    file_path = Column(String(500))
    file_size = Column(BigInteger)
    checksum = Column(String(64))
    created_at = Column(TIMESTAMP, default=datetime.utcnow)
    
    # リレーションシップ
    calculation = relationship("Calculation", back_populates="results")
    
    # インデックス
    __table_args__ = (
        Index('idx_calculation_type', 'calculation_id', 'result_type'),
    )

    def __repr__(self):
        return f"<Result(id='{self.id}', type='{self.result_type}', calculation_id='{self.calculation_id}')>"


class JobQueue(Base):
    """ジョブキュー管理テーブル"""
    __tablename__ = 'job_queue'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    calculation_id = Column(VARCHAR(36), ForeignKey('calculations.id'), nullable=False)
    priority = Column(Integer, default=5)
    status = Column(Enum(JobStatus), default=JobStatus.WAITING)
    assigned_worker = Column(String(100))
    created_at = Column(TIMESTAMP, default=datetime.utcnow)
    started_at = Column(TIMESTAMP, nullable=True)
    completed_at = Column(TIMESTAMP, nullable=True)
    
    # リレーションシップ
    calculation = relationship("Calculation", back_populates="job_queue")
    
    # インデックス
    __table_args__ = (
        Index('idx_status_priority', 'status', 'priority'),
    )

    def __repr__(self):
        return f"<JobQueue(id='{self.id}', status='{self.status.value}', priority={self.priority})>"


class MCPInteraction(Base):
    """MCP関連のインタラクション記録テーブル"""
    __tablename__ = 'mcp_interactions'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    instance_id = Column(VARCHAR(36), ForeignKey('instances.id'), nullable=True)
    interaction_type = Column(Enum(InteractionType))
    user_query = Column(Text)
    mcp_response = Column(JSON)
    accepted = Column(Boolean, default=False)
    created_at = Column(TIMESTAMP, default=datetime.utcnow)
    
    # リレーションシップ
    instance = relationship("Instance", back_populates="mcp_interactions")
    
    # インデックス
    __table_args__ = (
        Index('idx_instance_type', 'instance_id', 'interaction_type'),
    )

    def __repr__(self):
        return f"<MCPInteraction(id='{self.id}', type='{self.interaction_type.value}', accepted={self.accepted})>"


class MCPRecommendation(Base):
    """MCP推奨設定の記録テーブル"""
    __tablename__ = 'mcp_recommendations'
    
    id = Column(VARCHAR(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    calculation_id = Column(VARCHAR(36), ForeignKey('calculations.id'), nullable=True)
    recommendation_type = Column(String(50))
    recommendation_data = Column(JSON)
    confidence_score = Column(Numeric(3, 2))
    applied = Column(Boolean, default=False)
    created_at = Column(TIMESTAMP, default=datetime.utcnow)
    
    # リレーションシップ
    calculation = relationship("Calculation", back_populates="mcp_recommendations")
    
    # インデックス
    __table_args__ = (
        Index('idx_calculation_applied', 'calculation_id', 'applied'),
    )

    def __repr__(self):
        return f"<MCPRecommendation(id='{self.id}', type='{self.recommendation_type}', applied={self.applied})>"