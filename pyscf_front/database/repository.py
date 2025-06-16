"""
PySCF_Front データベースリポジトリ
データベース操作を抽象化するDAOパターン実装
"""
from typing import List, Optional, Dict, Any
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy.exc import SQLAlchemyError
from loguru import logger

from .models import (
    Instance, Molecule, Calculation, Result, JobQueue,
    MCPInteraction, MCPRecommendation,
    InstanceStatus, CalculationStatus, JobStatus
)
from .connection import get_db_session


class BaseRepository:
    """ベースリポジトリクラス"""
    
    def __init__(self, session: Optional[Session] = None):
        self.session = session or get_db_session()
        self._should_close_session = session is None
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):  # type: ignore
        if self._should_close_session:
            self.session.close()
    
    def commit(self):
        """変更をコミット"""
        try:
            self.session.commit()
        except SQLAlchemyError as e:
            self.session.rollback()
            raise e
    
    def rollback(self):
        """変更をロールバック"""
        self.session.rollback()


class InstanceRepository(BaseRepository):
    """インスタンス関連のデータベース操作"""
    
    def create(self, name: str, description: Optional[str] = None, 
               user_id: Optional[str] = None, project_id: Optional[str] = None) -> Instance:
        """新しいインスタンスを作成"""
        instance = Instance(
            name=name,
            description=description,
            user_id=user_id,
            project_id=project_id
        )
        self.session.add(instance)
        self.commit()
        logger.info(f"Created instance: {instance.id}")
        return instance
    
    def get_by_id(self, instance_id: str) -> Optional[Instance]:
        """IDでインスタンスを取得"""
        return self.session.query(Instance).filter(Instance.id == instance_id).first()
    
    def get_all(self, user_id: Optional[str] = None) -> List[Instance]:
        """すべてのインスタンスを取得"""
        query = self.session.query(Instance)
        if user_id:
            query = query.filter(Instance.user_id == user_id)
        return query.order_by(Instance.created_at.desc()).all()
    
    def update_status(self, instance_id: str, status: InstanceStatus) -> bool:
        """インスタンスのステータスを更新"""
        instance = self.get_by_id(instance_id)
        if instance:
            instance.status = status  # type: ignore
            self.commit()
            logger.info(f"Updated instance {instance_id} status to {status.value}")
            return True
        return False
    
    def delete(self, instance_id: str) -> bool:
        """インスタンスを削除"""
        instance = self.get_by_id(instance_id)
        if instance:
            self.session.delete(instance)
            self.commit()
            logger.info(f"Deleted instance: {instance_id}")
            return True
        return False


class MoleculeRepository(BaseRepository):
    """分子関連のデータベース操作"""
    
    def create(self, instance_id: str, name: str, formula: Optional[str] = None,
               geometry_data: Optional[Dict[str, Any]] = None, charge: int = 0,
               multiplicity: int = 1) -> Molecule:
        """新しい分子を作成"""
        molecule = Molecule(
            instance_id=instance_id,
            name=name,
            formula=formula,
            geometry_data=geometry_data,
            charge=charge,
            multiplicity=multiplicity
        )
        self.session.add(molecule)
        self.commit()
        logger.info(f"Created molecule: {molecule.id}")
        return molecule
    
    def get_by_instance(self, instance_id: str) -> List[Molecule]:
        """インスタンスIDで分子を取得"""
        return self.session.query(Molecule).filter(
            Molecule.instance_id == instance_id
        ).all()
    
    def get_by_id(self, molecule_id: str) -> Optional[Molecule]:
        """IDで分子を取得"""
        return self.session.query(Molecule).filter(Molecule.id == molecule_id).first()
    
    def update_geometry(self, molecule_id: str, geometry_data: Dict[str, Any]) -> bool:
        """分子の形状データを更新"""
        molecule = self.get_by_id(molecule_id)
        if molecule:
            molecule.geometry_data = geometry_data  # type: ignore
            self.commit()
            return True
        return False


class CalculationRepository(BaseRepository):
    """計算関連のデータベース操作"""
    
    def create(self, instance_id: str, method: str, basis_set: str,
               parameters: Optional[Dict[str, Any]] = None) -> Calculation:
        """新しい計算を作成"""
        calculation = Calculation(
            instance_id=instance_id,
            method=method,
            basis_set=basis_set,
            parameters=parameters
        )
        self.session.add(calculation)
        self.commit()
        logger.info(f"Created calculation: {calculation.id}")
        return calculation
    
    def get_by_instance(self, instance_id: str) -> List[Calculation]:
        """インスタンスIDで計算を取得"""
        return self.session.query(Calculation).filter(
            Calculation.instance_id == instance_id
        ).order_by(Calculation.start_time.desc()).all()
    
    def get_by_id(self, calculation_id: str) -> Optional[Calculation]:
        """IDで計算を取得"""
        return self.session.query(Calculation).filter(
            Calculation.id == calculation_id
        ).first()
    
    def update_status(self, calculation_id: str, status: CalculationStatus,
                     error_message: Optional[str] = None) -> bool:
        """計算のステータスを更新"""
        calculation = self.get_by_id(calculation_id)
        if calculation:
            calculation.status = status  # type: ignore
            if error_message:
                calculation.error_message = error_message  # type: ignore
            self.commit()
            logger.info(f"Updated calculation {calculation_id} status to {status.value}")
            return True
        return False
    
    def get_by_status(self, status: CalculationStatus) -> List[Calculation]:
        """ステータスで計算を取得"""
        return self.session.query(Calculation).filter(
            Calculation.status == status
        ).all()
    
    def update_start_time(self, calculation_id: str, start_time: datetime) -> bool:
        """計算開始時刻を更新"""
        calculation = self.get_by_id(calculation_id)
        if calculation:
            # SQLAlchemyのupdate文を使用してフィールドを更新
            self.session.query(Calculation).filter(
                Calculation.id == calculation_id
            ).update({'start_time': start_time})
            self.commit()
            return True
        return False
    
    def update_end_time(self, calculation_id: str, end_time: datetime) -> bool:
        """計算終了時刻を更新"""
        calculation = self.get_by_id(calculation_id)
        if calculation:
            # SQLAlchemyのupdate文を使用してフィールドを更新
            self.session.query(Calculation).filter(
                Calculation.id == calculation_id
            ).update({'end_time': end_time})
            self.commit()
            return True
        return False


class ResultRepository(BaseRepository):
    """結果関連のデータベース操作"""
    
    def create(self, calculation_id: str, result_type: str,
               result_data: Optional[Dict[str, Any]] = None,
               file_path: Optional[str] = None) -> Result:
        """新しい結果を作成"""
        result = Result(
            calculation_id=calculation_id,
            result_type=result_type,
            result_data=result_data,
            file_path=file_path
        )
        self.session.add(result)
        self.commit()
        logger.info(f"Created result: {result.id}")
        return result
    
    def get_by_calculation(self, calculation_id: str) -> List[Result]:
        """計算IDで結果を取得"""
        return self.session.query(Result).filter(
            Result.calculation_id == calculation_id
        ).order_by(Result.created_at.desc()).all()
    
    def get_by_type(self, calculation_id: str, result_type: str) -> Optional[Result]:
        """計算IDと結果タイプで結果を取得"""
        return self.session.query(Result).filter(
            Result.calculation_id == calculation_id,
            Result.result_type == result_type
        ).first()


class JobQueueRepository(BaseRepository):
    """ジョブキュー関連のデータベース操作"""
    
    def create(self, calculation_id: str, priority: int = 5) -> JobQueue:
        """新しいジョブを作成"""
        job = JobQueue(
            calculation_id=calculation_id,
            priority=priority
        )
        self.session.add(job)
        self.commit()
        logger.info(f"Created job: {job.id}")
        return job
    
    def get_next_job(self) -> Optional[JobQueue]:
        """次に実行すべきジョブを取得"""
        return self.session.query(JobQueue).filter(
            JobQueue.status == JobStatus.WAITING
        ).order_by(JobQueue.priority.desc(), JobQueue.created_at.asc()).first()
    
    def update_status(self, job_id: str, status: JobStatus,
                     assigned_worker: Optional[str] = None) -> bool:
        """ジョブのステータスを更新"""
        job = self.session.query(JobQueue).filter(JobQueue.id == job_id).first()
        if job:
            job.status = status  # type: ignore
            if assigned_worker:
                job.assigned_worker = assigned_worker  # type: ignore
            self.commit()
            return True
        return False
    
    def update_completed_time(self, job_id: str, completed_time: datetime) -> bool:
        """ジョブ完了時刻を更新"""
        job = self.session.query(JobQueue).filter(JobQueue.id == job_id).first()
        if job:
            # SQLAlchemyのupdate文を使用してフィールドを更新
            self.session.query(JobQueue).filter(
                JobQueue.id == job_id
            ).update({'completed_at': completed_time})
            self.commit()
            return True
        return False


class MCPRepository(BaseRepository):
    """MCP関連のデータベース操作"""
    
    def create_interaction(self, instance_id: str, interaction_type: str,
                          user_query: str, mcp_response: Dict[str, Any]) -> MCPInteraction:
        """新しいMCPインタラクションを作成"""
        interaction = MCPInteraction(
            instance_id=instance_id,
            interaction_type=interaction_type,
            user_query=user_query,
            mcp_response=mcp_response
        )
        self.session.add(interaction)
        self.commit()
        return interaction
    
    def create_recommendation(self, calculation_id: str, recommendation_type: str,
                            recommendation_data: Dict[str, Any],
                            confidence_score: float) -> MCPRecommendation:
        """新しいMCP推奨を作成"""
        recommendation = MCPRecommendation(
            calculation_id=calculation_id,
            recommendation_type=recommendation_type,
            recommendation_data=recommendation_data,
            confidence_score=confidence_score
        )
        self.session.add(recommendation)
        self.commit()
        return recommendation


# ユーティリティ関数
def create_complete_instance(name: str, description: Optional[str] = None,
                           molecule_data: Optional[Dict[str, Any]] = None,
                           calculation_config: Optional[Dict[str, Any]] = None) -> str:
    """分子と計算設定を含む完全なインスタンスを作成"""
    with InstanceRepository() as instance_repo:
        # インスタンス作成
        instance = instance_repo.create(name, description)
        
        # 分子作成
        if molecule_data:
            with MoleculeRepository(instance_repo.session) as mol_repo:
                mol_repo.create(
                    instance_id=str(instance.id),
                    name=molecule_data.get('name', 'Untitled'),
                    formula=molecule_data.get('formula'),
                    geometry_data=molecule_data.get('geometry_data'),
                    charge=molecule_data.get('charge', 0),
                    multiplicity=molecule_data.get('multiplicity', 1)
                )
        
        # 計算設定作成
        if calculation_config:
            with CalculationRepository(instance_repo.session) as calc_repo:
                calc_repo.create(
                    instance_id=str(instance.id),
                    method=calculation_config['method'],
                    basis_set=calculation_config['basis_set'],
                    parameters=calculation_config.get('parameters')
                )
        
        return str(instance.id)