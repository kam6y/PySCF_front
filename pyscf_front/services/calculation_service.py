"""
計算管理サービス
既存のCalculationJobクラスとデータベースを統合
"""
from typing import List, Dict, Any, Optional
from datetime import datetime, timezone
from loguru import logger

from pyscf_front.core.calculation_engine_unified import CalculationJob, CalculationStatus as CoreCalculationStatus
from pyscf_front.database import (
    CalculationRepository, ResultRepository, JobQueueRepository, InstanceRepository,
    Calculation as DBCalculation, CalculationStatus, JobStatus
)
from pyscf_front.database.models import JobQueue


class CalculationService:
    """計算管理サービスクラス"""
    
    def __init__(self):
        self.calculation_repo = CalculationRepository()
        self.result_repo = ResultRepository()
        self.job_repo = JobQueueRepository()
        self.instance_repo = InstanceRepository()
    
    def create_calculation(self, instance_id: str, method: str, basis_set: str,
                          parameters: Optional[Dict[str, Any]] = None,
                          priority: int = 5) -> str:
        """新しい計算を作成してジョブキューに追加"""
        try:
            with self.calculation_repo as calc_repo:
                # 計算設定を作成
                calculation = calc_repo.create(
                    instance_id=instance_id,
                    method=method,
                    basis_set=basis_set,
                    parameters=parameters or {}
                )
                
                # ジョブキューに追加
                with JobQueueRepository(calc_repo.session) as job_repo:
                    job = job_repo.create(
                        calculation_id=str(calculation.id),
                        priority=priority
                    )
                
                # インスタンスステータスを更新
                from pyscf_front.database import InstanceStatus
                with InstanceRepository(calc_repo.session) as inst_repo:
                    inst_repo.update_status(instance_id, InstanceStatus.READY)
                
                logger.info(f"Created calculation: {calculation.id}")
                return str(calculation.id)
                
        except Exception as e:
            logger.error(f"Failed to create calculation: {e}")
            raise
    
    def get_calculation_by_id(self, calculation_id: str) -> Optional[Dict[str, Any]]:
        """計算IDで計算情報を取得"""
        try:
            with self.calculation_repo as repo:
                calculation = repo.get_by_id(calculation_id)
                if calculation:
                    return self._convert_db_calculation_to_dict(calculation)
                return None
                
        except Exception as e:
            logger.error(f"Failed to get calculation {calculation_id}: {e}")
            return None
    
    def get_calculations_by_instance(self, instance_id: str) -> List[Dict[str, Any]]:
        """インスタンスの計算一覧を取得"""
        try:
            with self.calculation_repo as repo:
                calculations = repo.get_by_instance(instance_id)
                return [self._convert_db_calculation_to_dict(calc) for calc in calculations]
                
        except Exception as e:
            logger.error(f"Failed to get calculations for instance {instance_id}: {e}")
            return []
    
    def start_calculation(self, calculation_id: str, worker_id: str = "default") -> bool:
        """計算を開始状態に更新"""
        try:
            with self.calculation_repo as calc_repo:
                # 計算ステータスを更新
                success = calc_repo.update_status(calculation_id, CalculationStatus.RUNNING)
                if not success:
                    return False
                
                # 開始時刻を更新
                calculation = calc_repo.get_by_id(calculation_id)
                if calculation:
                    calc_repo.update_start_time(calculation_id, datetime.now(timezone.utc))
                    calc_repo.commit()
                
                # ジョブキューのステータスを更新
                with JobQueueRepository(calc_repo.session) as job_repo:
                    # 対応するジョブを見つけて更新
                    job = calc_repo.session.query(JobQueue).filter_by(
                        calculation_id=calculation_id
                    ).first()
                    if job:
                        job_repo.update_status(str(job.id), JobStatus.RUNNING, worker_id)
                
                logger.info(f"Started calculation: {calculation_id}")
                return True
                
        except Exception as e:
            logger.error(f"Failed to start calculation {calculation_id}: {e}")
            return False
    
    def complete_calculation(self, calculation_id: str, results: Dict[str, Any]) -> bool:
        """計算を完了状態に更新し、結果を保存"""
        try:
            with self.calculation_repo as calc_repo:
                # 計算ステータスを更新
                success = calc_repo.update_status(calculation_id, CalculationStatus.COMPLETED)
                if not success:
                    return False
                
                # 終了時刻を更新
                calculation = calc_repo.get_by_id(calculation_id)
                if calculation:
                    calc_repo.update_end_time(calculation_id, datetime.now(timezone.utc))
                    calc_repo.commit()
                
                # 結果を保存
                with ResultRepository(calc_repo.session) as result_repo:
                    for result_type, result_data in results.items():
                        result_repo.create(
                            calculation_id=calculation_id,
                            result_type=result_type,
                            result_data=result_data
                        )
                
                # ジョブキューのステータスを更新
                with JobQueueRepository(calc_repo.session) as job_repo:
                    job = calc_repo.session.query(JobQueue).filter_by(
                        calculation_id=calculation_id
                    ).first()
                    if job:
                        job_repo.update_status(str(job.id), JobStatus.COMPLETED)
                        job_repo.update_completed_time(str(job.id), datetime.now(timezone.utc))
                        calc_repo.commit()
                
                logger.info(f"Completed calculation: {calculation_id}")
                return True
                
        except Exception as e:
            logger.error(f"Failed to complete calculation {calculation_id}: {e}")
            return False
    
    def fail_calculation(self, calculation_id: str, error_message: str) -> bool:
        """計算を失敗状態に更新"""
        try:
            with self.calculation_repo as calc_repo:
                success = calc_repo.update_status(
                    calculation_id, 
                    CalculationStatus.FAILED,
                    error_message
                )
                
                if success:
                    # ジョブキューのステータスを更新
                    with JobQueueRepository(calc_repo.session) as job_repo:
                        job = calc_repo.session.query(JobQueue).filter_by(
                            calculation_id=calculation_id
                        ).first()
                        if job:
                            job_repo.update_status(str(job.id), JobStatus.FAILED)
                    
                    logger.info(f"Failed calculation: {calculation_id}")
                
                return success
                
        except Exception as e:
            logger.error(f"Failed to update calculation failure {calculation_id}: {e}")
            return False
    
    def get_next_pending_calculation(self) -> Optional[Dict[str, Any]]:
        """次に実行すべき計算を取得"""
        try:
            with self.job_repo as repo:
                job = repo.get_next_job()
                if job:
                    # 対応する計算情報を取得
                    with CalculationRepository(repo.session) as calc_repo:
                        calculation = calc_repo.get_by_id(str(job.calculation_id))
                        if calculation:
                            return {
                                'job_id': job.id,
                                'calculation_id': calculation.id,
                                'instance_id': calculation.instance_id,
                                'method': calculation.method,
                                'basis_set': calculation.basis_set,
                                'parameters': calculation.parameters or {},
                                'priority': job.priority
                            }
                return None
                
        except Exception as e:
            logger.error(f"Failed to get next pending calculation: {e}")
            return None
    
    def get_calculation_results(self, calculation_id: str) -> Dict[str, Any]:
        """計算結果を取得"""
        try:
            with self.result_repo as repo:
                results = repo.get_by_calculation(calculation_id)
                
                result_dict = {}
                for result in results:
                    result_dict[result.result_type] = {
                        'data': result.result_data,
                        'file_path': result.file_path,
                        'created_at': result.created_at
                    }
                
                return result_dict
                
        except Exception as e:
            logger.error(f"Failed to get calculation results {calculation_id}: {e}")
            return {}
    
    def get_calculation_history(self, instance_id: Optional[str] = None,
                              limit: int = 100) -> List[Dict[str, Any]]:
        """計算履歴を取得"""
        try:
            with self.calculation_repo as repo:
                if instance_id:
                    calculations = repo.get_by_instance(instance_id)
                else:
                    # すべての計算を取得（制限付き）
                    calculations = repo.session.query(DBCalculation).order_by(
                        DBCalculation.start_time.desc()
                    ).limit(limit).all()
                
                history = []
                for calc in calculations:
                    calc_dict = self._convert_db_calculation_to_dict(calc)
                    
                    # インスタンス名と分子名を取得
                    with InstanceRepository(repo.session) as inst_repo:
                        instance = inst_repo.get_by_id(str(calc.instance_id))
                        if instance:
                            calc_dict['instance_name'] = instance.name
                            # 分子名を取得
                            from pyscf_front.services.molecule_service import MoleculeService
                            mol_service = MoleculeService()
                            molecule = mol_service.load_molecule_from_instance(str(calc.instance_id))
                            calc_dict['molecule_name'] = molecule.name if molecule else 'Unknown'
                        else:
                            calc_dict['instance_name'] = 'Unknown'
                            calc_dict['molecule_name'] = 'Unknown'
                    
                    # 結果の概要を追加
                    results = self.get_calculation_results(str(calc.id))
                    calc_dict['result_summary'] = {
                        'result_types': list(results.keys()),
                        'result_count': len(results)
                    }
                    
                    # 作成日時を文字列に変換
                    if calc.start_time:
                        calc_dict['created_at'] = calc.start_time.strftime('%Y-%m-%d %H:%M:%S')
                    else:
                        calc_dict['created_at'] = 'Unknown'
                    
                    # 完了日時を追加
                    if calc.end_time:
                        calc_dict['completed_at'] = calc.end_time.strftime('%Y-%m-%d %H:%M:%S')
                    
                    # エネルギー結果を取得
                    if 'total_energy' in results:
                        energy_data = results['total_energy']
                        if isinstance(energy_data, dict) and 'data' in energy_data:
                            calc_dict['total_energy'] = energy_data['data']
                        else:
                            calc_dict['total_energy'] = energy_data
                    
                    history.append(calc_dict)
                
                return history
                
        except Exception as e:
            logger.error(f"Failed to get calculation history: {e}")
            return []
    
    def convert_to_calculation_job(self, calculation_dict: Dict[str, Any]) -> CalculationJob:
        """データベースの計算情報をCalculationJobオブジェクトに変換"""
        from pyscf_front.services.molecule_service import MoleculeService
        
        # 分子情報を取得
        mol_service = MoleculeService()
        molecule = mol_service.load_molecule_from_instance(calculation_dict['instance_id'])
        
        if not molecule:
            raise ValueError(f"Cannot load molecule for instance {calculation_dict['instance_id']}")
        
        # CalculationJobオブジェクトを作成
        job = CalculationJob(
            id=calculation_dict['id'],
            molecule=molecule,
            method=calculation_dict['method'],
            basis_set=calculation_dict['basis_set'],
            charge=molecule.charge,
            multiplicity=molecule.multiplicity,
            parameters=calculation_dict.get('parameters', {}),
            status=self._convert_db_status_to_core_status(calculation_dict['status'])
        )
        
        return job
    
    def _convert_db_calculation_to_dict(self, calculation: DBCalculation) -> Dict[str, Any]:
        """データベース計算オブジェクトを辞書に変換"""
        return {
            'id': calculation.id,
            'instance_id': calculation.instance_id,
            'method': calculation.method,
            'basis_set': calculation.basis_set,
            'parameters': calculation.parameters or {},
            'status': calculation.status.value,
            'start_time': calculation.start_time,
            'end_time': calculation.end_time,
            'error_message': calculation.error_message,
            'max_iterations': calculation.max_iterations,
            'convergence_criteria': calculation.convergence_criteria or {}
        }
    
    def _convert_db_status_to_core_status(self, db_status: str) -> CoreCalculationStatus:
        """データベースステータスをコアステータスに変換"""
        status_mapping = {
            'pending': CoreCalculationStatus.PENDING,
            'running': CoreCalculationStatus.RUNNING,
            'completed': CoreCalculationStatus.COMPLETED,
            'failed': CoreCalculationStatus.FAILED,
            'cancelled': CoreCalculationStatus.CANCELLED
        }
        return status_mapping.get(db_status, CoreCalculationStatus.PENDING)