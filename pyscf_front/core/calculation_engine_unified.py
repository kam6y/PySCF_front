"""
統合された計算エンジン（プラグインシステム対応）
"""
import uuid
import time
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List
from enum import Enum

from PySide6.QtCore import QObject, Signal, QRunnable, QThreadPool
from loguru import logger

from pyscf_front.core.molecule import Molecule


class CalculationStatus(Enum):
    """計算ステータス"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class CalculationJob:
    """計算ジョブ"""
    id: str
    molecule: Molecule
    method: str = "B3LYP"
    basis_set: str = "6-31G"
    charge: int = 0
    multiplicity: int = 1
    parameters: Dict[str, Any] = field(default_factory=dict)
    status: CalculationStatus = CalculationStatus.PENDING
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    completed_at: Optional[float] = None
    results: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class CalculationSignals(QObject):
    """計算関連シグナル"""
    started = Signal(str)  # job_id
    progress = Signal(str, float, str)  # job_id, progress, message
    completed = Signal(str, dict)  # job_id, results
    failed = Signal(str, str)  # job_id, error_message


from pyscf_front.plugins.base import get_plugin_manager, PluginType
from pyscf_front.plugins.methods import PySCFMethodPlugin
from pyscf_front.plugins.basis_sets import PySCFBasisSetPlugin


class UnifiedCalculationWorker(QRunnable):
    """プラグインシステム対応の計算ワーカー"""
    
    def __init__(self, job: CalculationJob, signals: CalculationSignals, 
                 calculation_service = None,
                 calculation_id: Optional[str] = None):
        super().__init__()
        self.job = job
        self.signals = signals
        self.calculation_service = calculation_service
        self.calculation_id = calculation_id
        
        self.plugin_manager = get_plugin_manager()
    
    def run(self):
        """計算実行"""
        try:
            # 計算開始をデータベースに記録
            if self.calculation_service and self.calculation_id:
                self.calculation_service.start_calculation(self.calculation_id, "unified_worker")
            
            self.job.status = CalculationStatus.RUNNING
            self.job.started_at = time.time()
            self.signals.started.emit(self.job.id)
            
            # 手法プラグインの取得
            self.signals.progress.emit(self.job.id, 0.05, "計算手法を準備中...")
            method_plugin = self.plugin_manager.find_method_plugin(self.job.method)
            if method_plugin is None:
                raise RuntimeError(f"No plugin found for method: {self.job.method}")
            
            # 基底関数プラグインの取得
            self.signals.progress.emit(self.job.id, 0.1, "基底関数を準備中...")
            basis_plugin = self.plugin_manager.find_basis_plugin(self.job.basis_set)
            if basis_plugin is None:
                raise RuntimeError(f"No plugin found for basis set: {self.job.basis_set}")
            
            # 分子構築
            self.signals.progress.emit(self.job.id, 0.15, "分子構造を構築中...")
            mol_obj = self._build_pyscf_molecule(basis_plugin)
            
            # 計算オブジェクト作成
            self.signals.progress.emit(self.job.id, 0.2, "計算オブジェクトを作成中...")
            calculator = method_plugin.create_calculator(
                self.job.method,
                mol=mol_obj,
                **self.job.parameters
            )
            
            # 計算実行
            self.signals.progress.emit(self.job.id, 0.3, "量子化学計算を実行中...")
            energy = calculator.kernel()
            
            # 結果処理
            self.signals.progress.emit(self.job.id, 0.9, "結果を処理中...")
            results = self._process_results(calculator, energy)
            
            self.job.results = results
            self.job.status = CalculationStatus.COMPLETED
            self.job.completed_at = time.time()
            
            # 結果をデータベースに保存
            if self.calculation_service and self.calculation_id:
                self.calculation_service.complete_calculation(self.calculation_id, results)
            
            self.signals.progress.emit(self.job.id, 1.0, "計算完了")
            self.signals.completed.emit(self.job.id, results)
            
        except Exception as e:
            error_msg = str(e)
            logger.error(f"Calculation failed for job {self.job.id}: {error_msg}")
            self.job.status = CalculationStatus.FAILED
            self.job.error = error_msg
            
            # 失敗をデータベースに記録
            if self.calculation_service and self.calculation_id:
                self.calculation_service.fail_calculation(self.calculation_id, error_msg)
            
            self.signals.failed.emit(self.job.id, error_msg)
    
    def _build_pyscf_molecule(self, basis_plugin):
        """PySCF分子オブジェクトを構築"""
        try:
            import pyscf
            from pyscf import gto
        except ImportError:
            raise RuntimeError("PySCF is not available")
        
        # 原子座標の構築
        atom_string = ""
        atom_symbols = []
        for atom in self.job.molecule.atoms:
            atom_string += f"{atom.symbol} {atom.x} {atom.y} {atom.z}; "
            atom_symbols.append(atom.symbol)
        
        # 基底関数の妥当性チェック
        if not basis_plugin.validate_basis_for_atoms(self.job.basis_set, atom_symbols):
            raise ValueError(f"Basis set {self.job.basis_set} is not suitable for atoms {atom_symbols}")
        
        # 分子オブジェクト作成
        mol = gto.Mole()
        mol.atom = atom_string.rstrip('; ')
        mol.basis = self.job.basis_set
        mol.charge = self.job.charge
        mol.spin = self.job.multiplicity - 1
        mol.verbose = self.job.parameters.get('verbose', 4)
        mol.max_memory = self.job.parameters.get('max_memory', 4000)
        
        mol.build()
        
        return mol
    
    def _process_results(self, calculator, energy):
        """計算結果を処理"""
        results = {
            'total_energy': float(energy),
            'converged': bool(calculator.converged),
            'calculation_time': time.time() - (self.job.started_at or time.time()),
            'method': self.job.method,
            'basis_set': self.job.basis_set,
            'molecule_name': self.job.molecule.name
        }
        
        # 追加の結果情報
        try:
            # 軌道エネルギー
            if hasattr(calculator, 'mo_energy') and calculator.mo_energy is not None:
                mo_energies = calculator.mo_energy.tolist()
                results['mo_energies'] = mo_energies
                
                # HOMO/LUMO
                nelectron = calculator.mol.nelectron
                if nelectron > 0:
                    homo_idx = nelectron // 2 - 1
                    if homo_idx >= 0 and homo_idx < len(mo_energies):
                        results['homo_energy'] = mo_energies[homo_idx]
                    
                    lumo_idx = nelectron // 2
                    if lumo_idx < len(mo_energies):
                        results['lumo_energy'] = mo_energies[lumo_idx]
                        
                        if 'homo_energy' in results:
                            results['homo_lumo_gap'] = results['lumo_energy'] - results['homo_energy']
            
            # 双極子モーメント
            if hasattr(calculator, 'dip_moment'):
                try:
                    dipole = calculator.dip_moment()
                    if dipole is not None:
                        results['dipole_moment'] = float((dipole[0]**2 + dipole[1]**2 + dipole[2]**2)**0.5)
                        results['dipole_components'] = [float(x) for x in dipole]
                except:
                    pass
            
            # Mulliken電荷
            if hasattr(calculator, 'mulliken_pop'):
                try:
                    pop_and_charges = calculator.mulliken_pop()
                    if len(pop_and_charges) > 1:
                        charges = pop_and_charges[1]
                        if charges is not None:
                            results['atomic_charges'] = [float(x) for x in charges]
                except:
                    pass
            
            # 反復回数
            if hasattr(calculator, 'niter'):
                results['iterations'] = int(calculator.niter)
            
        except Exception as e:
            logger.warning(f"Failed to extract additional results: {e}")
        
        return results


class UnifiedCalculationEngine:
    """統合計算エンジン（プラグインシステム対応）"""
    
    def __init__(self, use_database: bool = True):
        self.use_database = use_database
        self.jobs: Dict[str, CalculationJob] = {}
        self.signals = CalculationSignals()
        self.thread_pool = QThreadPool()
        
        # サービス層（データベース使用時のみ）
        if self.use_database:
            try:
                from pyscf_front.services import CalculationService, MoleculeService, InstanceService
                self.calculation_service = CalculationService()
                self.molecule_service = MoleculeService()
                self.instance_service = InstanceService()
                
                # データベース接続テスト
                self._test_database_connection()
                
                # マッピング辞書
                self.job_to_calculation_mapping: Dict[str, str] = {}
                self.job_to_instance_mapping: Dict[str, str] = {}
                
                logger.info("Database services initialized successfully")
                
            except Exception as e:
                logger.warning(f"Database initialization failed: {e}")
                logger.info("Falling back to memory-only mode")
                self.use_database = False
                self.calculation_service = None
                self.molecule_service = None
                self.instance_service = None
        else:
            self.calculation_service = None
            self.molecule_service = None
            self.instance_service = None
        
        # プラグインマネージャーの初期化
        self.plugin_manager = get_plugin_manager()
        self._initialize_plugins()
        
        logger.info(f"UnifiedCalculationEngine initialized (database: {self.use_database})")
    
    def _test_database_connection(self):
        """データベース接続をテスト"""
        try:
            if self.instance_service:
                # 実際のデータベース接続テスト
                self.instance_service.get_all_instances(include_details=False)
                logger.info("Database connection test successful")
        except Exception as e:
            logger.error(f"Database connection test failed: {e}")
            raise
    
    def _initialize_plugins(self):
        """プラグインの初期化"""
        try:
            # PySCF手法プラグインの登録
            pyscf_methods = PySCFMethodPlugin()
            if self.plugin_manager.register_plugin(pyscf_methods):
                self.plugin_manager.initialize_plugin(PluginType.METHOD, "PySCFMethods")
                logger.info("PySCF methods plugin registered and initialized")
            
            # PySCF基底関数プラグインの登録
            pyscf_basis = PySCFBasisSetPlugin()
            if self.plugin_manager.register_plugin(pyscf_basis):
                self.plugin_manager.initialize_plugin(PluginType.BASIS_SET, "PySCFBasisSets")
                logger.info("PySCF basis sets plugin registered and initialized")
                
        except Exception as e:
            logger.error(f"Failed to initialize plugins: {e}")
    
    def submit_calculation(
        self,
        molecule: Molecule,
        method: str = "B3LYP",
        basis_set: str = "6-31G(d)",
        charge: Optional[int] = None,
        multiplicity: Optional[int] = None,
        parameters: Optional[Dict[str, Any]] = None,
        priority: int = 5
    ) -> str:
        """計算を投入（メモリのみ）"""
        
        # プラグインの可用性チェック
        if not self.plugin_manager.find_method_plugin(method):
            raise ValueError(f"Method {method} is not supported by any plugin")
        
        if not self.plugin_manager.find_basis_plugin(basis_set):
            raise ValueError(f"Basis set {basis_set} is not supported by any plugin")
        
        # 電荷とスピン多重度の設定
        if charge is None:
            charge = molecule.charge
        if multiplicity is None:
            multiplicity = molecule.multiplicity
        
        # ジョブID生成
        job_id = str(uuid.uuid4())
        
        # ジョブ作成
        job = CalculationJob(
            id=job_id,
            molecule=molecule,
            method=method,
            basis_set=basis_set,
            charge=charge,
            multiplicity=multiplicity,
            parameters=parameters or {}
        )
        
        self.jobs[job_id] = job
        
        # ワーカー作成と投入
        worker = UnifiedCalculationWorker(job, self.signals)
        self.thread_pool.start(worker)
        
        logger.info(f"Submitted calculation job: {job_id}")
        return job_id
    
    def submit_calculation_with_persistence(
        self,
        molecule: Molecule,
        method: str = "B3LYP",
        basis_set: str = "6-31G(d)",
        charge: Optional[int] = None,
        multiplicity: Optional[int] = None,
        parameters: Optional[Dict[str, Any]] = None,
        instance_name: Optional[str] = None,
        priority: int = 5
    ) -> str:
        """計算を投入（データベース永続化付き）"""
        
        if not self.use_database:
            raise RuntimeError("Database persistence is not enabled")
        
        # プラグインの可用性チェック
        if not self.plugin_manager.find_method_plugin(method):
            raise ValueError(f"Method {method} is not supported by any plugin")
        
        if not self.plugin_manager.find_basis_plugin(basis_set):
            raise ValueError(f"Basis set {basis_set} is not supported by any plugin")
        
        # インスタンス作成（分子を含む）
        if instance_name is None:
            instance_name = f"{molecule.name}_{method}_{basis_set}_{uuid.uuid4().hex[:8]}"
        
        instance_id = self.instance_service.create_complete_instance(
            name=instance_name,
            description=f"Calculation: {method}/{basis_set}",
            molecule_data={'molecule_object': molecule},
            calculation_config={
                'method': method,
                'basis_set': basis_set,
                'parameters': parameters or {},
                'priority': priority
            }
        )
        
        # 作成された計算IDを取得
        calculations = self.calculation_service.get_calculations_by_instance(instance_id)
        if not calculations:
            raise RuntimeError("Failed to create calculation in database")
        
        calculation_id = calculations[0]['id']
        
        # 従来のジョブIDを生成
        job_id = str(uuid.uuid4())
        
        # 電荷とスピン多重度の設定
        if charge is None:
            charge = molecule.charge
        if multiplicity is None:
            multiplicity = molecule.multiplicity
        
        # ジョブ作成
        job = CalculationJob(
            id=job_id,
            molecule=molecule,
            method=method,
            basis_set=basis_set,
            charge=charge,
            multiplicity=multiplicity,
            parameters=parameters or {}
        )
        
        self.jobs[job_id] = job
        
        # マッピングを保存
        self.job_to_calculation_mapping[job_id] = calculation_id
        self.job_to_instance_mapping[job_id] = instance_id
        
        # 統合ワーカーを作成して投入
        worker = UnifiedCalculationWorker(
            job, self.signals, self.calculation_service, calculation_id
        )
        self.thread_pool.start(worker)
        
        logger.info(
            f"Submitted persistent calculation job: {job_id} "
            f"(instance: {instance_id}, calculation: {calculation_id})"
        )
        return job_id
    
    def get_available_methods(self) -> List[str]:
        """利用可能な計算手法の一覧"""
        return self.plugin_manager.get_available_methods()
    
    def get_available_basis_sets(self) -> List[str]:
        """利用可能な基底関数の一覧"""
        return self.plugin_manager.get_available_basis_sets()
    
    def get_job_status(self, job_id: str) -> Optional[CalculationJob]:
        """ジョブステータスを取得"""
        return self.jobs.get(job_id)
    
    def cancel_job(self, job_id: str) -> bool:
        """ジョブをキャンセル"""
        if job_id in self.jobs:
            self.jobs[job_id].status = CalculationStatus.CANCELLED
            return True
        return False
    
    def list_jobs(self) -> List[Dict[str, Any]]:
        """ジョブリストを取得"""
        job_list = []
        for job in self.jobs.values():
            job_info = {
                'id': job.id,
                'molecule_name': job.molecule.name,
                'method': job.method,
                'basis_set': job.basis_set,
                'status': job.status.value,
                'progress': 0.0,  # 詳細な進捗は別途実装
                'started_at': job.started_at,
                'completed_at': job.completed_at
            }
            
            if job.status == CalculationStatus.COMPLETED and job.results:
                job_info['total_energy'] = job.results.get('total_energy')
                job_info['converged'] = job.results.get('converged')
            
            job_list.append(job_info)
        
        return job_list
    
    def clear_completed_jobs(self):
        """完了したジョブをクリア"""
        completed_jobs = [
            job_id for job_id, job in self.jobs.items()
            if job.status in [CalculationStatus.COMPLETED, CalculationStatus.FAILED, CalculationStatus.CANCELLED]
        ]
        
        for job_id in completed_jobs:
            del self.jobs[job_id]
            # マッピングもクリア
            if self.use_database:
                self.job_to_calculation_mapping.pop(job_id, None)
                self.job_to_instance_mapping.pop(job_id, None)
        
        logger.info(f"Cleared {len(completed_jobs)} completed jobs")
    
    def estimate_calculation_time(self, molecule: Molecule, method: str, basis_set: str) -> str:
        """計算時間を推定"""
        method_plugin = self.plugin_manager.find_method_plugin(method)
        basis_plugin = self.plugin_manager.find_basis_plugin(basis_set)
        
        if not method_plugin or not basis_plugin:
            return "推定不可"
        
        try:
            atom_symbols = [atom.symbol for atom in molecule.atoms]
            basis_size = basis_plugin.get_basis_size_estimate(basis_set, atom_symbols)
            molecule_size = len(molecule.atoms)
            
            time_seconds = method_plugin.estimate_time(method, molecule_size, basis_size)
            
            if time_seconds < 60:
                return f"約{time_seconds:.0f}秒"
            elif time_seconds < 3600:
                return f"約{time_seconds/60:.0f}分"
            else:
                return f"約{time_seconds/3600:.1f}時間"
        
        except Exception as e:
            logger.error(f"Time estimation failed: {e}")
            return "推定不可"
    
    # データベース関連メソッド（use_database=Trueの場合のみ有効）
    def get_calculation_history(self, instance_id: Optional[str] = None) -> List[Dict[str, Any]]:
        """計算履歴を取得"""
        if not self.use_database:
            return []
        return self.calculation_service.get_calculation_history(instance_id)
    
    def get_instance_calculations(self, instance_id: str) -> List[Dict[str, Any]]:
        """インスタンスの計算一覧を取得"""
        if not self.use_database:
            return []
        return self.calculation_service.get_calculations_by_instance(instance_id)
    
    def list_all_instances(self) -> List[Dict[str, Any]]:
        """全インスタンスの一覧を取得"""
        if not self.use_database:
            return []
        return self.instance_service.get_all_instances(include_details=False)
    
    def get_instance_details(self, instance_id: str) -> Optional[Dict[str, Any]]:
        """インスタンスの詳細情報を取得"""
        if not self.use_database:
            return None
        return self.instance_service.get_instance_details(instance_id)
    
    def delete_instance(self, instance_id: str) -> bool:
        """インスタンスを削除"""
        if not self.use_database:
            return False
        return self.instance_service.delete_instance(instance_id)
    
    def resume_interrupted_calculations(self) -> List[str]:
        """中断された計算を再開"""
        if not self.use_database:
            return []
        
        # 実装は IntegratedCalculationEngine と同様
        # 詳細省略...
        return []
    
    def finalize(self):
        """エンジンの終了処理"""
        self.plugin_manager.finalize_all()
        logger.info("UnifiedCalculationEngine finalized")


# エイリアス（既存コードとの互換性のため）
CalculationEngine = UnifiedCalculationEngine