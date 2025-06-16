"""
統合計算エンジンのテスト
"""
import pytest
from unittest.mock import Mock, patch, MagicMock
import uuid

from pyscf_front.core.calculation_engine_unified import UnifiedCalculationEngine, UnifiedCalculationWorker
from pyscf_front.core.calculation_engine import CalculationJob, CalculationStatus, CalculationSignals
from pyscf_front.plugins.base import PluginManager, PluginType


class TestUnifiedCalculationEngine:
    """UnifiedCalculationEngine のテスト"""
    
    @pytest.fixture
    def mock_plugins(self):
        """プラグインのモック設定"""
        with patch('pyscf_front.core.calculation_engine_unified.get_plugin_manager') as mock_manager_func:
            mock_manager = Mock(spec=PluginManager)
            mock_manager_func.return_value = mock_manager
            
            # 手法プラグインのモック
            mock_method_plugin = Mock()
            mock_method_plugin.get_supported_methods.return_value = ['HF', 'B3LYP', 'PBE']
            mock_method_plugin.create_calculator.return_value = Mock()
            mock_method_plugin.estimate_time.return_value = 120.0
            
            # 基底関数プラグインのモック
            mock_basis_plugin = Mock()
            mock_basis_plugin.get_supported_basis_sets.return_value = ['STO-3G', '6-31G', '6-31G(d)']
            mock_basis_plugin.validate_basis_for_atoms.return_value = True
            mock_basis_plugin.get_basis_size_estimate.return_value = 50
            
            # プラグインマネージャーの設定
            mock_manager.get_available_methods.return_value = ['HF', 'B3LYP', 'PBE']
            mock_manager.get_available_basis_sets.return_value = ['STO-3G', '6-31G', '6-31G(d)']
            mock_manager.find_method_plugin.return_value = mock_method_plugin
            mock_manager.find_basis_plugin.return_value = mock_basis_plugin
            mock_manager.register_plugin.return_value = True
            mock_manager.initialize_plugin.return_value = True
            
            yield {
                'manager': mock_manager,
                'method_plugin': mock_method_plugin,
                'basis_plugin': mock_basis_plugin
            }
    
    def test_engine_initialization_database_required(self, mock_plugins):
        """データベース必須でのエンジン初期化テスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            engine = UnifiedCalculationEngine()
            
            assert engine.calculation_service is not None
            assert engine.molecule_service is not None
            assert engine.instance_service is not None
            assert isinstance(engine.jobs, dict)
            assert len(engine.jobs) == 0
    
    def test_database_initialization_failure(self, mock_plugins):
        """データベース初期化失敗テスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService', side_effect=Exception("DB Connection Failed")):
            
            # データベース初期化の失敗で RuntimeError が発生することを確認
            with pytest.raises(RuntimeError, match="Failed to initialize SQLite database"):
                UnifiedCalculationEngine()
    
    def test_available_methods_and_basis_sets(self, mock_plugins):
        """利用可能な手法と基底関数の取得テスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            engine = UnifiedCalculationEngine()
            
            methods = engine.get_available_methods()
            basis_sets = engine.get_available_basis_sets()
            
            assert methods == ['HF', 'B3LYP', 'PBE']
            assert basis_sets == ['STO-3G', '6-31G', '6-31G(d)']
            
            # プラグインマネージャーが呼ばれたことを確認
            mock_plugins['manager'].get_available_methods.assert_called()
            mock_plugins['manager'].get_available_basis_sets.assert_called()
    
    def test_calculation_submission_with_persistence(self, mock_plugins, sample_molecule):
        """永続化計算投入テスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.instance_service.InstanceService') as mock_inst_svc, \
             patch('pyscf_front.services.calculation_service.CalculationService') as mock_calc_svc, \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # サービスのモック設定
            mock_inst_service = Mock()
            mock_inst_service.create_complete_instance.return_value = "test-instance-id"
            mock_inst_svc.return_value = mock_inst_service
            
            mock_calc_service = Mock()
            mock_calc_service.get_calculations_by_instance.return_value = [
                {'id': 'test-calc-id', 'status': 'pending'}
            ]
            mock_calc_svc.return_value = mock_calc_service
            
            engine = UnifiedCalculationEngine()
            
            # 計算投入
            job_id = engine.submit_calculation(
                molecule=sample_molecule,
                method="HF",
                basis_set="STO-3G",
                instance_name="Test Calculation"
            )
            
            # 検証
            assert job_id is not None
            assert isinstance(job_id, str)
            assert job_id in engine.jobs
            assert job_id in engine.job_to_calculation_mapping
            assert job_id in engine.job_to_instance_mapping
            
            job = engine.jobs[job_id]
            assert job.molecule == sample_molecule
            assert job.method == "HF"
            assert job.basis_set == "STO-3G"
            assert job.status == CalculationStatus.PENDING
            
            # プラグインが呼ばれたことを確認
            mock_plugins['manager'].find_method_plugin.assert_called_with("HF")
            mock_plugins['manager'].find_basis_plugin.assert_called_with("STO-3G")
            
            # サービスが呼ばれたことを確認
            mock_inst_service.create_complete_instance.assert_called_once()
            mock_calc_service.get_calculations_by_instance.assert_called_once()
    
    
    def test_unsupported_method_handling(self, mock_plugins, sample_molecule):
        """サポートされていない手法のハンドリングテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # 手法プラグインが見つからない場合
            mock_plugins['manager'].find_method_plugin.return_value = None
            
            engine = UnifiedCalculationEngine()
            
            with pytest.raises(ValueError, match="Method INVALID_METHOD is not supported"):
                engine.submit_calculation(
                    molecule=sample_molecule,
                    method="INVALID_METHOD",
                    basis_set="STO-3G"
                )
    
    def test_unsupported_basis_set_handling(self, mock_plugins, sample_molecule):
        """サポートされていない基底関数のハンドリングテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # 基底関数プラグインが見つからない場合
            mock_plugins['manager'].find_basis_plugin.return_value = None
            
            engine = UnifiedCalculationEngine()
            
            with pytest.raises(ValueError, match="Basis set INVALID_BASIS is not supported"):
                engine.submit_calculation(
                    molecule=sample_molecule,
                    method="HF",
                    basis_set="INVALID_BASIS"
                )
    
    def test_job_status_and_management(self, mock_plugins, sample_molecule):
        """ジョブステータスと管理のテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.instance_service.InstanceService') as mock_inst_svc, \
             patch('pyscf_front.services.calculation_service.CalculationService') as mock_calc_svc, \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # サービスのモック設定
            mock_inst_service = Mock()
            mock_inst_service.create_complete_instance.return_value = "test-instance-id"
            mock_inst_svc.return_value = mock_inst_service
            
            mock_calc_service = Mock()
            mock_calc_service.get_calculations_by_instance.return_value = [
                {'id': 'test-calc-id', 'status': 'pending'}
            ]
            mock_calc_svc.return_value = mock_calc_service
            
            engine = UnifiedCalculationEngine()
            
            # 計算投入
            job_id = engine.submit_calculation(
                molecule=sample_molecule,
                method="HF",
                basis_set="STO-3G"
            )
            
            # ステータス取得
            job = engine.get_job_status(job_id)
            assert job is not None
            assert job.status == CalculationStatus.PENDING
            
            # ジョブリスト取得
            job_list = engine.list_jobs()
            assert len(job_list) == 1
            assert job_list[0]['id'] == job_id
            assert job_list[0]['method'] == "HF"
            
            # ジョブキャンセル
            success = engine.cancel_job(job_id)
            assert success is True
            assert engine.jobs[job_id].status == CalculationStatus.CANCELLED
    
    def test_time_estimation(self, mock_plugins, sample_molecule):
        """計算時間推定のテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            engine = UnifiedCalculationEngine()
            
            time_str = engine.estimate_calculation_time(
                molecule=sample_molecule,
                method="HF",
                basis_set="STO-3G"
            )
            
            assert isinstance(time_str, str)
            assert any(unit in time_str for unit in ["秒", "分", "時間"])
            
            # プラグインメソッドが呼ばれたことを確認
            mock_plugins['method_plugin'].estimate_time.assert_called()
            mock_plugins['basis_plugin'].get_basis_size_estimate.assert_called()
    
    def test_completed_jobs_cleanup(self, mock_plugins, sample_molecule):
        """完了ジョブのクリーンアップテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.instance_service.InstanceService') as mock_inst_svc, \
             patch('pyscf_front.services.calculation_service.CalculationService') as mock_calc_svc, \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # サービスのモック設定
            mock_inst_service = Mock()
            mock_inst_service.create_complete_instance.return_value = "test-instance-id"
            mock_inst_svc.return_value = mock_inst_service
            
            mock_calc_service = Mock()
            mock_calc_service.get_calculations_by_instance.return_value = [
                {'id': 'test-calc-id', 'status': 'pending'}
            ]
            mock_calc_svc.return_value = mock_calc_service
            
            engine = UnifiedCalculationEngine()
            
            # 複数のジョブを作成
            job_ids = []
            for i in range(3):
                job_id = engine.submit_calculation(
                    molecule=sample_molecule,
                    method="HF",
                    basis_set="STO-3G"
                )
                job_ids.append(job_id)
            
            # いくつかのジョブを完了状態に
            engine.jobs[job_ids[0]].status = CalculationStatus.COMPLETED
            engine.jobs[job_ids[1]].status = CalculationStatus.FAILED
            
            # クリーンアップ実行
            engine.clear_completed_jobs()
            
            # 完了・失敗したジョブが削除され、実行中のジョブは残っている
            assert job_ids[0] not in engine.jobs  # 完了済み削除
            assert job_ids[1] not in engine.jobs  # 失敗済み削除
            assert job_ids[2] in engine.jobs      # 実行中は残存


class TestUnifiedCalculationWorker:
    """UnifiedCalculationWorker のテスト"""
    
    @pytest.fixture
    def mock_worker_environment(self, sample_molecule):
        """ワーカーテスト環境のモック"""
        job = CalculationJob(
            id="test-worker-job",
            molecule=sample_molecule,
            method="HF",
            basis_set="STO-3G"
        )
        
        signals = CalculationSignals()
        
        # プラグインマネージャーのモック
        mock_plugin_manager = Mock()
        
        mock_method_plugin = Mock()
        mock_basis_plugin = Mock()
        
        mock_plugin_manager.find_method_plugin.return_value = mock_method_plugin
        mock_plugin_manager.find_basis_plugin.return_value = mock_basis_plugin
        
        return {
            'job': job,
            'signals': signals,
            'plugin_manager': mock_plugin_manager,
            'method_plugin': mock_method_plugin,
            'basis_plugin': mock_basis_plugin
        }
    
    def test_worker_initialization(self, mock_worker_environment):
        """ワーカー初期化テスト"""
        env = mock_worker_environment
        
        worker = UnifiedCalculationWorker(env['job'], env['signals'])
        
        assert worker.job == env['job']
        assert worker.signals == env['signals']
        assert worker.calculation_service is None
        assert worker.calculation_id is None
    
    def test_worker_with_database_service(self, mock_worker_environment):
        """データベースサービス付きワーカーテスト"""
        env = mock_worker_environment
        mock_calc_service = Mock()
        
        worker = UnifiedCalculationWorker(
            env['job'], 
            env['signals'],
            calculation_service=mock_calc_service,
            calculation_id="test-calc-id"
        )
        
        assert worker.calculation_service == mock_calc_service
        assert worker.calculation_id == "test-calc-id"
    
    @patch('pyscf_front.core.calculation_engine_unified.get_plugin_manager')
    def test_worker_plugin_not_found_error(self, mock_get_manager, mock_worker_environment):
        """プラグインが見つからない場合のエラーテスト"""
        env = mock_worker_environment
        
        # プラグインが見つからない設定
        mock_manager = Mock()
        mock_manager.find_method_plugin.return_value = None
        mock_get_manager.return_value = mock_manager
        
        worker = UnifiedCalculationWorker(env['job'], env['signals'])
        
        # シグナルのモック
        error_signal_mock = Mock()
        env['signals'].failed.connect(error_signal_mock)
        
        worker.run()
        
        # エラーシグナルが発信されることを確認
        error_signal_mock.assert_called_once()
        args = error_signal_mock.call_args[0]
        assert "No plugin found for method" in args[1]


class TestUnifiedEngineIntegration:
    """統合エンジンの統合テスト"""
    
    def test_end_to_end_calculation_workflow(self, sample_molecule, mock_pyscf):
        """エンドツーエンドの計算ワークフローテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin') as mock_method_cls, \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin') as mock_basis_cls, \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # プラグインのモック設定
            mock_method_plugin = Mock()
            mock_method_plugin.get_supported_methods.return_value = ['HF']
            mock_method_plugin.create_calculator.return_value = mock_pyscf['mf']
            mock_method_plugin.estimate_time.return_value = 60.0
            mock_method_cls.return_value = mock_method_plugin
            
            mock_basis_plugin = Mock()
            mock_basis_plugin.get_supported_basis_sets.return_value = ['STO-3G']
            mock_basis_plugin.validate_basis_for_atoms.return_value = True
            mock_basis_plugin.get_basis_size_estimate.return_value = 20
            mock_basis_cls.return_value = mock_basis_plugin
            
            # エンジン作成
            engine = UnifiedCalculationEngine()
            
            # 利用可能な手法と基底関数の確認
            methods = engine.get_available_methods()
            basis_sets = engine.get_available_basis_sets()
            
            assert 'HF' in methods
            assert 'STO-3G' in basis_sets
            
            # 計算投入
            job_id = engine.submit_calculation(
                molecule=sample_molecule,
                method="HF",
                basis_set="STO-3G"
            )
            
            # ジョブが作成されたことを確認
            assert job_id in engine.jobs
            job = engine.jobs[job_id]
            assert job.status == CalculationStatus.PENDING
            
            # 時間推定
            time_estimate = engine.estimate_calculation_time(
                sample_molecule, "HF", "STO-3G"
            )
            assert "分" in time_estimate or "秒" in time_estimate
    
    def test_database_only_mode(self, sample_molecule):
        """データベース専用モードテスト"""
        with patch('pyscf_front.core.calculation_engine_unified.PySCFMethodPlugin'), \
             patch('pyscf_front.core.calculation_engine_unified.PySCFBasisSetPlugin'), \
             patch('pyscf_front.services.calculation_service.CalculationService'), \
             patch('pyscf_front.services.molecule_service.MoleculeService'), \
             patch('pyscf_front.services.instance_service.InstanceService'), \
             patch.object(UnifiedCalculationEngine, '_test_database_connection'):
            
            # データベース専用エンジン
            db_engine = UnifiedCalculationEngine()
            
            # APIが利用可能
            methods = db_engine.get_available_methods()
            basis_sets = db_engine.get_available_basis_sets()
            
            assert isinstance(methods, list)
            assert isinstance(basis_sets, list)
            
            # 基本的な計算投入機能
            assert hasattr(db_engine, 'submit_calculation')
            assert hasattr(db_engine, 'get_job_status')
            assert hasattr(db_engine, 'list_jobs')
            
            # データベース機能
            assert hasattr(db_engine, 'get_calculation_history')
            assert hasattr(db_engine, 'list_all_instances')
            assert hasattr(db_engine, 'get_instance_details')


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])