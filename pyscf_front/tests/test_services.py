"""
サービス層のユニットテスト
"""
import pytest

from pyscf_front.services import MoleculeService, CalculationService, InstanceService
from pyscf_front.database import (
    InstanceRepository, MoleculeRepository, CalculationRepository, ResultRepository, JobQueueRepository,
)


class TestMoleculeService:
    """MoleculeServiceのテスト"""
    
    def test_create_molecule_instance(self, test_session, sample_molecule):
        """分子インスタンスの作成テスト"""
        # Arrange
        service = MoleculeService()
        service.molecule_repo = MoleculeRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        # Act
        instance_id = service.create_molecule_instance(
            name="Test Water Molecule",
            molecule=sample_molecule,
            description="Test water molecule for unit testing"
        )
        
        # Assert
        assert instance_id is not None
        assert isinstance(instance_id, str)
        
        # データベースから確認
        with InstanceRepository(test_session) as repo:
            instance = repo.get_by_id(instance_id)
            assert instance is not None
            assert instance.name == "Test Water Molecule"  # type: ignore
            assert instance.description == "Test water molecule for unit testing"  # type: ignore
    
    def test_load_molecule_from_instance(self, test_session, sample_molecule):
        """インスタンスから分子の読み込みテスト"""
        # Arrange
        service = MoleculeService()
        service.molecule_repo = MoleculeRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        # 分子インスタンスを作成
        instance_id = service.create_molecule_instance(
            name="Test Molecule",
            molecule=sample_molecule
        )
        
        # Act
        loaded_molecule = service.load_molecule_from_instance(instance_id)
        
        # Assert
        assert loaded_molecule is not None
        assert loaded_molecule.name == sample_molecule.name
        assert loaded_molecule.charge == sample_molecule.charge
        assert loaded_molecule.multiplicity == sample_molecule.multiplicity
        assert len(loaded_molecule.atoms) == len(sample_molecule.atoms)
        
        # 原子データの確認
        for original, loaded in zip(sample_molecule.atoms, loaded_molecule.atoms):
            assert original.symbol == loaded.symbol
            assert abs(original.x - loaded.x) < 1e-6
            assert abs(original.y - loaded.y) < 1e-6
            assert abs(original.z - loaded.z) < 1e-6
    
    def test_get_molecule_instances(self, test_session, sample_molecule, sample_molecule_methane):
        """分子インスタンス一覧の取得テスト"""
        # Arrange
        service = MoleculeService()
        service.molecule_repo = MoleculeRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        # 複数の分子インスタンスを作成
        instance1_id = service.create_molecule_instance("Water", sample_molecule)
        instance2_id = service.create_molecule_instance("Methane", sample_molecule_methane)
        
        # Act
        instances = service.get_molecule_instances()
        
        # Assert
        assert len(instances) == 2
        
        # インスタンス情報の確認
        instance_ids = [inst['instance_id'] for inst in instances]
        assert instance1_id in instance_ids
        assert instance2_id in instance_ids
        
        # 分子情報の確認
        water_instance = next(inst for inst in instances if inst['name'] == 'Water')
        assert water_instance['molecule']['name'] == 'H2O'
        
        methane_instance = next(inst for inst in instances if inst['name'] == 'Methane')
        assert methane_instance['molecule']['name'] == 'CH4'
    
    def test_delete_molecule_instance(self, test_session, sample_molecule):
        """分子インスタンスの削除テスト"""
        # Arrange
        service = MoleculeService()
        service.molecule_repo = MoleculeRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        instance_id = service.create_molecule_instance("Test Molecule", sample_molecule)
        
        # Act
        success = service.delete_molecule_instance(instance_id)
        
        # Assert
        assert success is True
        
        # 削除されたことを確認
        with InstanceRepository(test_session) as repo:
            instance = repo.get_by_id(instance_id)
            assert instance is None


class TestCalculationService:
    """CalculationServiceのテスト"""
    
    def test_create_calculation(self, test_session):
        """計算の作成テスト"""
        # Arrange
        service = CalculationService()
        service.calculation_repo = CalculationRepository(test_session)
        service.job_repo = JobQueueRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        # テスト用インスタンスを作成
        with InstanceRepository(test_session) as repo:
            instance = repo.create("Test Instance")
        
        # Act
        calculation_id = service.create_calculation(
            instance_id=str(instance.id),  # type: ignore
            method="B3LYP",
            basis_set="6-31G(d)",
            parameters={'max_iter': 100},
            priority=8
        )
        
        # Assert
        assert calculation_id is not None
        calc_info = service.get_calculation_by_id(calculation_id)
        assert calc_info is not None
        assert calc_info['method'] == "B3LYP"
        assert calc_info['basis_set'] == "6-31G(d)"
        assert calc_info['parameters']['max_iter'] == 100
        assert calc_info['status'] == 'pending'
    
    def test_calculation_lifecycle(self, test_session):
        """計算のライフサイクルテスト"""
        # Arrange
        service = CalculationService()
        service.calculation_repo = CalculationRepository(test_session)
        service.job_repo = JobQueueRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        with InstanceRepository(test_session) as repo:
            instance = repo.create("Test Instance")
        
        calculation_id = service.create_calculation(
            instance_id=str(instance.id),  # type: ignore
            method="HF",
            basis_set="STO-3G"
        )
        
        # Act & Assert - 計算開始
        success = service.start_calculation(calculation_id, "worker_01")
        assert success is True
        
        calc_info = service.get_calculation_by_id(calculation_id)
        assert calc_info is not None
        assert calc_info['status'] == 'running'
        
        # Act & Assert - 計算完了
        results = {
            'total_energy': -1.123456,
            'converged': True,
            'calculation_time': 45.2
        }
        success = service.complete_calculation(calculation_id, results)
        assert success is True
        
        calc_info = service.get_calculation_by_id(calculation_id)
        assert calc_info is not None
        assert calc_info['status'] == 'completed'
        
        # 結果確認
        service.result_repo = ResultRepository(test_session)
        saved_results = service.get_calculation_results(calculation_id)
        assert 'total_energy' in saved_results
        assert saved_results['total_energy']['data'] == -1.123456
    
    def test_get_next_pending_calculation(self, test_session):
        """次の待機中計算の取得テスト"""
        # Arrange
        service = CalculationService()
        service.calculation_repo = CalculationRepository(test_session)
        service.job_repo = JobQueueRepository(test_session)
        service.instance_repo = InstanceRepository(test_session)
        
        with InstanceRepository(test_session) as repo:
            instance = repo.create("Test Instance")
        
        # 優先度の異なる計算を作成
        calc1_id = service.create_calculation(
            instance_id=str(instance.id),  # type: ignore
            method="HF",
            basis_set="STO-3G",
            priority=3
        )
        
        calc2_id = service.create_calculation(
            instance_id=str(instance.id),  # type: ignore
            method="B3LYP",
            basis_set="6-31G",
            priority=8  # より高い優先度
        )
        
        # Act
        next_calc = service.get_next_pending_calculation()
        
        # Assert - 優先度の高い計算が返されること
        assert next_calc is not None
        assert next_calc['calculation_id'] == calc2_id
        assert next_calc['priority'] == 8


class TestInstanceService:
    """InstanceServiceのテスト"""
    
    def test_create_complete_instance(self, test_session, sample_molecule):
        """完全なインスタンスの作成テスト"""
        # Arrange
        service = InstanceService()
        service.instance_repo = InstanceRepository(test_session)
        service.molecule_service.instance_repo = InstanceRepository(test_session)
        service.molecule_service.molecule_repo = MoleculeRepository(test_session)
        service.calculation_service.calculation_repo = CalculationRepository(test_session)
        service.calculation_service.job_repo = JobQueueRepository(test_session)
        service.calculation_service.instance_repo = InstanceRepository(test_session)
        
        # Act
        instance_id = service.create_complete_instance(
            name="Complete Water Analysis",
            description="Complete analysis of water molecule",
            molecule_data={
                'molecule_object': sample_molecule
            },
            calculation_config={
                'method': 'B3LYP',
                'basis_set': '6-31G(d)',
                'priority': 7
            }
        )
        
        # Assert
        assert instance_id is not None
        
        # インスタンス詳細を確認
        details = service.get_instance_details(instance_id)
        assert details is not None
        assert details['name'] == "Complete Water Analysis"
        assert details['molecule']['name'] == 'H2O'
        assert len(details['calculations']) == 1
        assert details['calculations'][0]['method'] == 'B3LYP'
        assert details['calculations'][0]['basis_set'] == '6-31G(d)'
    
    def test_get_instance_details(self, test_session, sample_molecule):
        """インスタンス詳細の取得テスト"""
        # Arrange
        service = InstanceService()
        service.instance_repo = InstanceRepository(test_session)
        service.molecule_service.instance_repo = InstanceRepository(test_session)
        service.molecule_service.molecule_repo = MoleculeRepository(test_session)
        service.calculation_service.calculation_repo = CalculationRepository(test_session)
        
        instance_id = service.create_complete_instance(
            name="Test Instance",
            molecule_data={'molecule_object': sample_molecule}
        )
        
        # Act
        details = service.get_instance_details(instance_id)
        
        # Assert
        assert details is not None
        assert details['instance_id'] == instance_id
        assert details['name'] == "Test Instance"
        assert 'molecule' in details
        assert 'calculations' in details
        assert 'statistics' in details
        
        # 分子情報の確認
        molecule_info = details['molecule']
        assert molecule_info['name'] == 'H2O'
        assert molecule_info['charge'] == 0
        assert molecule_info['multiplicity'] == 1
        assert molecule_info['atom_count'] == 3
    
    def test_search_instances(self, test_session, sample_molecule, sample_molecule_methane):
        """インスタンス検索テスト"""
        # Arrange
        service = InstanceService()
        service.instance_repo = InstanceRepository(test_session)
        service.molecule_service.instance_repo = InstanceRepository(test_session)
        service.molecule_service.molecule_repo = MoleculeRepository(test_session)
        
        # 複数のインスタンスを作成
        water_id = service.create_complete_instance(
            name="Water Analysis",
            description="Analysis of water molecule",
            molecule_data={'molecule_object': sample_molecule}
        )
        
        methane_id = service.create_complete_instance(
            name="Methane Study",
            description="Study of methane molecule",
            molecule_data={'molecule_object': sample_molecule_methane}
        )
        
        # Act & Assert - 名前で検索
        results = service.search_instances("water")
        assert len(results) == 1
        assert results[0]['instance_id'] == water_id
        
        # Act & Assert - 説明で検索
        results = service.search_instances("study")
        assert len(results) == 1
        assert results[0]['instance_id'] == methane_id
        
        # Act & Assert - 存在しない検索語
        results = service.search_instances("benzene")
        assert len(results) == 0