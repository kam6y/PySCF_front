"""
プラグインシステムのテスト
"""
import pytest
from unittest.mock import Mock, patch
from typing import Dict, Any, List

from pyscf_front.plugins.base import (
    PluginManager, PluginInterface, CalculationMethodPlugin, BasisSetPlugin,
    PluginInfo, PluginType, get_plugin_manager, register_plugin
)
from pyscf_front.plugins.methods.pyscf_methods import PySCFMethodPlugin
from pyscf_front.plugins.basis_sets.pyscf_basis_sets import PySCFBasisSetPlugin


class TestPluginManager:
    """PluginManager のテスト"""
    
    def test_plugin_manager_initialization(self):
        """プラグインマネージャーの初期化テスト"""
        manager = PluginManager()
        
        assert len(manager.plugins) == len(PluginType)
        for plugin_type in PluginType:
            assert plugin_type in manager.plugins
            assert isinstance(manager.plugins[plugin_type], dict)
        
        assert len(manager.initialized_plugins) == 0
    
    def test_plugin_registration(self):
        """プラグインの登録テスト"""
        manager = PluginManager()
        
        # モックプラグイン作成
        mock_plugin = Mock(spec=PluginInterface)
        mock_plugin.info = PluginInfo(
            name="TestPlugin",
            version="1.0.0",
            description="Test plugin",
            author="Test Author",
            plugin_type=PluginType.METHOD
        )
        
        # 登録
        success = manager.register_plugin(mock_plugin)
        assert success is True
        
        # 確認
        assert "TestPlugin" in manager.plugins[PluginType.METHOD]
        assert manager.plugins[PluginType.METHOD]["TestPlugin"] == mock_plugin
    
    def test_duplicate_plugin_registration(self):
        """重複プラグイン登録のテスト"""
        manager = PluginManager()
        
        # 同じ名前のプラグインを2つ作成
        mock_plugin1 = Mock(spec=PluginInterface)
        mock_plugin1.info = PluginInfo(
            name="DuplicatePlugin",
            version="1.0.0",
            description="First plugin",
            author="Author1",
            plugin_type=PluginType.METHOD
        )
        
        mock_plugin2 = Mock(spec=PluginInterface)
        mock_plugin2.info = PluginInfo(
            name="DuplicatePlugin",
            version="2.0.0", 
            description="Second plugin",
            author="Author2",
            plugin_type=PluginType.METHOD
        )
        
        # 最初の登録は成功
        assert manager.register_plugin(mock_plugin1) is True
        
        # 重複登録は失敗
        assert manager.register_plugin(mock_plugin2) is False
        
        # 最初のプラグインが残っている
        assert manager.plugins[PluginType.METHOD]["DuplicatePlugin"] == mock_plugin1
    
    def test_plugin_retrieval(self):
        """プラグインの取得テスト"""
        manager = PluginManager()
        
        mock_plugin = Mock(spec=PluginInterface)
        mock_plugin.info = PluginInfo(
            name="RetrievalTest",
            version="1.0.0",
            description="Test",
            author="Test",
            plugin_type=PluginType.BASIS_SET
        )
        
        manager.register_plugin(mock_plugin)
        
        # 取得テスト
        retrieved = manager.get_plugin(PluginType.BASIS_SET, "RetrievalTest")
        assert retrieved == mock_plugin
        
        # 存在しないプラグインの取得
        none_plugin = manager.get_plugin(PluginType.BASIS_SET, "NonExistent")
        assert none_plugin is None
    
    def test_plugin_unregistration(self):
        """プラグインの登録解除テスト"""
        manager = PluginManager()
        
        mock_plugin = Mock(spec=PluginInterface)
        mock_plugin.info = PluginInfo(
            name="UnregisterTest",
            version="1.0.0", 
            description="Test",
            author="Test",
            plugin_type=PluginType.ANALYSIS
        )
        mock_plugin.finalize.return_value = True
        
        # 登録
        manager.register_plugin(mock_plugin)
        manager.initialized_plugins.add("UnregisterTest")
        
        # 登録解除
        success = manager.unregister_plugin(PluginType.ANALYSIS, "UnregisterTest")
        assert success is True
        
        # 確認
        assert "UnregisterTest" not in manager.plugins[PluginType.ANALYSIS]
        assert "UnregisterTest" not in manager.initialized_plugins
        mock_plugin.finalize.assert_called_once()


class TestPySCFMethodPlugin:
    """PySCFMethodPlugin のテスト"""
    
    @pytest.fixture
    def mock_pyscf_plugin(self):
        """PySCFプラグインのモック"""
        with patch('pyscf_front.plugins.methods.pyscf_methods.PYSCF_AVAILABLE', True):
            plugin = PySCFMethodPlugin()
            plugin._initialized = True
            return plugin
    
    def test_plugin_info(self):
        """プラグイン情報のテスト"""
        plugin = PySCFMethodPlugin()
        info = plugin.info
        
        assert info.name == "PySCFMethods"
        assert info.plugin_type == PluginType.METHOD
        assert "pyscf" in info.dependencies[0].lower()
    
    def test_supported_methods(self, mock_pyscf_plugin):
        """サポートする手法のテスト"""
        methods = mock_pyscf_plugin.get_supported_methods()
        
        assert isinstance(methods, list)
        assert len(methods) > 0
        assert "HF" in methods
        assert "B3LYP" in methods
        assert "PBE" in methods
    
    def test_parameter_validation(self, mock_pyscf_plugin):
        """パラメータ検証のテスト"""
        # 有効なパラメータ
        valid_params = {
            'max_cycles': 50,
            'conv_threshold': 1e-8,
            'verbose': True
        }
        assert mock_pyscf_plugin.validate_parameters("HF", valid_params) is True
        
        # 無効なパラメータ
        invalid_params = {
            'max_cycles': -1,  # 負の値
            'conv_threshold': 1.0  # 大きすぎる
        }
        assert mock_pyscf_plugin.validate_parameters("HF", invalid_params) is False
    
    def test_time_estimation(self, mock_pyscf_plugin):
        """時間推定のテスト"""
        # 小さい分子
        time1 = mock_pyscf_plugin.estimate_time("HF", 3, 20)
        assert isinstance(time1, float)
        assert time1 > 0
        
        # より複雑な計算
        time2 = mock_pyscf_plugin.estimate_time("B3LYP", 10, 100)
        assert time2 > time1  # より時間がかかるはず
    
    def test_method_info_retrieval(self, mock_pyscf_plugin):
        """手法情報取得のテスト"""
        hf_info = mock_pyscf_plugin.get_method_info("HF")
        assert isinstance(hf_info, dict)
        assert hf_info['type'] == 'scf'
        assert 'description' in hf_info
        
        b3lyp_info = mock_pyscf_plugin.get_method_info("B3LYP")
        assert b3lyp_info['type'] == 'dft'
        assert b3lyp_info['xc'] == 'b3lyp'
    
    def test_feature_support(self, mock_pyscf_plugin):
        """機能サポートのテスト"""
        # 全手法で勾配計算サポート
        assert mock_pyscf_plugin.supports_feature("HF", "gradient") is True
        assert mock_pyscf_plugin.supports_feature("B3LYP", "gradient") is True
        
        # ハイブリッド汎関数の確認
        assert mock_pyscf_plugin.supports_feature("B3LYP", "hybrid") is True
        assert mock_pyscf_plugin.supports_feature("PBE", "hybrid") is False
        
        # 分散補正
        assert mock_pyscf_plugin.supports_feature("wB97X-D", "dispersion") is True
        assert mock_pyscf_plugin.supports_feature("B3LYP", "dispersion") is False


class TestPySCFBasisSetPlugin:
    """PySCFBasisSetPlugin のテスト"""
    
    @pytest.fixture
    def mock_basis_plugin(self):
        """基底関数プラグインのモック"""
        with patch('pyscf_front.plugins.basis_sets.pyscf_basis_sets.PYSCF_AVAILABLE', True):
            plugin = PySCFBasisSetPlugin()
            plugin._initialized = True
            return plugin
    
    def test_plugin_info(self):
        """プラグイン情報のテスト"""
        plugin = PySCFBasisSetPlugin()
        info = plugin.info
        
        assert info.name == "PySCFBasisSets"
        assert info.plugin_type == PluginType.BASIS_SET
        assert "pyscf" in info.dependencies[0].lower()
    
    def test_supported_basis_sets(self, mock_basis_plugin):
        """サポートする基底関数のテスト"""
        basis_sets = mock_basis_plugin.get_supported_basis_sets()
        
        assert isinstance(basis_sets, list)
        assert len(basis_sets) > 0
        assert "STO-3G" in basis_sets
        assert "6-31G" in basis_sets
        assert "6-31G(d)" in basis_sets
        assert "cc-pVDZ" in basis_sets
    
    def test_basis_validation_for_atoms(self, mock_basis_plugin):
        """原子に対する基底関数検証のテスト"""
        # 一般的な原子（H, C, N, O）でSTO-3G
        common_atoms = ['H', 'C', 'N', 'O']
        assert mock_basis_plugin.validate_basis_for_atoms("STO-3G", common_atoms) is True
        
        # 水素のみ
        hydrogen_only = ['H']
        assert mock_basis_plugin.validate_basis_for_atoms("6-31G", hydrogen_only) is True
        
        # 存在しない基底関数
        assert mock_basis_plugin.validate_basis_for_atoms("NONEXISTENT", common_atoms) is False
    
    def test_basis_info_retrieval(self, mock_basis_plugin):
        """基底関数情報取得のテスト"""
        sto3g_info = mock_basis_plugin.get_basis_info("STO-3G")
        assert isinstance(sto3g_info, dict)
        assert sto3g_info['type'] == 'minimal'
        assert sto3g_info['quality'] == 1
        
        aug_cc_pvdz_info = mock_basis_plugin.get_basis_info("aug-cc-pVDZ")
        assert aug_cc_pvdz_info['type'] == 'augmented_correlation_consistent'
        assert aug_cc_pvdz_info['quality'] > sto3g_info['quality']
    
    def test_basis_size_estimation(self, mock_basis_plugin):
        """基底関数サイズ推定のテスト"""
        atoms = ['H', 'H', 'O']  # 水分子
        
        sto3g_size = mock_basis_plugin.get_basis_size_estimate("STO-3G", atoms)
        assert isinstance(sto3g_size, int)
        assert sto3g_size > 0
        
        larger_basis_size = mock_basis_plugin.get_basis_size_estimate("6-31G(d)", atoms)
        assert larger_basis_size > sto3g_size
    
    def test_basis_recommendations(self, mock_basis_plugin):
        """基底関数推奨のテスト"""
        # 小さい分子でのHF計算
        hf_small = mock_basis_plugin.recommend_basis_for_method("HF", 5)
        assert isinstance(hf_small, list)
        assert len(hf_small) > 0
        
        # 大きい分子でのDFT計算
        dft_large = mock_basis_plugin.recommend_basis_for_method("B3LYP", 50)
        assert isinstance(dft_large, list)
        
        # より計算負荷の低い基底関数が推奨されるはず
        assert any("STO-3G" in rec or "6-31G" in rec for rec in dft_large)
    
    def test_basis_quality_comparison(self, mock_basis_plugin):
        """基底関数品質比較のテスト"""
        quality_sto3g = mock_basis_plugin.get_basis_quality("STO-3G")
        quality_631gd = mock_basis_plugin.get_basis_quality("6-31G(d)")
        quality_cc_pvtz = mock_basis_plugin.get_basis_quality("cc-pVTZ")
        
        assert quality_sto3g < quality_631gd < quality_cc_pvtz
    
    def test_basis_family(self, mock_basis_plugin):
        """基底関数ファミリーのテスト"""
        assert mock_basis_plugin.get_basis_family("6-31G") == "Pople"
        assert mock_basis_plugin.get_basis_family("cc-pVDZ") == "correlation-consistent"
        assert mock_basis_plugin.get_basis_family("def2-SVP") == "def2"


class TestPluginSystemIntegration:
    """プラグインシステム統合テスト"""
    
    def test_global_plugin_manager(self):
        """グローバルプラグインマネージャーのテスト"""
        manager1 = get_plugin_manager()
        manager2 = get_plugin_manager()
        
        # 同じインスタンスが返される
        assert manager1 is manager2
    
    def test_plugin_registration_via_global_function(self):
        """グローバル関数経由でのプラグイン登録テスト"""
        mock_plugin = Mock(spec=PluginInterface)
        mock_plugin.info = PluginInfo(
            name="GlobalTest",
            version="1.0.0",
            description="Test",
            author="Test",
            plugin_type=PluginType.ENGINE
        )
        
        # グローバル関数で登録
        success = register_plugin(mock_plugin)
        assert success is True
        
        # グローバルマネージャーで確認
        manager = get_plugin_manager()
        retrieved = manager.get_plugin(PluginType.ENGINE, "GlobalTest")
        assert retrieved == mock_plugin
    
    def test_method_and_basis_plugin_integration(self):
        """手法と基底関数プラグインの統合テスト"""
        with patch('pyscf_front.plugins.methods.pyscf_methods.PYSCF_AVAILABLE', True), \
             patch('pyscf_front.plugins.basis_sets.pyscf_basis_sets.PYSCF_AVAILABLE', True):
            
            manager = PluginManager()
            
            # プラグインを登録
            method_plugin = PySCFMethodPlugin()
            basis_plugin = PySCFBasisSetPlugin()
            
            manager.register_plugin(method_plugin)
            manager.register_plugin(basis_plugin)
            
            manager.initialize_plugin(PluginType.METHOD, "PySCFMethods")
            manager.initialize_plugin(PluginType.BASIS_SET, "PySCFBasisSets")
            
            # 統合機能のテスト
            methods = manager.get_available_methods()
            basis_sets = manager.get_available_basis_sets()
            
            assert len(methods) > 0
            assert len(basis_sets) > 0
            assert "HF" in methods
            assert "STO-3G" in basis_sets
            
            # プラグイン検索のテスト
            found_method_plugin = manager.find_method_plugin("B3LYP")
            found_basis_plugin = manager.find_basis_plugin("6-31G")
            
            assert found_method_plugin is not None
            assert found_basis_plugin is not None
            assert isinstance(found_method_plugin, PySCFMethodPlugin)
            assert isinstance(found_basis_plugin, PySCFBasisSetPlugin)
    
    def test_plugin_unavailable_handling(self):
        """プラグイン利用不可時のハンドリングテスト"""
        manager = PluginManager()
        
        # 利用不可能なプラグイン
        mock_plugin = Mock(spec=PluginInterface)
        mock_plugin.info = PluginInfo(
            name="UnavailablePlugin",
            version="1.0.0",
            description="Test",
            author="Test",
            plugin_type=PluginType.METHOD
        )
        mock_plugin.is_available.return_value = False
        
        manager.register_plugin(mock_plugin)
        
        # 利用可能な手法リストに含まれない
        methods = manager.get_available_methods()
        assert len(methods) == 0  # 利用可能なプラグインがないため
        
        # プラグイン検索でも見つからない
        found = manager.find_method_plugin("SomeMethod")
        assert found is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])