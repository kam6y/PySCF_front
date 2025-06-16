"""
プラグインシステムの基底クラス
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass
from enum import Enum
from loguru import logger


class PluginType(Enum):
    """プラグインタイプ"""
    METHOD = "method"
    BASIS_SET = "basis_set"
    ANALYSIS = "analysis"
    ENGINE = "engine"


@dataclass
class PluginInfo:
    """プラグイン情報"""
    name: str
    version: str
    description: str
    author: str
    plugin_type: PluginType
    dependencies: List[str] = None
    
    def __post_init__(self):
        if self.dependencies is None:
            self.dependencies = []


class PluginInterface(ABC):
    """プラグインの基底インターフェース"""
    
    @property
    @abstractmethod
    def info(self) -> PluginInfo:
        """プラグイン情報を返す"""
        pass
    
    @abstractmethod
    def initialize(self, config: Dict[str, Any] = None) -> bool:
        """プラグインを初期化する"""
        pass
    
    @abstractmethod
    def finalize(self) -> bool:
        """プラグインを終了処理する"""
        pass
    
    @abstractmethod
    def is_available(self) -> bool:
        """プラグインが利用可能かチェック"""
        pass


class CalculationMethodPlugin(PluginInterface):
    """計算手法プラグインの基底クラス"""
    
    @abstractmethod
    def get_supported_methods(self) -> List[str]:
        """サポートする計算手法のリストを返す"""
        pass
    
    @abstractmethod
    def create_calculator(self, method: str, **kwargs) -> Any:
        """指定された手法の計算オブジェクトを作成"""
        pass
    
    @abstractmethod
    def validate_parameters(self, method: str, parameters: Dict[str, Any]) -> bool:
        """パラメータの妥当性をチェック"""
        pass
    
    @abstractmethod
    def estimate_time(self, method: str, molecule_size: int, basis_size: int) -> float:
        """計算時間を推定（秒単位）"""
        pass


class BasisSetPlugin(PluginInterface):
    """基底関数プラグインの基底クラス"""
    
    @abstractmethod
    def get_supported_basis_sets(self) -> List[str]:
        """サポートする基底関数のリストを返す"""
        pass
    
    @abstractmethod
    def get_basis_info(self, basis_name: str) -> Dict[str, Any]:
        """基底関数の詳細情報を返す"""
        pass
    
    @abstractmethod
    def validate_basis_for_atoms(self, basis_name: str, atom_symbols: List[str]) -> bool:
        """指定された原子に対して基底関数が利用可能かチェック"""
        pass


class AnalysisPlugin(PluginInterface):
    """解析プラグインの基底クラス"""
    
    @abstractmethod
    def get_supported_analyses(self) -> List[str]:
        """サポートする解析手法のリストを返す"""
        pass
    
    @abstractmethod
    def analyze(self, analysis_type: str, calculation_results: Dict[str, Any], **kwargs) -> Dict[str, Any]:
        """解析を実行"""
        pass
    
    @abstractmethod
    def format_results(self, analysis_results: Dict[str, Any], format_type: str = "text") -> str:
        """解析結果をフォーマットして返す"""
        pass


class CalculationEnginePlugin(PluginInterface):
    """計算エンジンプラグインの基底クラス"""
    
    @abstractmethod
    def get_engine_capabilities(self) -> Dict[str, Any]:
        """エンジンの機能を返す"""
        pass
    
    @abstractmethod
    def create_worker(self, job_config: Dict[str, Any]) -> Any:
        """計算ワーカーを作成"""
        pass
    
    @abstractmethod
    def supports_method_basis_combination(self, method: str, basis_set: str) -> bool:
        """手法と基底関数の組み合わせがサポートされているかチェック"""
        pass


class PluginManager:
    """プラグイン管理クラス"""
    
    def __init__(self):
        self.plugins: Dict[PluginType, Dict[str, PluginInterface]] = {
            plugin_type: {} for plugin_type in PluginType
        }
        self.initialized_plugins: set = set()
    
    def register_plugin(self, plugin: PluginInterface) -> bool:
        """プラグインを登録"""
        try:
            plugin_type = plugin.info.plugin_type
            plugin_name = plugin.info.name
            
            if plugin_name in self.plugins[plugin_type]:
                logger.debug(f"Plugin '{plugin_name}' of type '{plugin_type.value}' already registered, skipping")
                return True  # 既に登録済みの場合は成功として扱う
            
            self.plugins[plugin_type][plugin_name] = plugin
            return True
            
        except Exception as e:
            logger.error(f"Failed to register plugin: {e}")
            return False
    
    def unregister_plugin(self, plugin_type: PluginType, plugin_name: str) -> bool:
        """プラグインの登録を解除"""
        try:
            if plugin_name in self.plugins[plugin_type]:
                plugin = self.plugins[plugin_type][plugin_name]
                if plugin_name in self.initialized_plugins:
                    plugin.finalize()
                    self.initialized_plugins.remove(plugin_name)
                del self.plugins[plugin_type][plugin_name]
                return True
            return False
            
        except Exception as e:
            print(f"Failed to unregister plugin: {e}")
            return False
    
    def get_plugin(self, plugin_type: PluginType, plugin_name: str) -> Optional[PluginInterface]:
        """プラグインを取得"""
        return self.plugins[plugin_type].get(plugin_name)
    
    def get_plugins_by_type(self, plugin_type: PluginType) -> Dict[str, PluginInterface]:
        """指定されたタイプのプラグインをすべて取得"""
        return self.plugins[plugin_type].copy()
    
    def initialize_plugin(self, plugin_type: PluginType, plugin_name: str, config: Dict[str, Any] = None) -> bool:
        """プラグインを初期化"""
        plugin = self.get_plugin(plugin_type, plugin_name)
        if plugin is None:
            return False
        
        try:
            if plugin.initialize(config or {}):
                self.initialized_plugins.add(plugin_name)
                return True
            return False
            
        except Exception as e:
            print(f"Failed to initialize plugin '{plugin_name}': {e}")
            return False
    
    def is_plugin_available(self, plugin_type: PluginType, plugin_name: str) -> bool:
        """プラグインが利用可能かチェック"""
        plugin = self.get_plugin(plugin_type, plugin_name)
        if plugin is None:
            return False
        
        return plugin.is_available()
    
    def get_available_methods(self) -> List[str]:
        """利用可能な計算手法の一覧を取得"""
        methods = []
        for plugin in self.plugins[PluginType.METHOD].values():
            if plugin.is_available():
                methods.extend(plugin.get_supported_methods())
        return sorted(list(set(methods)))
    
    def get_available_basis_sets(self) -> List[str]:
        """利用可能な基底関数の一覧を取得"""
        basis_sets = []
        for plugin in self.plugins[PluginType.BASIS_SET].values():
            if plugin.is_available():
                basis_sets.extend(plugin.get_supported_basis_sets())
        return sorted(list(set(basis_sets)))
    
    def get_available_analyses(self) -> List[str]:
        """利用可能な解析手法の一覧を取得"""
        analyses = []
        for plugin in self.plugins[PluginType.ANALYSIS].values():
            if plugin.is_available():
                analyses.extend(plugin.get_supported_analyses())
        return sorted(list(set(analyses)))
    
    def find_method_plugin(self, method: str) -> Optional[CalculationMethodPlugin]:
        """指定された手法をサポートするプラグインを見つける"""
        for plugin in self.plugins[PluginType.METHOD].values():
            if isinstance(plugin, CalculationMethodPlugin) and plugin.is_available():
                if method in plugin.get_supported_methods():
                    return plugin
        return None
    
    def find_basis_plugin(self, basis_set: str) -> Optional[BasisSetPlugin]:
        """指定された基底関数をサポートするプラグインを見つける"""
        for plugin in self.plugins[PluginType.BASIS_SET].values():
            if isinstance(plugin, BasisSetPlugin) and plugin.is_available():
                if basis_set in plugin.get_supported_basis_sets():
                    return plugin
        return None
    
    def find_analysis_plugin(self, analysis_type: str) -> Optional[AnalysisPlugin]:
        """指定された解析手法をサポートするプラグインを見つける"""
        for plugin in self.plugins[PluginType.ANALYSIS].values():
            if isinstance(plugin, AnalysisPlugin) and plugin.is_available():
                if analysis_type in plugin.get_supported_analyses():
                    return plugin
        return None
    
    def finalize_all(self):
        """すべてのプラグインを終了処理"""
        for plugin_type in PluginType:
            for plugin_name, plugin in self.plugins[plugin_type].items():
                if plugin_name in self.initialized_plugins:
                    try:
                        plugin.finalize()
                    except Exception as e:
                        logger.error(f"Error finalizing plugin '{plugin_name}': {e}")
        
        self.initialized_plugins.clear()


# グローバルプラグインマネージャー
_global_plugin_manager = None


def get_plugin_manager() -> PluginManager:
    """グローバルプラグインマネージャーを取得"""
    global _global_plugin_manager
    if _global_plugin_manager is None:
        _global_plugin_manager = PluginManager()
    return _global_plugin_manager


def register_plugin(plugin: PluginInterface) -> bool:
    """プラグインをグローバルマネージャーに登録"""
    return get_plugin_manager().register_plugin(plugin)