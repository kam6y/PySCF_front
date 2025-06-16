"""
PySCF計算手法プラグイン
"""
import math
from typing import Dict, Any, List, Optional
from loguru import logger

from ..base import CalculationMethodPlugin, PluginInfo, PluginType

try:
    import pyscf
    from pyscf import gto, scf, dft
    PYSCF_AVAILABLE = True
    logger.info(f"PySCF {pyscf.__version__} loaded for plugin")
except ImportError:
    logger.warning("PySCF not available for plugin")
    PYSCF_AVAILABLE = False
    pyscf = None
    gto = scf = dft = None


class PySCFMethodPlugin(CalculationMethodPlugin):
    """PySCF計算手法プラグイン"""
    
    def __init__(self):
        self._initialized = False
        self._config = {}
        
        # サポートする手法の定義
        self._supported_methods = {
            'HF': {
                'type': 'scf',
                'class': 'RHF',  # または UHF
                'description': 'Hartree-Fock method',
                'complexity': 1.0
            },
            'B3LYP': {
                'type': 'dft',
                'class': 'RKS',  # または UKS
                'xc': 'b3lyp',
                'description': 'B3LYP density functional theory',
                'complexity': 2.0
            },
            'PBE': {
                'type': 'dft',
                'class': 'RKS',
                'xc': 'pbe',
                'description': 'PBE density functional theory',
                'complexity': 1.8
            },
            'PBE0': {
                'type': 'dft',
                'class': 'RKS',
                'xc': 'pbe0',
                'description': 'PBE0 hybrid density functional theory',
                'complexity': 2.2
            },
            'M06': {
                'type': 'dft',
                'class': 'RKS',
                'xc': 'm06',
                'description': 'M06 meta-GGA functional',
                'complexity': 2.5
            },
            'wB97X-D': {
                'type': 'dft',
                'class': 'RKS',
                'xc': 'wb97x-d',
                'description': 'wB97X-D range-separated functional with dispersion',
                'complexity': 3.0
            }
        }
    
    @property
    def info(self) -> PluginInfo:
        """プラグイン情報"""
        return PluginInfo(
            name="PySCFMethods",
            version="1.0.0",
            description="PySCF-based quantum chemistry calculation methods",
            author="PySCF_Front Team",
            plugin_type=PluginType.METHOD,
            dependencies=["pyscf>=2.0.0"]
        )
    
    def initialize(self, config: Dict[str, Any] = None) -> bool:
        """プラグイン初期化"""
        try:
            if not PYSCF_AVAILABLE:
                logger.error("PySCF is not available for initialization")
                return False
            
            self._config = config or {}
            
            # PySCFのバージョンチェック
            pyscf_version = tuple(map(int, pyscf.__version__.split('.')[:2]))
            if pyscf_version < (2, 0):
                logger.warning(f"PySCF version {pyscf.__version__} may not be fully supported")
            
            # 設定の適用
            if 'verbose' in self._config:
                # PySCFの詳細度設定
                pass
            
            self._initialized = True
            logger.info("PySCF methods plugin initialized successfully")
            return True
            
        except Exception as e:
            logger.error(f"Failed to initialize PySCF methods plugin: {e}")
            return False
    
    def finalize(self) -> bool:
        """プラグイン終了処理"""
        try:
            self._initialized = False
            self._config = {}
            logger.info("PySCF methods plugin finalized")
            return True
        except Exception as e:
            logger.error(f"Failed to finalize PySCF methods plugin: {e}")
            return False
    
    def is_available(self) -> bool:
        """プラグインが利用可能かチェック"""
        return PYSCF_AVAILABLE and self._initialized
    
    def get_supported_methods(self) -> List[str]:
        """サポートする計算手法のリスト"""
        if not self.is_available():
            return []
        
        return list(self._supported_methods.keys())
    
    def create_calculator(self, method: str, **kwargs) -> Any:
        """計算オブジェクトを作成"""
        if not self.is_available():
            raise RuntimeError("PySCF methods plugin is not available")
        
        if method not in self._supported_methods:
            raise ValueError(f"Unsupported method: {method}")
        
        try:
            method_info = self._supported_methods[method]
            mol = kwargs.get('mol')
            
            if mol is None:
                raise ValueError("Molecule object (mol) is required")
            
            # 計算オブジェクトの作成
            if method_info['type'] == 'scf':
                if mol.spin == 0:
                    calc = scf.RHF(mol)
                else:
                    calc = scf.UHF(mol)
            
            elif method_info['type'] == 'dft':
                if mol.spin == 0:
                    calc = dft.RKS(mol)
                else:
                    calc = dft.UKS(mol)
                
                # 交換相関汎関数の設定
                calc.xc = method_info['xc']
            
            else:
                raise ValueError(f"Unknown method type: {method_info['type']}")
            
            # 追加パラメータの設定
            if 'max_cycles' in kwargs:
                calc.max_cycle = kwargs['max_cycles']
            
            if 'conv_threshold' in kwargs:
                calc.conv_tol = kwargs['conv_threshold']
            
            if 'verbose' in kwargs:
                calc.verbose = kwargs['verbose']
            
            return calc
            
        except Exception as e:
            logger.error(f"Failed to create calculator for method {method}: {e}")
            raise
    
    def validate_parameters(self, method: str, parameters: Dict[str, Any]) -> bool:
        """パラメータの妥当性チェック"""
        if method not in self._supported_methods:
            return False
        
        try:
            # 基本的なパラメータチェック
            if 'max_cycles' in parameters:
                max_cycles = parameters['max_cycles']
                if not isinstance(max_cycles, int) or max_cycles < 1 or max_cycles > 1000:
                    return False
            
            if 'conv_threshold' in parameters:
                conv_tol = parameters['conv_threshold']
                if not isinstance(conv_tol, (int, float)) or conv_tol <= 0 or conv_tol > 1e-3:
                    return False
            
            if 'verbose' in parameters:
                verbose = parameters['verbose']
                if not isinstance(verbose, (int, bool)):
                    return False
            
            return True
            
        except Exception as e:
            logger.error(f"Parameter validation failed for method {method}: {e}")
            return False
    
    def estimate_time(self, method: str, molecule_size: int, basis_size: int) -> float:
        """計算時間推定（秒単位）"""
        if method not in self._supported_methods:
            return 0.0
        
        method_info = self._supported_methods[method]
        complexity = method_info['complexity']
        
        # 基本時間の推定（経験式）
        # basis_size^3 * molecule_size * complexity_factor
        base_time = (basis_size ** 2.5) * molecule_size * complexity / 1000.0
        
        # 最小時間と最大時間の制限
        estimated_time = max(1.0, min(base_time, 3600.0 * 24))  # 1秒から24時間
        
        return estimated_time
    
    def get_method_info(self, method: str) -> Dict[str, Any]:
        """手法の詳細情報を取得"""
        if method not in self._supported_methods:
            return {}
        
        return self._supported_methods[method].copy()
    
    def get_default_parameters(self, method: str) -> Dict[str, Any]:
        """手法のデフォルトパラメータを取得"""
        if method not in self._supported_methods:
            return {}
        
        method_info = self._supported_methods[method]
        
        defaults = {
            'max_cycles': 50,
            'conv_threshold': 1e-8,
            'verbose': 4
        }
        
        # 手法固有のデフォルト設定
        if method_info['type'] == 'dft':
            defaults.update({
                'grids_level': 3,
                'max_memory': 4000  # MB
            })
        
        return defaults
    
    def supports_feature(self, method: str, feature: str) -> bool:
        """手法が特定の機能をサポートするかチェック"""
        if method not in self._supported_methods:
            return False
        
        method_info = self._supported_methods[method]
        
        feature_support = {
            'gradient': True,  # 全手法で勾配計算可能
            'hessian': method_info['type'] in ['scf', 'dft'],
            'dispersion': method in ['wB97X-D'],
            'meta_gga': method in ['M06'],
            'range_separated': method in ['wB97X-D'],
            'hybrid': method in ['B3LYP', 'PBE0', 'M06', 'wB97X-D']
        }
        
        return feature_support.get(feature, False)
    
    def get_citation(self, method: str) -> str:
        """手法の引用情報を取得"""
        citations = {
            'HF': "Hartree, D. R. (1928). The Wave Mechanics of an Atom with a Non-Coulomb Central Field.",
            'B3LYP': "Becke, A. D. (1993). Density‐functional thermochemistry. III.",
            'PBE': "Perdew, J. P., Burke, K., & Ernzerhof, M. (1996). Generalized gradient approximation.",
            'PBE0': "Adamo, C., & Barone, V. (1999). Toward reliable density functional methods.",
            'M06': "Zhao, Y., & Truhlar, D. G. (2008). The M06 suite of density functionals.",
            'wB97X-D': "Chai, J. D., & Head-Gordon, M. (2008). Long-range corrected hybrid density functionals."
        }
        
        return citations.get(method, "No citation available")