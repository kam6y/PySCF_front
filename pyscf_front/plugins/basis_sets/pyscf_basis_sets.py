"""
PySCF基底関数プラグイン
"""
from typing import Dict, Any, List
from loguru import logger

from ..base import BasisSetPlugin, PluginInfo, PluginType

try:
    import pyscf
    from pyscf import gto
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False
    pyscf = gto = None


class PySCFBasisSetPlugin(BasisSetPlugin):
    """PySCF基底関数プラグイン"""
    
    def __init__(self):
        self._initialized = False
        self._config = {}
        
        # サポートする基底関数の定義
        self._supported_basis_sets = {
            # 最小基底関数
            'STO-3G': {
                'type': 'minimal',
                'description': 'STO-3G minimal basis set',
                'quality': 1,
                'atoms': 'H-Ne',  # 水素からネオンまで
                'size_factor': 1.0
            },
            
            # 分割価基底関数
            '3-21G': {
                'type': 'split_valence',
                'description': '3-21G split-valence basis set',
                'quality': 2,
                'atoms': 'H-Ar',
                'size_factor': 2.5
            },
            '6-31G': {
                'type': 'split_valence',
                'description': '6-31G split-valence basis set',
                'quality': 3,
                'atoms': 'H-Kr',
                'size_factor': 3.0
            },
            '6-31G(d)': {
                'type': 'split_valence_polarized',
                'description': '6-31G with d-polarization functions',
                'quality': 4,
                'atoms': 'H-Kr',
                'size_factor': 4.0
            },
            '6-31G(d,p)': {
                'type': 'split_valence_polarized',
                'description': '6-31G with d,p-polarization functions',
                'quality': 4,
                'atoms': 'H-Kr',
                'size_factor': 4.5
            },
            '6-31+G(d)': {
                'type': 'augmented_split_valence',
                'description': '6-31G(d) with diffuse functions',
                'quality': 5,
                'atoms': 'H-Kr',
                'size_factor': 5.0
            },
            '6-31+G(d,p)': {
                'type': 'augmented_split_valence',
                'description': '6-31G(d,p) with diffuse functions',
                'quality': 5,
                'atoms': 'H-Kr',
                'size_factor': 5.5
            },
            '6-31++G(d,p)': {
                'type': 'augmented_split_valence',
                'description': '6-31G(d,p) with diffuse functions on all atoms',
                'quality': 6,
                'atoms': 'H-Kr',
                'size_factor': 6.0
            },
            
            # Pople基底関数（大きめ）
            '6-311G': {
                'type': 'triple_zeta',
                'description': '6-311G triple-zeta basis set',
                'quality': 5,
                'atoms': 'H-Kr',
                'size_factor': 6.0
            },
            '6-311G(d,p)': {
                'type': 'triple_zeta_polarized',
                'description': '6-311G with d,p-polarization functions',
                'quality': 6,
                'atoms': 'H-Kr',
                'size_factor': 7.0
            },
            '6-311+G(d,p)': {
                'type': 'augmented_triple_zeta',
                'description': '6-311G(d,p) with diffuse functions',
                'quality': 7,
                'atoms': 'H-Kr',
                'size_factor': 8.0
            },
            '6-311++G(d,p)': {
                'type': 'augmented_triple_zeta',
                'description': '6-311G(d,p) with diffuse functions on all atoms',
                'quality': 8,
                'atoms': 'H-Kr',
                'size_factor': 8.5
            },
            
            # correlation-consistent基底関数
            'cc-pVDZ': {
                'type': 'correlation_consistent',
                'description': 'Correlation-consistent double-zeta',
                'quality': 5,
                'atoms': 'H-Ar',
                'size_factor': 5.0
            },
            'cc-pVTZ': {
                'type': 'correlation_consistent',
                'description': 'Correlation-consistent triple-zeta',
                'quality': 7,
                'atoms': 'H-Ar',
                'size_factor': 10.0
            },
            'cc-pVQZ': {
                'type': 'correlation_consistent',
                'description': 'Correlation-consistent quadruple-zeta',
                'quality': 9,
                'atoms': 'H-Ar',
                'size_factor': 20.0
            },
            'aug-cc-pVDZ': {
                'type': 'augmented_correlation_consistent',
                'description': 'Augmented correlation-consistent double-zeta',
                'quality': 6,
                'atoms': 'H-Ar',
                'size_factor': 7.0
            },
            'aug-cc-pVTZ': {
                'type': 'augmented_correlation_consistent',
                'description': 'Augmented correlation-consistent triple-zeta',
                'quality': 8,
                'atoms': 'H-Ar',
                'size_factor': 15.0
            },
            
            # 密度汎関数用基底関数
            'def2-SVP': {
                'type': 'def2',
                'description': 'def2 split-valence polarized',
                'quality': 4,
                'atoms': 'H-Kr',
                'size_factor': 4.0
            },
            'def2-SVPD': {
                'type': 'def2',
                'description': 'def2-SVP with diffuse functions',
                'quality': 5,
                'atoms': 'H-Kr',
                'size_factor': 5.0
            },
            'def2-TZVP': {
                'type': 'def2',
                'description': 'def2 triple-zeta valence polarized',
                'quality': 6,
                'atoms': 'H-Kr',
                'size_factor': 7.0
            },
            'def2-TZVPD': {
                'type': 'def2',
                'description': 'def2-TZVP with diffuse functions',
                'quality': 7,
                'atoms': 'H-Kr',
                'size_factor': 8.0
            }
        }
        
        # 原子記号のマッピング
        self._atom_ranges = {
            'H-He': ['H', 'He'],
            'H-Ne': ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'],
            'H-Ar': ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                     'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
            'H-Kr': ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                     'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                     'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
                     'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']
        }
    
    @property
    def info(self) -> PluginInfo:
        """プラグイン情報"""
        return PluginInfo(
            name="PySCFBasisSets",
            version="1.0.0",
            description="PySCF-based basis set library",
            author="PySCF_Front Team",
            plugin_type=PluginType.BASIS_SET,
            dependencies=["pyscf>=2.0.0"]
        )
    
    def initialize(self, config: Dict[str, Any] = None) -> bool:
        """プラグイン初期化"""
        try:
            if not PYSCF_AVAILABLE:
                logger.error("PySCF is not available for basis set plugin")
                return False
            
            self._config = config or {}
            self._initialized = True
            logger.info("PySCF basis sets plugin initialized successfully")
            return True
            
        except Exception as e:
            logger.error(f"Failed to initialize PySCF basis sets plugin: {e}")
            return False
    
    def finalize(self) -> bool:
        """プラグイン終了処理"""
        try:
            self._initialized = False
            self._config = {}
            logger.info("PySCF basis sets plugin finalized")
            return True
        except Exception as e:
            logger.error(f"Failed to finalize PySCF basis sets plugin: {e}")
            return False
    
    def is_available(self) -> bool:
        """プラグインが利用可能かチェック"""
        return PYSCF_AVAILABLE and self._initialized
    
    def get_supported_basis_sets(self) -> List[str]:
        """サポートする基底関数のリスト"""
        if not self.is_available():
            return []
        
        return list(self._supported_basis_sets.keys())
    
    def get_basis_info(self, basis_name: str) -> Dict[str, Any]:
        """基底関数の詳細情報"""
        if basis_name not in self._supported_basis_sets:
            return {}
        
        return self._supported_basis_sets[basis_name].copy()
    
    def validate_basis_for_atoms(self, basis_name: str, atom_symbols: List[str]) -> bool:
        """指定された原子に対して基底関数が利用可能かチェック"""
        if not self.is_available():
            return False
        
        if basis_name not in self._supported_basis_sets:
            return False
        
        basis_info = self._supported_basis_sets[basis_name]
        atom_range = basis_info['atoms']
        
        if atom_range not in self._atom_ranges:
            return False
        
        supported_atoms = self._atom_ranges[atom_range]
        
        # すべての原子がサポートされているかチェック
        for atom in atom_symbols:
            if atom not in supported_atoms:
                return False
        
        return True
    
    def get_basis_size_estimate(self, basis_name: str, atom_symbols: List[str]) -> int:
        """基底関数のサイズ推定"""
        if not self.validate_basis_for_atoms(basis_name, atom_symbols):
            return 0
        
        basis_info = self._supported_basis_sets[basis_name]
        size_factor = basis_info['size_factor']
        
        # 簡易的なサイズ推定
        # 実際の基底関数サイズは原子の種類によって大きく異なる
        size_per_atom = {
            'H': 1, 'He': 2,
            'Li': 5, 'Be': 5, 'B': 5, 'C': 5, 'N': 5, 'O': 5, 'F': 5, 'Ne': 5,
            'Na': 9, 'Mg': 9, 'Al': 9, 'Si': 9, 'P': 9, 'S': 9, 'Cl': 9, 'Ar': 9
        }
        
        total_size = 0
        for atom in atom_symbols:
            base_size = size_per_atom.get(atom, 5)  # デフォルトは5
            total_size += int(base_size * size_factor)
        
        return total_size
    
    def get_basis_quality(self, basis_name: str) -> int:
        """基底関数の品質レベル（1-10）"""
        if basis_name not in self._supported_basis_sets:
            return 0
        
        return self._supported_basis_sets[basis_name]['quality']
    
    def recommend_basis_for_method(self, method: str, molecule_size: int) -> List[str]:
        """手法と分子サイズに応じた推奨基底関数"""
        recommendations = []
        
        # 手法に応じた推奨
        if method in ['HF']:
            # Hartree-Fock用推奨
            if molecule_size <= 10:
                recommendations = ['6-31G(d)', '6-311G(d,p)', 'cc-pVDZ']
            elif molecule_size <= 30:
                recommendations = ['6-31G', '6-31G(d)', 'def2-SVP']
            else:
                recommendations = ['STO-3G', '3-21G', '6-31G']
        
        elif method in ['B3LYP', 'PBE', 'PBE0']:
            # DFT用推奨
            if molecule_size <= 10:
                recommendations = ['def2-TZVP', '6-311+G(d,p)', 'cc-pVTZ']
            elif molecule_size <= 30:
                recommendations = ['def2-SVP', '6-31G(d)', '6-31+G(d)']
            else:
                recommendations = ['6-31G', 'def2-SVP', '3-21G']
        
        elif method in ['M06', 'wB97X-D']:
            # メタGGAやrange-separated汎関数用
            if molecule_size <= 10:
                recommendations = ['def2-TZVPD', 'aug-cc-pVTZ', '6-311++G(d,p)']
            elif molecule_size <= 30:
                recommendations = ['def2-SVPD', '6-31+G(d,p)', 'aug-cc-pVDZ']
            else:
                recommendations = ['6-31+G(d)', 'def2-SVP', '6-31G(d)']
        
        else:
            # デフォルト推奨
            recommendations = ['6-31G(d)', 'def2-SVP', 'cc-pVDZ']
        
        # 利用可能な基底関数のみを返す
        available_basis = self.get_supported_basis_sets()
        return [basis for basis in recommendations if basis in available_basis]
    
    def get_basis_family(self, basis_name: str) -> str:
        """基底関数ファミリーを取得"""
        if basis_name not in self._supported_basis_sets:
            return "unknown"
        
        basis_type = self._supported_basis_sets[basis_name]['type']
        
        if 'correlation_consistent' in basis_type:
            return "correlation-consistent"
        elif 'def2' in basis_type:
            return "def2"
        elif basis_name.startswith('6-3'):
            return "Pople"
        else:
            return "other"
    
    def compare_basis_sets(self, basis1: str, basis2: str) -> Dict[str, Any]:
        """2つの基底関数を比較"""
        if basis1 not in self._supported_basis_sets or basis2 not in self._supported_basis_sets:
            return {}
        
        info1 = self._supported_basis_sets[basis1]
        info2 = self._supported_basis_sets[basis2]
        
        return {
            'quality_comparison': {
                basis1: info1['quality'],
                basis2: info2['quality'],
                'better': basis1 if info1['quality'] > info2['quality'] else basis2
            },
            'size_comparison': {
                basis1: info1['size_factor'],
                basis2: info2['size_factor'],
                'smaller': basis1 if info1['size_factor'] < info2['size_factor'] else basis2
            },
            'type_comparison': {
                basis1: info1['type'],
                basis2: info2['type']
            },
            'atom_coverage': {
                basis1: info1['atoms'],
                basis2: info2['atoms']
            }
        }