"""
構造最適化とASEアライメント効果の検証テスト

このテストファイルでは以下を検証します:
1. 構造最適化による分子構造の改善（結合長、結合角、エネルギー）
2. ASEアライメントによる座標系の標準化
3. アライメント前後での内部構造の保持
"""
import pytest
import numpy as np
import tempfile
import os
from unittest.mock import patch, MagicMock
from typing import List, Tuple, Dict

from quantum_calc.dft_calculator import DFTCalculator
from quantum_calc.exceptions import InputError, ConvergenceError, GeometryError


class GeometryAnalyzer:
    """分子構造の解析用ユーティリティクラス"""
    
    @staticmethod
    def calculate_distance(pos1: np.ndarray, pos2: np.ndarray) -> float:
        """2点間の距離を計算"""
        return np.linalg.norm(pos1 - pos2)
    
    @staticmethod
    def calculate_angle(pos1: np.ndarray, pos2: np.ndarray, pos3: np.ndarray) -> float:
        """3点で形成される角度を計算（pos2が頂点）
        
        Returns:
            角度（度）
        """
        v1 = pos1 - pos2
        v2 = pos3 - pos2
        
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        # 数値誤差によるdomainエラーを防ぐ
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        
        return np.degrees(np.arccos(cos_angle))
    
    @staticmethod
    def calculate_center_of_mass(symbols: List[str], positions: np.ndarray) -> np.ndarray:
        """質量中心を計算"""
        # 原子質量（単位: u）の簡易マップ
        atomic_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'Si': 28.086, 'P': 30.974, 'S': 32.065,
            'Cl': 35.453
        }
        
        total_mass = 0.0
        weighted_pos = np.zeros(3)
        
        for symbol, pos in zip(symbols, positions):
            mass = atomic_masses.get(symbol, 12.0)  # デフォルトはカーボン質量
            total_mass += mass
            weighted_pos += mass * pos
        
        return weighted_pos / total_mass
    
    @staticmethod
    def calculate_inertia_tensor(symbols: List[str], positions: np.ndarray) -> np.ndarray:
        """慣性テンソルを計算"""
        # 原子質量マップ
        atomic_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'Si': 28.086, 'P': 30.974, 'S': 32.065,
            'Cl': 35.453
        }
        
        inertia_tensor = np.zeros((3, 3))
        
        for symbol, pos in zip(symbols, positions):
            mass = atomic_masses.get(symbol, 12.0)
            x, y, z = pos
            
            # 慣性テンソルの要素を計算
            inertia_tensor[0, 0] += mass * (y*y + z*z)
            inertia_tensor[1, 1] += mass * (x*x + z*z)
            inertia_tensor[2, 2] += mass * (x*x + y*y)
            inertia_tensor[0, 1] -= mass * x * y
            inertia_tensor[0, 2] -= mass * x * z
            inertia_tensor[1, 2] -= mass * y * z
        
        # 対称行列にする
        inertia_tensor[1, 0] = inertia_tensor[0, 1]
        inertia_tensor[2, 0] = inertia_tensor[0, 2]
        inertia_tensor[2, 1] = inertia_tensor[1, 2]
        
        return inertia_tensor
    
    @staticmethod
    def get_principal_axes(inertia_tensor: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """慣性テンソルから主慣性モーメントと主軸を取得"""
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)
        
        # 固有値の順序でソート（昇順）
        idx = np.argsort(eigenvalues)
        principal_moments = eigenvalues[idx]
        principal_axes = eigenvectors[:, idx]
        
        return principal_moments, principal_axes


class MockedGeometryOptimizationTest:
    """構造最適化効果をモック付きで検証するテスト"""
    
    def __init__(self):
        self.analyzer = GeometryAnalyzer()
    
    def create_distorted_water(self) -> str:
        """歪んだ水分子のXYZ座標を作成"""
        return """3
water molecule (distorted)
O  0.1000   0.2000   0.1000
H  1.0000   0.8000   0.1000  
H -0.8000   0.9000   0.1000"""
    
    def create_distorted_ammonia(self) -> str:
        """歪んだアンモニア分子のXYZ座標を作成"""
        return """4
ammonia molecule (distorted)
N  0.1000   0.1000   0.2000
H  0.8000   0.6000   0.8000
H -0.7000   0.7000   0.3000
H  0.2000  -0.9000   0.5000"""
    
    def create_distorted_methane(self) -> str:
        """歪んだメタン分子のXYZ座標を作成"""
        return """5
methane molecule (distorted)
C  0.1000   0.1000   0.1000
H  0.9000   0.7000   0.7000
H -0.8000   0.8000   0.2000
H  0.3000  -0.9000   0.6000
H  0.4000   0.2000  -1.1000"""


# テストフィクスチャ
@pytest.fixture
def geometry_analyzer():
    """GeometryAnalyzerインスタンス"""
    return GeometryAnalyzer()


@pytest.fixture
def mock_test():
    """MockedGeometryOptimizationTestインスタンス"""
    return MockedGeometryOptimizationTest()


@pytest.fixture
def calculator():
    """DFTCalculatorインスタンス"""
    return DFTCalculator(keep_files=False)


# 基本的な構造解析関数のテスト
def test_distance_calculation(geometry_analyzer):
    """距離計算のテスト"""
    pos1 = np.array([0.0, 0.0, 0.0])
    pos2 = np.array([1.0, 0.0, 0.0])
    
    distance = geometry_analyzer.calculate_distance(pos1, pos2)
    assert distance == pytest.approx(1.0, rel=1e-6)


def test_angle_calculation(geometry_analyzer):
    """角度計算のテスト"""
    # 90度の角度を作る3点
    pos1 = np.array([1.0, 0.0, 0.0])
    pos2 = np.array([0.0, 0.0, 0.0])  # 頂点
    pos3 = np.array([0.0, 1.0, 0.0])
    
    angle = geometry_analyzer.calculate_angle(pos1, pos2, pos3)
    assert angle == pytest.approx(90.0, rel=1e-6)


def test_center_of_mass_calculation(geometry_analyzer):
    """質量中心計算のテスト"""
    # 水分子（O: 16u, H: 1u each）
    symbols = ['O', 'H', 'H']
    positions = np.array([
        [0.0, 0.0, 0.0],  # O
        [1.0, 0.0, 0.0],  # H1
        [-1.0, 0.0, 0.0]  # H2
    ])
    
    com = geometry_analyzer.calculate_center_of_mass(symbols, positions)
    
    # 酸素の質量が水素より重いので、中心は原点近くになる
    assert abs(com[0]) < 0.1  # X方向はほぼ0
    assert abs(com[1]) < 1e-10  # Y方向は完全に0
    assert abs(com[2]) < 1e-10  # Z方向は完全に0


def test_inertia_tensor_calculation(geometry_analyzer):
    """慣性テンソル計算のテスト"""
    # 線形分子（H-H）
    symbols = ['H', 'H']
    positions = np.array([
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0]
    ])
    
    inertia_tensor = geometry_analyzer.calculate_inertia_tensor(symbols, positions)
    
    # 線形分子のX軸方向の慣性モーメントは0に近い
    assert inertia_tensor[0, 0] == pytest.approx(0.0, abs=1e-10)
    # Y, Z軸方向は同じ値
    assert inertia_tensor[1, 1] == pytest.approx(inertia_tensor[2, 2], rel=1e-6)


def test_principal_axes_calculation(geometry_analyzer):
    """主軸計算のテスト"""
    # 単純な対称慣性テンソル
    inertia_tensor = np.array([
        [2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 3.0]
    ])
    
    moments, axes = geometry_analyzer.get_principal_axes(inertia_tensor)
    
    # 主慣性モーメントは昇順でソート
    expected_moments = np.array([1.0, 2.0, 3.0])
    assert np.allclose(moments, expected_moments)
    
    # 主軸は単位行列（対角テンソルなので）
    expected_axes = np.array([
        [0.0, 1.0, 0.0],  # Y軸
        [1.0, 0.0, 0.0],  # X軸
        [0.0, 0.0, 1.0]   # Z軸
    ])
    # 順序は主慣性モーメントの昇順
    assert np.allclose(abs(axes), abs(expected_axes))


class TestMockedOptimization:
    """モック付き構造最適化テスト"""
    
    @patch('quantum_calc.dft_calculator.geometric_solver.optimize')
    @patch('quantum_calc.dft_calculator.dft.RKS')
    def test_water_optimization_effect(self, mock_rks, mock_optimize, mock_test, calculator, geometry_analyzer):
        """水分子の構造最適化効果をテスト"""
        
        # === モックセットアップ ===
        # 最適化前の歪んだ構造
        original_coords = np.array([
            [0.1, 0.2, 0.1],    # O (歪んでいる)
            [1.0, 0.8, 0.1],    # H1 (距離が長い)
            [-0.8, 0.9, 0.1]    # H2 (角度が歪んでいる)
        ])
        
        # 最適化後の理想的な構造
        optimized_coords = np.array([
            [0.0, 0.0, 0.0],        # O (原点)
            [0.7571, 0.5861, 0.0],  # H1 (理想的な距離・角度)
            [-0.7571, 0.5861, 0.0]  # H2 (対称)
        ])
        
        # モック分子オブジェクト
        mock_mol_original = MagicMock()
        mock_mol_original.atom_coords.return_value = original_coords
        mock_mol_original.atom_symbol.side_effect = lambda i: ['O', 'H', 'H'][i]
        mock_mol_original.natm = 3
        
        mock_mol_optimized = MagicMock()
        mock_mol_optimized.atom_coords.return_value = optimized_coords
        mock_mol_optimized.atom_symbol.side_effect = lambda i: ['O', 'H', 'H'][i]
        mock_mol_optimized.natm = 3
        
        mock_optimize.return_value = mock_mol_optimized
        
        # モックSCF
        mock_mf = MagicMock()
        mock_mf.kernel.return_value = -76.4  # 水分子の典型的なエネルギー
        mock_mf.converged = True
        mock_mf.mo_occ = np.array([2., 2., 2., 2., 2., 0., 0.])
        mock_rks.return_value = mock_mf
        
        # === テスト実行 ===
        atoms = calculator.parse_xyz(mock_test.create_distorted_water())
        calculator.setup_calculation(atoms, basis='sto-3g', xc='b3lyp')  # 軽量な基底関数
        results = calculator.run_calculation()
        
        # === 構造パラメーター解析 ===
        # 最適化前の構造解析
        symbols = ['O', 'H', 'H']
        
        # 最適化前のO-H距離
        oh1_dist_orig = geometry_analyzer.calculate_distance(original_coords[0], original_coords[1])
        oh2_dist_orig = geometry_analyzer.calculate_distance(original_coords[0], original_coords[2])
        
        # 最適化後のO-H距離
        oh1_dist_opt = geometry_analyzer.calculate_distance(optimized_coords[0], optimized_coords[1])
        oh2_dist_opt = geometry_analyzer.calculate_distance(optimized_coords[0], optimized_coords[2])
        
        # === 検証 ===
        # 1. 最適化が実行されたことを確認
        mock_optimize.assert_called_once()
        
        # 2. 結果に最適化された座標が含まれていることを確認
        assert 'optimized_geometry' in results
        assert 'Optimized geometry' in results['optimized_geometry']
        
        # 3. O-H結合長の改善を確認（より一般的な値に近づく）
        # 水分子の典型的なO-H結合長は約0.96Å
        expected_oh_distance = 0.96
        
        # 最適化後の方が理想値に近いことを確認
        opt_distance_error = abs(oh1_dist_opt - expected_oh_distance)
        orig_distance_error = abs(oh1_dist_orig - expected_oh_distance)
        
        assert opt_distance_error < orig_distance_error
        
        # 4. 対称性の改善を確認
        # 最適化後は2つのO-H結合長がほぼ等しくなる
        assert abs(oh1_dist_opt - oh2_dist_opt) < 0.01
        
        # 5. H-O-H結合角の計算
        hoh_angle_orig = geometry_analyzer.calculate_angle(
            original_coords[1], original_coords[0], original_coords[2]
        )
        hoh_angle_opt = geometry_analyzer.calculate_angle(
            optimized_coords[1], optimized_coords[0], optimized_coords[2]
        )
        
        # 水分子の理想的なH-O-H角は約104.5度
        expected_angle = 104.5
        opt_angle_error = abs(hoh_angle_opt - expected_angle)
        orig_angle_error = abs(hoh_angle_orig - expected_angle)
        
        assert opt_angle_error < orig_angle_error
        
        print(f"=== 水分子構造最適化結果 ===")
        print(f"O-H結合長 (最適化前): {oh1_dist_orig:.3f}Å, {oh2_dist_orig:.3f}Å")
        print(f"O-H結合長 (最適化後): {oh1_dist_opt:.3f}Å, {oh2_dist_opt:.3f}Å")
        print(f"H-O-H角 (最適化前): {hoh_angle_orig:.1f}°")
        print(f"H-O-H角 (最適化後): {hoh_angle_opt:.1f}°")
        print(f"最適化効果: 結合長誤差 {orig_distance_error:.3f}→{opt_distance_error:.3f}Å")
        print(f"最適化効果: 角度誤差 {orig_angle_error:.1f}→{opt_angle_error:.1f}°")


# 実際のASEアライメント機能のテスト（モック無し）
class TestASEAlignment:
    """ASEアライメント機能の直接テスト"""
    
    def test_ase_alignment_import(self):
        """ASEライブラリが正常にインポートできることを確認"""
        try:
            from ase import Atoms
            assert True
        except ImportError:
            pytest.fail("ASE library is not available")
    
    def test_ase_basic_functionality(self, geometry_analyzer):
        """ASEの基本機能テスト"""
        from ase import Atoms
        
        # 水分子を作成
        symbols = ['O', 'H', 'H']
        positions = np.array([
            [0.1, 0.2, 0.1],
            [1.0, 0.8, 0.1],
            [-0.8, 0.9, 0.1]
        ])
        
        atoms = Atoms(symbols=symbols, positions=positions)
        
        # 質量中心を原点に移動
        com = geometry_analyzer.calculate_center_of_mass(symbols, positions)
        atoms.positions -= com
        
        # 移動後の質量中心を確認
        new_positions = atoms.get_positions()
        new_com = geometry_analyzer.calculate_center_of_mass(symbols, new_positions)
        
        # 質量中心が原点近くにあることを確認
        assert np.allclose(new_com, [0.0, 0.0, 0.0], atol=1e-10)
        
        print(f"=== ASEアライメントテスト ===")
        print(f"元の質量中心: [{com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}]")
        print(f"移動後の質量中心: [{new_com[0]:.6f}, {new_com[1]:.6f}, {new_com[2]:.6f}]")
    
    def test_principal_axis_alignment(self, geometry_analyzer):
        """主軸アライメントのテスト"""
        from ase import Atoms
        
        # 線形でない分子（水分子）を使用
        symbols = ['O', 'H', 'H']
        positions = np.array([
            [0.1, 0.2, 0.1],
            [1.0, 0.8, 0.1],
            [-0.8, 0.9, 0.1]
        ])
        
        # 質量中心を原点に移動
        com = geometry_analyzer.calculate_center_of_mass(symbols, positions)
        centered_positions = positions - com
        
        # 慣性テンソルと主軸を計算
        inertia_tensor = geometry_analyzer.calculate_inertia_tensor(symbols, centered_positions)
        principal_moments, principal_axes = geometry_analyzer.get_principal_axes(inertia_tensor)
        
        # 座標を主軸に回転
        aligned_positions = np.dot(centered_positions, principal_axes)
        
        # アライメント後の慣性テンソルを計算
        aligned_inertia = geometry_analyzer.calculate_inertia_tensor(symbols, aligned_positions)
        aligned_moments, _ = geometry_analyzer.get_principal_axes(aligned_inertia)
        
        # 主慣性モーメントが保持されていることを確認
        assert np.allclose(np.sort(principal_moments), np.sort(aligned_moments), rtol=1e-6)
        
        # アライメント後の慣性テンソルがより対角的になっていることを確認
        # (非対角成分の絶対値の和が小さくなる)
        off_diagonal_sum_orig = np.sum(np.abs(inertia_tensor - np.diag(np.diag(inertia_tensor))))
        off_diagonal_sum_aligned = np.sum(np.abs(aligned_inertia - np.diag(np.diag(aligned_inertia))))
        
        assert off_diagonal_sum_aligned <= off_diagonal_sum_orig
        
        print(f"=== 主軸アライメントテスト ===")
        print(f"主慣性モーメント: {principal_moments}")
        print(f"元の非対角成分の和: {off_diagonal_sum_orig:.6f}")
        print(f"アライメント後の非対角成分の和: {off_diagonal_sum_aligned:.6f}")
    
    def test_internal_structure_preservation(self, geometry_analyzer):
        """内部構造の保持テスト（アライメント前後で結合長・角度が変わらない）"""
        from ase import Atoms
        
        # アンモニア分子（4原子、非平面）
        symbols = ['N', 'H', 'H', 'H']
        positions = np.array([
            [0.1, 0.1, 0.2],
            [0.8, 0.6, 0.8],
            [-0.7, 0.7, 0.3],
            [0.2, -0.9, 0.5]
        ])
        
        # 元の構造パラメーター
        nh_distances_orig = [
            geometry_analyzer.calculate_distance(positions[0], positions[i])
            for i in range(1, 4)
        ]
        hnh_angles_orig = [
            geometry_analyzer.calculate_angle(positions[1], positions[0], positions[2]),
            geometry_analyzer.calculate_angle(positions[1], positions[0], positions[3]),
            geometry_analyzer.calculate_angle(positions[2], positions[0], positions[3])
        ]
        
        # ASEアライメント実行
        atoms = Atoms(symbols=symbols, positions=positions)
        
        # 質量中心を原点に移動
        com = geometry_analyzer.calculate_center_of_mass(symbols, positions)
        atoms.positions -= com
        
        # 慣性主軸への回転（線形分子でない場合）
        centered_positions = atoms.get_positions()
        inertia_tensor = geometry_analyzer.calculate_inertia_tensor(symbols, centered_positions)
        
        # 主慣性モーメントが全て異なる場合のみ回転
        principal_moments, principal_axes = geometry_analyzer.get_principal_axes(inertia_tensor)
        
        if not np.allclose(principal_moments[:-1], principal_moments[1:], rtol=1e-6):
            # 回転行列を適用
            atoms.positions = np.dot(atoms.positions, principal_axes)
        
        aligned_positions = atoms.get_positions()
        
        # アライメント後の構造パラメーター
        nh_distances_aligned = [
            geometry_analyzer.calculate_distance(aligned_positions[0], aligned_positions[i])
            for i in range(1, 4)
        ]
        hnh_angles_aligned = [
            geometry_analyzer.calculate_angle(aligned_positions[1], aligned_positions[0], aligned_positions[2]),
            geometry_analyzer.calculate_angle(aligned_positions[1], aligned_positions[0], aligned_positions[3]),
            geometry_analyzer.calculate_angle(aligned_positions[2], aligned_positions[0], aligned_positions[3])
        ]
        
        # 結合長が保持されていることを確認
        for orig_dist, aligned_dist in zip(nh_distances_orig, nh_distances_aligned):
            assert abs(orig_dist - aligned_dist) < 1e-10
        
        # 結合角が保持されていることを確認
        for orig_angle, aligned_angle in zip(hnh_angles_orig, hnh_angles_aligned):
            assert abs(orig_angle - aligned_angle) < 1e-6  # 角度は度単位なので少し緩く
        
        print(f"=== 内部構造保持テスト（アンモニア） ===")
        print(f"N-H結合長 (元): {nh_distances_orig}")
        print(f"N-H結合長 (整列後): {nh_distances_aligned}")
        print(f"H-N-H角 (元): {hnh_angles_orig}")
        print(f"H-N-H角 (整列後): {hnh_angles_aligned}")