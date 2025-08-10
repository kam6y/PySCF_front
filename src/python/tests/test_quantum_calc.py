import pytest
import numpy as np
from unittest.mock import patch, MagicMock

from quantum_calc.dft_calculator import DFTCalculator
from quantum_calc.exceptions import InputError, ConvergenceError, GeometryError

@pytest.fixture
def water_xyz():
    return """3
water molecule
O  0.0000   0.0000   0.0000
H  0.7571   0.5861   0.0000  
H -0.7571   0.5861   0.0000"""

@pytest.fixture
def calculator():
    # keep_files=Trueにすることで、テスト実行時に一時ディレクトリが削除されないようにし、
    # 失敗時のデバッグがしやすくなる（通常はFalseで良い）
    return DFTCalculator(keep_files=False)

def test_parse_valid_xyz(calculator, water_xyz):
    """Test parsing a valid XYZ string."""
    atoms = calculator.parse_xyz(water_xyz)
    assert len(atoms) == 3
    assert atoms[0][0] == 'O'
    assert atoms[1][1] == [0.7571, 0.5861, 0.0]

def test_parse_invalid_xyz(calculator):
    """Test parsing invalid XYZ strings."""
    with pytest.raises(ValueError, match="insufficient lines"):
        calculator.parse_xyz("1\nC 0 0 0")
    with pytest.raises(ValueError, match="invalid coordinates"):
        calculator.parse_xyz("1\nComment\nC 0 0 a")

@patch('quantum_calc.dft_calculator.gto')
def test_setup_calculation(mock_gto, calculator, water_xyz):
    """Test the setup of a DFT calculation."""
    atoms = calculator.parse_xyz(water_xyz)
    calculator.setup_calculation(
        atoms, basis='6-31g', xc='b3lyp', charge=0, spin=0, max_cycle=100
    )
    
    assert calculator.mol is not None
    assert calculator.mf is not None
    mock_gto.M.assert_called_once()
    assert calculator.mf.xc == 'b3lyp'
    assert calculator.mf.max_cycle == 100

def test_setup_with_bad_geometry(calculator):
    """Test setup with invalid atom coordinates."""
    with pytest.raises(GeometryError):
        atoms = [['C', [0.0, 0.0]]] # 2D coordinates instead of 3D
        calculator.setup_calculation(atoms)

@patch('quantum_calc.dft_calculator.geometric_solver.optimize')
def test_run_successful_calculation(mock_optimize, calculator, water_xyz):
    """Test a successful calculation run, mocking PySCF solvers."""
    # --- Mock Setup ---
    mock_mol = MagicMock()
    mock_mol.atom_coords.return_value = np.array([[0,0,0], [0,0,1], [0,1,0]])
    mock_mol.atom_symbol.side_effect = ['O', 'H', 'H']
    mock_mol.natm = 3
    mock_optimize.return_value = mock_mol

    mock_mf = MagicMock()
    mock_mf.kernel.return_value = -76.4 # Mocked SCF energy
    mock_mf.converged = True
    mock_mf.mo_occ = np.array([2., 2., 2., 2., 2., 0., 0.]) # 5 occupied, 2 virtual
    
    # --- Test Execution ---
    with patch('quantum_calc.dft_calculator.dft.RKS', return_value=mock_mf):
        atoms = calculator.parse_xyz(water_xyz)
        calculator.setup_calculation(atoms)
        results = calculator.run_calculation()

    # --- Assertions ---
    mock_optimize.assert_called_once()
    assert mock_mf.kernel.call_count == 1
    assert results['converged'] is True
    assert results['scf_energy'] == pytest.approx(-76.4)
    assert results['homo_index'] == 4 # Index of the last occupied orbital
    assert results['lumo_index'] == 5 # Index of the first virtual orbital
    assert 'Optimized geometry' in results['optimized_geometry']

@patch('quantum_calc.dft_calculator.geometric_solver.optimize')
def test_run_convergence_error(mock_optimize, calculator, water_xyz):
    """Test handling of an SCF convergence failure."""
    mock_mol = MagicMock()
    mock_optimize.return_value = mock_mol

    mock_mf = MagicMock()
    mock_mf.kernel.return_value = -75.0
    mock_mf.converged = False # Simulate convergence failure

    with patch('quantum_calc.dft_calculator.dft.RKS', return_value=mock_mf):
        atoms = calculator.parse_xyz(water_xyz)
        calculator.setup_calculation(atoms)
        with pytest.raises(ConvergenceError, match="SCF calculation failed to converge"):
            calculator.run_calculation()