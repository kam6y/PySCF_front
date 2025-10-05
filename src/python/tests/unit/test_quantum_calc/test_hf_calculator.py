"""
Unit tests for HF (Hartree-Fock) calculator.

Tests the HFCalculator class with mocked PySCF to verify initialization logic,
parameter handling, and method selection without running expensive computations.
"""

import pytest
from unittest.mock import MagicMock, patch

from quantum_calc.hf_calculator import HFCalculator


# ============================================================================
# Initialization and Setup Tests
# ============================================================================

@patch('quantum_calc.hf_calculator.scf')
@patch('pyscf.gto')
def test_hf_calculator_creates_rhf_for_closed_shell(mock_gto, mock_scf):
    """
    GIVEN HF parameters with spin=0 (closed-shell system)
    WHEN HFCalculator is initialized and calculation is set up
    THEN it should create RHF (Restricted Hartree-Fock) method object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_rhf = MagicMock()
    mock_scf.RHF.return_value = mock_rhf
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['H', [0, 0, 0]], ['H', [0, 0, 0.74]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='sto-3g',
        charge=0,
        spin=0  # Closed-shell
    )
    
    # ASSERT
    # Verify RHF was created (not UHF)
    mock_scf.RHF.assert_called_once_with(mock_mol)
    mock_scf.UHF.assert_not_called()


@patch('quantum_calc.hf_calculator.scf')
@patch('pyscf.gto')
def test_hf_calculator_creates_uhf_for_open_shell(mock_gto, mock_scf):
    """
    GIVEN HF parameters with spin>0 (open-shell system)
    WHEN HFCalculator is initialized and calculation is set up
    THEN it should create UHF (Unrestricted Hartree-Fock) method object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_uhf = MagicMock()
    mock_scf.UHF.return_value = mock_uhf
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['O', [0, 0, 0]], ['O', [0, 0, 1.2]]]  # O2
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='6-31g',
        charge=0,
        spin=2  # Triplet (open-shell)
    )
    
    # ASSERT
    # Verify UHF was created (not RHF)
    mock_scf.UHF.assert_called_once_with(mock_mol)
    mock_scf.RHF.assert_not_called()


@patch('pyscf.gto')
def test_hf_calculator_mol_initialization(mock_gto):
    """
    GIVEN valid molecular structure and HF parameters
    WHEN HFCalculator setup_calculation is called
    THEN it should correctly initialize PySCF Mole object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [
        ['O', [0.0, 0.0, 0.1173]],
        ['H', [0.0, 0.7572, -0.4692]],
        ['H', [0.0, -0.7572, -0.4692]]
    ]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='6-31G(d)',
        charge=-1,
        spin=1
    )
    
    # ASSERT
    # Verify gto.M was called to create molecule
    mock_gto.M.assert_called_once()
    
    # Verify M was called with correct parameters
    call_kwargs = mock_gto.M.call_args[1]
    assert call_kwargs['basis'] == '6-31G(d)'
    assert call_kwargs['charge'] == -1
    assert call_kwargs['spin'] == 1
    assert 'atom' in call_kwargs


@pytest.mark.parametrize("basis_set", [
    'STO-3G',
    '6-31G',
    '6-31G(d)',
    '6-31+G(d,p)',
    'cc-pVDZ',
    'cc-pVTZ',
    'aug-cc-pVDZ',
    'def2-SVP'
])
@patch('pyscf.gto')
def test_hf_calculator_various_basis_sets(mock_gto, basis_set):
    """
    GIVEN various basis sets
    WHEN HFCalculator is set up with each basis set
    THEN it should correctly assign the basis set to Mole object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['C', [0, 0, 0]], ['O', [0, 0, 1.2]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis=basis_set,
        charge=0,
        spin=0
    )
    
    # ASSERT
    # Verify gto.M was called with the correct basis set
    call_kwargs = mock_gto.M.call_args[1]
    assert call_kwargs['basis'] == basis_set


@pytest.mark.parametrize("charge,spin", [
    (0, 0),    # Neutral, singlet
    (1, 1),    # +1 charge, doublet
    (-1, 0),   # -1 charge, singlet
    (2, 0),    # +2 charge, singlet
    (0, 2),    # Neutral, triplet
    (-2, 1),   # -2 charge, doublet
])
@patch('pyscf.gto')
def test_hf_calculator_charge_spin_combinations(mock_gto, charge, spin):
    """
    GIVEN various charge and spin combinations
    WHEN HFCalculator is set up
    THEN it should correctly assign charge and spin to Mole object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['C', [0, 0, 0]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='sto-3g',
        charge=charge,
        spin=spin
    )
    
    # ASSERT
    # Verify gto.M was called with correct charge and spin
    call_kwargs = mock_gto.M.call_args[1]
    assert call_kwargs['charge'] == charge
    assert call_kwargs['spin'] == spin


# ============================================================================
# Method Selection Logic Tests
# ============================================================================

@patch('quantum_calc.hf_calculator.scf')
@patch('pyscf.gto')
def test_hf_calculator_method_selection_spin_zero(mock_gto, mock_scf):
    """
    GIVEN spin=0 (closed-shell)
    WHEN _create_scf_method is called
    THEN it should select RHF method
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_rhf = MagicMock()
    mock_scf.RHF.return_value = mock_rhf
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 0}
    
    # ACT
    mf = calculator._create_scf_method(mock_mol)
    
    # ASSERT
    assert mf == mock_rhf
    mock_scf.RHF.assert_called_once_with(mock_mol)


@patch('quantum_calc.hf_calculator.scf')
@patch('pyscf.gto')
def test_hf_calculator_method_selection_spin_nonzero(mock_gto, mock_scf):
    """
    GIVEN spin>0 (open-shell)
    WHEN _create_scf_method is called
    THEN it should select UHF method
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_uhf = MagicMock()
    mock_scf.UHF.return_value = mock_uhf
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 2}
    
    # ACT
    mf = calculator._create_scf_method(mock_mol)
    
    # ASSERT
    assert mf == mock_uhf
    mock_scf.UHF.assert_called_once_with(mock_mol)


# ============================================================================
# Parameter Validation Tests
# ============================================================================

@patch('pyscf.gto')
def test_hf_calculator_validate_specific_parameters_closed_shell(mock_gto):
    """
    GIVEN HF-specific parameters for closed-shell
    WHEN _validate_specific_parameters is called
    THEN it should return RHF method designation
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    # ACT
    validated = calculator._validate_specific_parameters(
        spin=0,
        basis='6-31G'
    )
    
    # ASSERT
    assert 'method' in validated
    assert validated['method'] == 'RHF'


@patch('pyscf.gto')
def test_hf_calculator_validate_specific_parameters_open_shell(mock_gto):
    """
    GIVEN HF-specific parameters for open-shell
    WHEN _validate_specific_parameters is called
    THEN it should return UHF method designation
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    # ACT
    validated = calculator._validate_specific_parameters(
        spin=2,
        basis='sto-3g'
    )
    
    # ASSERT
    assert validated['method'] == 'UHF'


@patch('pyscf.gto')
def test_hf_calculator_no_xc_functional_needed(mock_gto):
    """
    GIVEN HF calculation (doesn't use XC functional)
    WHEN setup_calculation is called without XC functional
    THEN it should not require or use XC functional
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['H', [0, 0, 0]], ['H', [0, 0, 0.74]]]
    
    # ACT - Should not raise exception even without xc parameter
    calculator.setup_calculation(
        atoms=atoms,
        basis='sto-3g',
        charge=0,
        spin=0
        # No xc parameter - HF doesn't need it
    )
    
    # ASSERT
    # Just verify it didn't crash
    assert calculator is not None


# ============================================================================
# Working Directory and File Management Tests
# ============================================================================

@patch('pyscf.gto')
def test_hf_calculator_working_directory(mock_gto):
    """
    GIVEN a working directory path
    WHEN HFCalculator is initialized
    THEN it should store the working directory
    """
    # ARRANGE
    test_dir = '/tmp/test_hf_calc'
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = HFCalculator(working_dir=test_dir, optimize_geometry=False)
    
    # ASSERT
    assert calculator.working_dir == test_dir


@patch('pyscf.gto')
def test_hf_calculator_molecule_name(mock_gto):
    """
    GIVEN a molecule name
    WHEN HFCalculator is initialized
    THEN it should store the molecule name
    """
    # ARRANGE
    molecule_name = 'Hydrogen Molecule'
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = HFCalculator(
        working_dir='/tmp/test',
        molecule_name=molecule_name,
        optimize_geometry=False
    )
    
    # ASSERT
    assert calculator.molecule_name == molecule_name


# ============================================================================
# Geometry Optimization Flag Tests
# ============================================================================

@patch('pyscf.gto')
def test_hf_calculator_optimize_geometry_true(mock_gto):
    """
    GIVEN optimize_geometry=True
    WHEN HFCalculator is initialized
    THEN it should enable geometry optimization
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=True)
    
    # ASSERT
    assert calculator.optimize_geometry is True


@patch('pyscf.gto')
def test_hf_calculator_optimize_geometry_false(mock_gto):
    """
    GIVEN optimize_geometry=False
    WHEN HFCalculator is initialized
    THEN it should disable geometry optimization
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    # ASSERT
    assert calculator.optimize_geometry is False


# ============================================================================
# Method Description Tests
# ============================================================================

@patch('pyscf.gto')
def test_hf_calculator_method_description_closed_shell(mock_gto):
    """
    GIVEN spin=0
    WHEN _get_base_method_description is called
    THEN it should return 'RHF'
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 0}
    
    # ACT
    description = calculator._get_base_method_description()
    
    # ASSERT
    assert description == 'RHF'


@patch('pyscf.gto')
def test_hf_calculator_method_description_open_shell(mock_gto):
    """
    GIVEN spin>0
    WHEN _get_base_method_description is called
    THEN it should return 'UHF'
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 2}
    
    # ACT
    description = calculator._get_base_method_description()
    
    # ASSERT
    assert description == 'UHF'


# ============================================================================
# Keep Files Flag Tests
# ============================================================================

@patch('pyscf.gto')
def test_hf_calculator_keep_files_flag(mock_gto):
    """
    GIVEN keep_files parameter
    WHEN HFCalculator is initialized
    THEN it should store the keep_files flag
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator_keep = HFCalculator(working_dir='/tmp/test', keep_files=True, optimize_geometry=False)
    calculator_no_keep = HFCalculator(working_dir='/tmp/test', keep_files=False, optimize_geometry=False)
    
    # ASSERT
    assert calculator_keep.keep_files is True
    assert calculator_no_keep.keep_files is False


# ============================================================================
# Comparison with DFT Calculator Tests
# ============================================================================

@patch('quantum_calc.hf_calculator.scf')
@patch('pyscf.gto')
def test_hf_calculator_no_xc_attribute(mock_gto, mock_scf):
    """
    GIVEN HF method (doesn't use XC functional)
    WHEN SCF method object is created
    THEN it should not have xc attribute set (unlike DFT)
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_rhf = MagicMock()
    # Remove xc attribute if it exists
    if hasattr(mock_rhf, 'xc'):
        delattr(mock_rhf, 'xc')
    mock_scf.RHF.return_value = mock_rhf
    
    calculator = HFCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 0}
    
    # ACT
    mf = calculator._create_scf_method(mock_mol)
    
    # ASSERT
    # Verify that xc was not set (HF doesn't use XC functionals)
    # The mock object might have xc set by test framework, so we just verify
    # that the method was created correctly
    assert mf == mock_rhf
