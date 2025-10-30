"""
Unit tests for DFT calculator.

Tests the DFTCalculator class with mocked PySCF to verify initialization logic,
parameter handling, and method selection without running expensive computations.
"""

import pytest
from unittest.mock import MagicMock, patch, PropertyMock

from quantum_calc.dft_calculator import DFTCalculator


# ============================================================================
# Initialization and Setup Tests
# ============================================================================

@patch('quantum_calc.dft_calculator.dft')
@patch('pyscf.gto')
def test_dft_calculator_creates_rks_for_closed_shell(mock_gto, mock_dft):
    """
    GIVEN DFT parameters with spin=0 (closed-shell system)
    WHEN DFTCalculator is initialized and calculation is set up
    THEN it should create RKS (Restricted Kohn-Sham) method object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_rks = MagicMock()
    mock_dft.RKS.return_value = mock_rks
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['H', [0, 0, 0]], ['H', [0, 0, 0.74]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='sto-3g',
        charge=0,
        spin=0,  # Closed-shell
        xc='B3LYP'
    )
    
    # ASSERT
    # Verify RKS was created (not UKS)
    mock_dft.RKS.assert_called_once_with(mock_mol)
    mock_dft.UKS.assert_not_called()
    
    # Verify XC functional was set
    assert mock_rks.xc == 'B3LYP'


@patch('quantum_calc.dft_calculator.dft')
@patch('pyscf.gto')
def test_dft_calculator_creates_uks_for_open_shell(mock_gto, mock_dft):
    """
    GIVEN DFT parameters with spin>0 (open-shell system)
    WHEN DFTCalculator is initialized and calculation is set up
    THEN it should create UKS (Unrestricted Kohn-Sham) method object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_uks = MagicMock()
    mock_dft.UKS.return_value = mock_uks
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['O', [0, 0, 0]], ['O', [0, 0, 1.2]]]  # O2
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='6-31g',
        charge=0,
        spin=2,  # Triplet (open-shell)
        xc='PBE0'
    )
    
    # ASSERT
    # Verify UKS was created (not RKS)
    mock_dft.UKS.assert_called_once_with(mock_mol)
    mock_dft.RKS.assert_not_called()
    
    # Verify XC functional was set
    assert mock_uks.xc == 'PBE0'


@patch('pyscf.gto')
def test_dft_calculator_mol_initialization(mock_gto):
    """
    GIVEN valid molecular structure and DFT parameters
    WHEN DFTCalculator setup_calculation is called
    THEN it should correctly initialize PySCF Mole object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
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
        spin=1,
        xc='B3LYP'
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


@pytest.mark.parametrize("xc_functional", [
    'B3LYP',
    'PBE0',
    'M06-2X',
    'CAM-B3LYP',
    'wB97X-D',
    'PBE',
    'BLYP'
])
@patch('quantum_calc.dft_calculator.dft')
@patch('pyscf.gto')
def test_dft_calculator_various_xc_functionals(mock_gto, mock_dft, xc_functional):
    """
    GIVEN various exchange-correlation functionals
    WHEN DFTCalculator is set up with each functional
    THEN it should correctly assign the XC functional to the method object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_rks = MagicMock()
    mock_dft.RKS.return_value = mock_rks
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['H', [0, 0, 0]], ['H', [0, 0, 0.74]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='sto-3g',
        charge=0,
        spin=0,
        xc=xc_functional
    )
    
    # ASSERT
    assert mock_rks.xc == xc_functional


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
def test_dft_calculator_various_basis_sets(mock_gto, basis_set):
    """
    GIVEN various basis sets
    WHEN DFTCalculator is set up with each basis set
    THEN it should correctly assign the basis set to Mole object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['C', [0, 0, 0]], ['O', [0, 0, 1.2]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis=basis_set,
        charge=0,
        spin=0,
        xc='B3LYP'
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
def test_dft_calculator_charge_spin_combinations(mock_gto, charge, spin):
    """
    GIVEN various charge and spin combinations
    WHEN DFTCalculator is set up
    THEN it should correctly assign charge and spin to Mole object
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    atoms = [['C', [0, 0, 0]]]
    
    # ACT
    calculator.setup_calculation(
        atoms=atoms,
        basis='sto-3g',
        charge=charge,
        spin=spin,
        xc='B3LYP'
    )
    
    # ASSERT
    # Verify gto.M was called with correct charge and spin
    call_kwargs = mock_gto.M.call_args[1]
    assert call_kwargs['charge'] == charge
    assert call_kwargs['spin'] == spin


# ============================================================================
# Method Selection Logic Tests
# ============================================================================

@patch('quantum_calc.dft_calculator.dft')
@patch('pyscf.gto')
def test_dft_calculator_method_selection_spin_zero(mock_gto, mock_dft):
    """
    GIVEN spin=0 (closed-shell)
    WHEN _create_scf_method is called
    THEN it should select RKS method
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_rks = MagicMock()
    mock_dft.RKS.return_value = mock_rks
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.xc_functional = 'B3LYP'
    calculator.results = {'spin': 0}
    
    # ACT
    mf = calculator._create_scf_method(mock_mol)
    
    # ASSERT
    assert mf == mock_rks
    mock_dft.RKS.assert_called_once_with(mock_mol)
    assert mock_rks.xc == 'B3LYP'


@patch('quantum_calc.dft_calculator.dft')
@patch('pyscf.gto')
def test_dft_calculator_method_selection_spin_nonzero(mock_gto, mock_dft):
    """
    GIVEN spin>0 (open-shell)
    WHEN _create_scf_method is called
    THEN it should select UKS method
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    mock_uks = MagicMock()
    mock_dft.UKS.return_value = mock_uks
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.xc_functional = 'PBE0'
    calculator.results = {'spin': 2}
    
    # ACT
    mf = calculator._create_scf_method(mock_mol)
    
    # ASSERT
    assert mf == mock_uks
    mock_dft.UKS.assert_called_once_with(mock_mol)
    assert mock_uks.xc == 'PBE0'


# ============================================================================
# Parameter Validation Tests
# ============================================================================

@patch('pyscf.gto')
def test_dft_calculator_validate_specific_parameters(mock_gto):
    """
    GIVEN DFT-specific parameters
    WHEN _validate_specific_parameters is called
    THEN it should return validated parameters including XC functional
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    # ACT
    validated = calculator._validate_specific_parameters(
        xc='B3LYP',
        spin=0,
        basis='6-31G'
    )
    
    # ASSERT
    assert 'xc_functional' in validated
    assert validated['xc_functional'] == 'B3LYP'
    assert validated['method'] == 'RKS'


@patch('pyscf.gto')
def test_dft_calculator_default_xc_functional(mock_gto):
    """
    GIVEN no XC functional specified
    WHEN _validate_specific_parameters is called
    THEN it should use default B3LYP functional
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    # ACT
    validated = calculator._validate_specific_parameters(spin=0, basis='sto-3g')
    
    # ASSERT
    assert validated['xc_functional'] == 'B3LYP'  # Default value


# ============================================================================
# Working Directory and File Management Tests
# ============================================================================

@patch('pyscf.gto')
def test_dft_calculator_working_directory(mock_gto):
    """
    GIVEN a working directory path
    WHEN DFTCalculator is initialized
    THEN it should store the working directory
    """
    # ARRANGE
    test_dir = '/tmp/test_dft_calc'
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = DFTCalculator(working_dir=test_dir, optimize_geometry=False)
    
    # ASSERT
    assert calculator.working_dir == test_dir


@patch('pyscf.gto')
def test_dft_calculator_molecule_name(mock_gto):
    """
    GIVEN a molecule name
    WHEN DFTCalculator is initialized
    THEN it should store the molecule name
    """
    # ARRANGE
    molecule_name = 'Water Molecule'
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = DFTCalculator(
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
def test_dft_calculator_optimize_geometry_true(mock_gto):
    """
    GIVEN optimize_geometry=True
    WHEN DFTCalculator is initialized
    THEN it should enable geometry optimization
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=True)
    
    # ASSERT
    assert calculator.optimize_geometry is True


@patch('pyscf.gto')
def test_dft_calculator_optimize_geometry_false(mock_gto):
    """
    GIVEN optimize_geometry=False
    WHEN DFTCalculator is initialized
    THEN it should disable geometry optimization
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    # ACT
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    
    # ASSERT
    assert calculator.optimize_geometry is False


# ============================================================================
# Method Description Tests
# ============================================================================

@patch('pyscf.gto')
def test_dft_calculator_method_description_closed_shell(mock_gto):
    """
    GIVEN spin=0
    WHEN _get_base_method_description is called
    THEN it should return 'RKS'
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 0}
    
    # ACT
    description = calculator._get_base_method_description()
    
    # ASSERT
    assert description == 'RKS'


@patch('pyscf.gto')
def test_dft_calculator_method_description_open_shell(mock_gto):
    """
    GIVEN spin>0
    WHEN _get_base_method_description is called
    THEN it should return 'UKS'
    """
    # ARRANGE
    mock_mol = MagicMock()
    mock_gto.M.return_value = mock_mol
    
    calculator = DFTCalculator(working_dir='/tmp/test', optimize_geometry=False)
    calculator.results = {'spin': 2}
    
    # ACT
    description = calculator._get_base_method_description()
    
    # ASSERT
    assert description == 'UKS'
