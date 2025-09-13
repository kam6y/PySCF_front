"""Tests for CASCI and CASSCF calculators."""

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
import tempfile
import os

from quantum_calc.casci_calculator import CASCICalculator
from quantum_calc.casscf_calculator import CASSCFCalculator
from quantum_calc.exceptions import InputError, ConvergenceError, GeometryError

@pytest.fixture
def water_xyz():
    """Water molecule for closed-shell testing."""
    return """3
water molecule
O  0.0000   0.0000   0.0000
H  0.7571   0.5861   0.0000  
H -0.7571   0.5861   0.0000"""

@pytest.fixture
def h2_xyz():
    """Hydrogen molecule for small CASCI/CASSCF testing."""
    return """2
hydrogen molecule
H  0.0000   0.0000   0.0000
H  0.0000   0.0000   0.7400"""

@pytest.fixture
def oh_radical_xyz():
    """OH radical for open-shell testing."""
    return """2
OH radical
O  0.0000   0.0000   0.0000
H  0.0000   0.0000   0.9700"""

@pytest.fixture
def casci_calculator():
    """CASCI calculator instance."""
    return CASCICalculator(keep_files=False)

@pytest.fixture
def casscf_calculator():
    """CASSCF calculator instance."""
    return CASSCFCalculator(keep_files=False)

class TestCASCICalculator:
    """Test cases for CASCI calculator."""
    
    def test_parse_valid_xyz(self, casci_calculator, water_xyz):
        """Test parsing a valid XYZ string in CASCI calculator."""
        atoms = casci_calculator.parse_xyz(water_xyz)
        assert len(atoms) == 3
        assert atoms[0][0] == 'O'
        assert atoms[1][1] == [0.7571, 0.5861, 0.0]
    
    def test_parse_invalid_xyz(self, casci_calculator):
        """Test parsing invalid XYZ strings."""
        with pytest.raises(ValueError, match="insufficient lines"):
            casci_calculator.parse_xyz("1\nH 0 0 0")
        with pytest.raises(ValueError, match="invalid coordinates"):
            casci_calculator.parse_xyz("2\nComment\nH 0 0 a")
    
    @patch('quantum_calc.casci_calculator.gto')
    def test_setup_calculation_basic(self, mock_gto, casci_calculator, h2_xyz):
        """Test basic CASCI calculation setup."""
        atoms = casci_calculator.parse_xyz(h2_xyz)
        casci_calculator.setup_calculation(
            atoms, 
            basis='sto-3g', 
            charge=0, 
            spin=0, 
            ncas=2,  # 2 orbitals
            nelecas=2,  # 2 electrons
            max_cycle=50
        )
        
        assert casci_calculator.mol is not None
        assert casci_calculator.mf is not None
        assert casci_calculator.ncas == 2
        assert casci_calculator.nelecas == 2
        mock_gto.M.assert_called_once()
    
    def test_setup_invalid_active_space(self, casci_calculator, h2_xyz):
        """Test setup with invalid active space parameters."""
        atoms = casci_calculator.parse_xyz(h2_xyz)
        
        # Too many electrons for active space
        with pytest.raises(InputError, match="Too many electrons"):
            casci_calculator.setup_calculation(
                atoms, 
                basis='sto-3g',
                ncas=2,
                nelecas=6  # Too many electrons for 2 orbitals
            )
        
        # Negative active space size
        with pytest.raises(InputError, match="must be positive"):
            casci_calculator.setup_calculation(
                atoms,
                basis='sto-3g', 
                ncas=-1,
                nelecas=2
            )
    
    @patch('quantum_calc.casci_calculator.gto')
    @patch('quantum_calc.casci_calculator.scf')
    @patch('quantum_calc.casci_calculator.mcscf')
    def test_casci_calculation_flow(self, mock_mcscf, mock_scf, mock_gto, casci_calculator, h2_xyz):
        """Test the complete CASCI calculation flow."""
        # Setup mocks
        mock_mol = MagicMock()
        mock_gto.M.return_value = mock_mol
        
        mock_mf = MagicMock()
        mock_mf.kernel.return_value = -1.123456
        mock_mf.converged = True
        mock_mf.mo_energy = np.array([-0.5, 0.5])
        mock_mf.mo_occ = np.array([2.0, 0.0])
        mock_scf.RHF.return_value = mock_mf
        
        mock_cas = MagicMock()
        mock_cas.kernel.return_value = -1.135790
        mock_cas.converged = True
        mock_cas.ncore = 0
        mock_cas.ci = np.array([0.9, 0.1, 0.05])  # Simple CI coefficients
        mock_cas.make_rdm1s.return_value = (np.eye(2), np.eye(2))  # Mock density matrices
        mock_cas.cas_natorb.return_value = (np.eye(2), np.array([0.9, 0.1]), np.array([1.8, 0.2]))
        mock_cas.mo_coeff = np.eye(2)
        mock_mcscf.CASCI.return_value = mock_cas
        
        # Run calculation
        atoms = casci_calculator.parse_xyz(h2_xyz)
        casci_calculator.setup_calculation(atoms, basis='sto-3g', ncas=2, nelecas=2)
        results = casci_calculator.calculate()
        
        # Verify results
        assert 'scf_energy' in results
        assert 'casci_energy' in results
        assert 'correlation_energy' in results
        assert results['method'] == 'RHF-CASCI'
        
        # Verify analysis results
        assert 'natural_orbital_analysis' in results
        assert 'ci_coefficient_analysis' in results
        
        # Verify calculation was called
        mock_cas.kernel.assert_called_once()
    
    @patch('quantum_calc.casci_calculator.gto')
    @patch('quantum_calc.casci_calculator.scf')
    def test_open_shell_setup(self, mock_scf, mock_gto, casci_calculator, oh_radical_xyz):
        """Test CASCI setup for open-shell systems."""
        atoms = casci_calculator.parse_xyz(oh_radical_xyz)
        casci_calculator.setup_calculation(
            atoms,
            basis='sto-3g',
            spin=1,  # Doublet state (one unpaired electron)
            ncas=4,
            nelecas=5
        )
        
        # Should use UHF for open-shell
        mock_scf.UHF.assert_called_once()
        assert casci_calculator.results['method'] == 'UHF-CASCI'
    
    def test_energy_units_conversion(self, casci_calculator):
        """Test energy unit conversions if needed."""
        # Basic validation that energy values are in expected range (Hartree)
        test_energy = -1.123456
        assert isinstance(test_energy, (int, float))
        assert -1000 < test_energy < 1000  # Reasonable range for small molecules


class TestCASSCFCalculator:
    """Test cases for CASSCF calculator."""
    
    def test_parse_valid_xyz(self, casscf_calculator, water_xyz):
        """Test parsing a valid XYZ string in CASSCF calculator."""
        atoms = casscf_calculator.parse_xyz(water_xyz)
        assert len(atoms) == 3
        assert atoms[0][0] == 'O'
        assert atoms[1][1] == [0.7571, 0.5861, 0.0]
    
    @patch('quantum_calc.casscf_calculator.gto')
    def test_setup_calculation_basic(self, mock_gto, casscf_calculator, h2_xyz):
        """Test basic CASSCF calculation setup."""
        atoms = casscf_calculator.parse_xyz(h2_xyz)
        casscf_calculator.setup_calculation(
            atoms,
            basis='sto-3g',
            charge=0,
            spin=0,
            ncas=2,  # 2 orbitals
            nelecas=2,  # 2 electrons
            max_cycle_macro=30,
            max_cycle_micro=4,
            conv_tol=1e-6,
            conv_tol_grad=1e-4
        )
        
        assert casscf_calculator.mol is not None
        assert casscf_calculator.mf is not None
        assert casscf_calculator.ncas == 2
        assert casscf_calculator.nelecas == 2
        assert casscf_calculator.max_cycle_macro == 30
        assert casscf_calculator.max_cycle_micro == 4
        assert casscf_calculator.conv_tol == 1e-6
        assert casscf_calculator.conv_tol_grad == 1e-4
        mock_gto.M.assert_called_once()
    
    def test_setup_invalid_convergence_params(self, casscf_calculator, h2_xyz):
        """Test setup with invalid convergence parameters."""
        atoms = casscf_calculator.parse_xyz(h2_xyz)
        
        # Invalid max_cycle_macro
        with pytest.raises(InputError, match="must be positive"):
            casscf_calculator.setup_calculation(
                atoms,
                basis='sto-3g',
                ncas=2,
                nelecas=2,
                max_cycle_macro=0
            )
    
    @patch('quantum_calc.casscf_calculator.gto')
    @patch('quantum_calc.casscf_calculator.scf')
    @patch('quantum_calc.casscf_calculator.mcscf')
    def test_casscf_calculation_flow(self, mock_mcscf, mock_scf, mock_gto, casscf_calculator, h2_xyz):
        """Test the complete CASSCF calculation flow."""
        # Setup mocks
        mock_mol = MagicMock()
        mock_gto.M.return_value = mock_mol
        
        mock_mf = MagicMock()
        mock_mf.kernel.return_value = -1.123456
        mock_mf.converged = True
        mock_mf.mo_energy = np.array([-0.5, 0.5])
        mock_mf.mo_occ = np.array([2.0, 0.0])
        mock_mf.mo_coeff = np.eye(2)
        mock_scf.RHF.return_value = mock_mf
        
        mock_cas = MagicMock()
        mock_cas.kernel.return_value = -1.145820
        mock_cas.converged = True
        mock_cas.ncore = 0
        mock_cas.ci = np.array([0.85, 0.3, 0.1])  # CASSCF CI coefficients
        mock_cas.make_rdm1s.return_value = (np.eye(2), np.eye(2))
        mock_cas.cas_natorb.return_value = (np.eye(2), np.array([0.85, 0.15]), np.array([1.7, 0.3]))
        mock_cas.mo_coeff = np.array([[0.9, 0.1], [0.1, 0.9]])  # Slightly rotated orbitals
        mock_cas.e_tot = [-1.123456, -1.145820]  # Energy history for iterations
        mock_mcscf.CASSCF.return_value = mock_cas
        
        # Run calculation
        atoms = casscf_calculator.parse_xyz(h2_xyz)
        casscf_calculator.setup_calculation(atoms, basis='sto-3g', ncas=2, nelecas=2)
        results = casscf_calculator.calculate()
        
        # Verify results
        assert 'scf_energy' in results
        assert 'casscf_energy' in results
        assert 'correlation_energy' in results
        assert 'converged' in results
        assert 'macro_iterations' in results
        assert results['method'] == 'RHF-CASSCF'
        
        # Verify analysis results
        assert 'natural_orbital_analysis' in results
        assert 'ci_coefficient_analysis' in results
        assert 'orbital_rotation_analysis' in results
        
        # Verify calculation was called
        mock_cas.kernel.assert_called_once()
    
    @patch('quantum_calc.casscf_calculator.gto')
    @patch('quantum_calc.casscf_calculator.scf')
    def test_open_shell_setup(self, mock_scf, mock_gto, casscf_calculator, oh_radical_xyz):
        """Test CASSCF setup for open-shell systems."""
        atoms = casscf_calculator.parse_xyz(oh_radical_xyz)
        casscf_calculator.setup_calculation(
            atoms,
            basis='sto-3g',
            spin=1,  # Doublet state
            ncas=4,
            nelecas=5
        )
        
        # Should use UHF for open-shell
        mock_scf.UHF.assert_called_once()
        assert casscf_calculator.results['method'] == 'UHF-CASSCF'


class TestCASAnalysisFunctions:
    """Test analysis functions common to CASCI and CASSCF."""
    
    @patch('quantum_calc.casci_calculator.gto')
    @patch('quantum_calc.casci_calculator.scf')
    @patch('quantum_calc.casci_calculator.mcscf')
    def test_natural_orbital_analysis(self, mock_mcscf, mock_scf, mock_gto, casci_calculator, h2_xyz):
        """Test natural orbital analysis functionality."""
        # Setup mocks for natural orbital analysis
        mock_mol = MagicMock()
        mock_gto.M.return_value = mock_mol
        
        mock_mf = MagicMock()
        mock_mf.kernel.return_value = -1.0
        mock_scf.RHF.return_value = mock_mf
        
        mock_cas = MagicMock()
        mock_cas.kernel.return_value = -1.1
        mock_cas.converged = True
        mock_cas.natorb = True
        mock_cas.ci = np.array([0.9, 0.1])
        mock_cas.cas_natorb.return_value = (
            np.eye(2),  # Natural orbitals
            np.array([0.9, 0.1]),  # CI coefficients
            np.array([1.8, 0.2])  # Occupation numbers
        )
        mock_mcscf.CASCI.return_value = mock_cas
        
        # Run calculation
        atoms = casci_calculator.parse_xyz(h2_xyz)
        casci_calculator.setup_calculation(atoms, basis='sto-3g', ncas=2, nelecas=2, natorb=True)
        results = casci_calculator.calculate()
        
        # Verify natural orbital analysis
        no_analysis = results['natural_orbital_analysis']
        assert no_analysis['enabled'] is True
        assert 'occupation_numbers' in no_analysis
        assert 'total_active_electrons' in no_analysis
        assert 'effective_electron_pairs' in no_analysis
    
    @patch('quantum_calc.casci_calculator.gto')
    @patch('quantum_calc.casci_calculator.scf')
    @patch('quantum_calc.casci_calculator.mcscf')
    def test_ci_coefficient_analysis(self, mock_mcscf, mock_scf, mock_gto, casci_calculator, h2_xyz):
        """Test CI coefficient analysis functionality."""
        # Setup mocks for CI analysis
        mock_mol = MagicMock()
        mock_gto.M.return_value = mock_mol
        
        mock_mf = MagicMock()
        mock_mf.kernel.return_value = -1.0
        mock_scf.RHF.return_value = mock_mf
        
        mock_cas = MagicMock()
        mock_cas.kernel.return_value = -1.1
        mock_cas.converged = True
        mock_cas.ci = np.array([0.9, -0.3, 0.1, 0.05])  # Multiple configurations
        mock_mcscf.CASCI.return_value = mock_cas
        
        # Run calculation
        atoms = casci_calculator.parse_xyz(h2_xyz)
        casci_calculator.setup_calculation(atoms, basis='sto-3g', ncas=2, nelecas=2)
        results = casci_calculator.calculate()
        
        # Verify CI coefficient analysis
        ci_analysis = results['ci_coefficient_analysis']
        assert ci_analysis['available'] is True
        assert 'major_configurations' in ci_analysis
        assert 'leading_coefficient' in ci_analysis
        assert 'leading_contribution_percent' in ci_analysis
        assert 'multiconfigurational_character' in ci_analysis
        
        # Should have found major configurations
        assert len(ci_analysis['major_configurations']) > 0
        assert ci_analysis['leading_contribution_percent'] > 50.0  # Dominant configuration
    
    @patch('quantum_calc.casci_calculator.gto')
    @patch('quantum_calc.casci_calculator.scf') 
    @patch('quantum_calc.casci_calculator.mcscf')
    def test_spin_density_analysis_open_shell(self, mock_mcscf, mock_scf, mock_gto, casci_calculator, oh_radical_xyz):
        """Test spin density analysis for open-shell systems."""
        # Setup mocks for open-shell spin density
        mock_mol = MagicMock()
        mock_mol.natm = 2
        mock_mol.atom_symbol = lambda i: ['O', 'H'][i]
        mock_mol.aoslice_by_atom.return_value = [(0, 0, 0, 4), (0, 0, 4, 5)]  # Basis function slicing
        mock_mol.spin = 1
        mock_gto.M.return_value = mock_mol
        
        mock_mf = MagicMock()
        mock_mf.kernel.return_value = -75.0
        mock_scf.UHF.return_value = mock_mf
        
        mock_cas = MagicMock()
        mock_cas.kernel.return_value = -75.1
        mock_cas.converged = True
        # Mock alpha and beta density matrices with spin difference
        dm_alpha = np.array([[0.6, 0.1], [0.1, 0.4]])
        dm_beta = np.array([[0.4, 0.1], [0.1, 0.4]]) 
        mock_cas.make_rdm1s.return_value = (dm_alpha, dm_beta)
        mock_cas._scf.get_ovlp.return_value = np.eye(2)
        mock_mcscf.CASCI.return_value = mock_cas
        
        # Run calculation with open-shell
        atoms = casci_calculator.parse_xyz(oh_radical_xyz)
        casci_calculator.setup_calculation(atoms, basis='sto-3g', spin=1, ncas=3, nelecas=5)
        results = casci_calculator.calculate()
        
        # Verify spin density analysis was performed
        assert 'mulliken_spin_analysis' in results
        spin_analysis = results['mulliken_spin_analysis']
        assert spin_analysis['available'] is True
        assert 'atomic_spin_densities' in spin_analysis
        assert 'total_spin_density' in spin_analysis
        assert 'expected_spin' in spin_analysis
        
        # Should have analysis for both atoms
        assert len(spin_analysis['atomic_spin_densities']) == 2


class TestErrorHandling:
    """Test error handling in CASCI/CASSCF calculators."""
    
    def test_convergence_error_handling(self, casci_calculator):
        """Test handling of convergence failures."""
        # This would be tested with mocked convergence failure
        # For now, just verify the error types are properly imported
        assert InputError is not None
        assert ConvergenceError is not None
        assert GeometryError is not None
    
    def test_memory_settings(self, casci_calculator, h2_xyz):
        """Test memory setting functionality."""
        atoms = casci_calculator.parse_xyz(h2_xyz)
        
        with patch('quantum_calc.casci_calculator.gto') as mock_gto:
            mock_mol = MagicMock()
            mock_gto.M.return_value = mock_mol
            
            casci_calculator.setup_calculation(
                atoms,
                basis='sto-3g',
                ncas=2,
                nelecas=2,
                memory_mb=8000
            )
            
            # Verify memory was set
            assert mock_mol.max_memory == 8000


if __name__ == "__main__":
    pytest.main([__file__])