"""CASSCF calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from pyscf import gto, scf, mcscf

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects
from .config_manager import get_memory_for_method, get_max_cycle_macro, get_max_cycle_micro, get_ah_max_cycle

logger = logging.getLogger(__name__)


class CASSCFCalculator(BaseCalculator):
    """CASSCF calculator using PySCF for multi-configurational self-consistent field calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None, optimize_geometry: bool = False):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir, optimize_geometry)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None  # Can be RHF or UHF
        self.mycas: Optional[mcscf.casscf.CASSCF] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        # Set calculation method for memory management
        self.calculation_method = 'CASSCF'
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup CASSCF calculation using the base template method."""
        # Call the base template method which handles common setup
        super().setup_calculation(atoms, **kwargs)
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """Validate CASSCF-specific parameters."""
        # CASSCF-specific parameters
        ncas = kwargs.get('ncas', 4)  # Number of active space orbitals
        nelecas = kwargs.get('nelecas', 4)  # Number of active space electrons
        max_cycle_macro = kwargs.get('max_cycle_macro', get_max_cycle_macro())  # CASSCF macro iterations
        max_cycle_micro = kwargs.get('max_cycle_micro', get_max_cycle_micro())  # CI solver micro iterations
        analyze_nto = kwargs.get('analyze_nto', False)  # Natural transition orbitals
        natorb = kwargs.get('natorb', True)  # Transform to natural orbitals
        conv_tol = kwargs.get('conv_tol', 1e-6)  # Convergence tolerance
        conv_tol_grad = kwargs.get('conv_tol_grad', 1e-4)  # Gradient tolerance
        
        # AH solver parameters for improved convergence
        ah_conv_tol = kwargs.get('ah_conv_tol', 1e-12)  # AH solver convergence tolerance
        ah_max_cycle = kwargs.get('ah_max_cycle', get_ah_max_cycle())  # AH solver max iterations
        ah_lindep = kwargs.get('ah_lindep', 1e-14)  # AH linear dependence threshold
        
        # Validate and adjust CASSCF parameters using base class method
        ncas, nelecas = self._validate_active_space_parameters(ncas, nelecas)
        
        # Validate max_cycle_macro
        if max_cycle_macro <= 0:
            logger.warning(f"Invalid max_cycle_macro={max_cycle_macro}, adjusting to default 50")
            max_cycle_macro = 50
        if max_cycle_macro > 100:
            logger.warning(f"Many CASSCF macro cycles ({max_cycle_macro}) may require long computation time")
        
        # Store parameters for template method access
        self.ncas = ncas
        self.nelecas = nelecas
        self.max_cycle_macro = max_cycle_macro
        self.max_cycle_micro = max_cycle_micro
        self.analyze_nto = analyze_nto
        self.natorb = natorb
        self.conv_tol = conv_tol
        self.conv_tol_grad = conv_tol_grad
        self.ah_conv_tol = ah_conv_tol
        self.ah_max_cycle = ah_max_cycle
        self.ah_lindep = ah_lindep
        
        return {
            'ncas': ncas,
            'nelecas': nelecas,
            'max_cycle_macro': max_cycle_macro,
            'max_cycle_micro': max_cycle_micro,
            'analyze_nto': analyze_nto,
            'natorb': natorb,
            'conv_tol': conv_tol,
            'conv_tol_grad': conv_tol_grad,
            'ah_conv_tol': ah_conv_tol,
            'ah_max_cycle': ah_max_cycle,
            'ah_lindep': ah_lindep,
            'method': 'UHF-CASSCF' if kwargs.get('spin', 0) > 0 else 'RHF-CASSCF'
        }
    
    def _get_default_memory_mb(self) -> int:
        """Get default memory setting for CASSCF calculations."""
        return get_memory_for_method('CASSCF')
    
    def _get_calculation_method_name(self) -> str:
        """Get the name of the calculation method for logging."""
        return 'CASSCF'
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform CASSCF-specific calculation after SCF."""
        # Safely convert base SCF energy
        if base_energy is None:
            raise CalculationError("SCF energy is None - calculation may have failed")
        try:
            scf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert SCF energy to float: {e}")
        
        logger.info(f"Reference SCF completed. SCF energy: {scf_energy_float:.6f} hartree")
        
        # CASSCF calculation
        logger.info("Starting CASSCF calculation...")
        ncas = getattr(self, 'ncas', 4)
        nelecas = getattr(self, 'nelecas', 4)
        natorb = getattr(self, 'natorb', True)
        max_cycle_macro = getattr(self, 'max_cycle_macro', 50)
        max_cycle_micro = getattr(self, 'max_cycle_micro', 3)
        conv_tol = getattr(self, 'conv_tol', 1e-6)
        conv_tol_grad = getattr(self, 'conv_tol_grad', 1e-4)
        ah_conv_tol = getattr(self, 'ah_conv_tol', 1e-12)
        ah_max_cycle = getattr(self, 'ah_max_cycle', 30)
        ah_lindep = getattr(self, 'ah_lindep', 1e-14)
        
        # Create CASSCF object
        self.mycas = mcscf.CASSCF(self.mf, ncas, nelecas)
        
        # Set CASSCF parameters
        self.mycas.natorb = natorb
        self.mycas.max_cycle_macro = max_cycle_macro
        self.mycas.max_cycle_micro = max_cycle_micro
        self.mycas.conv_tol = conv_tol
        self.mycas.conv_tol_grad = conv_tol_grad
        
        # Set AH solver parameters for improved convergence
        self.mycas.ah_conv_tol = ah_conv_tol
        self.mycas.ah_max_cycle = ah_max_cycle
        self.mycas.ah_lindep = ah_lindep
        
        # Validate active space parameters against molecular orbital structure
        if hasattr(self.mf, 'mo_energy') and self.mf.mo_energy is not None:
            n_orb = len(self.mf.mo_energy)
            n_occupied = int(np.sum(self.mf.mo_occ > 0)) if hasattr(self.mf, 'mo_occ') else 0
            
            if ncas > n_orb:
                raise InputError(f"Active space size ({ncas}) exceeds total orbitals ({n_orb})")
            if ncas > n_occupied + 10:  # Reasonable check for too many virtual orbitals
                logger.warning(f"Large active space ({ncas} orbitals) includes many virtual orbitals. "
                             f"System has {n_occupied} occupied orbitals out of {n_orb} total.")
        
        logger.info(f"CASSCF setup: {ncas} active orbitals, {nelecas} active electrons")
        logger.info(f"Convergence criteria: energy = {conv_tol:.2e}, gradient = {conv_tol_grad:.2e}")
        logger.info(f"Maximum cycles: macro = {max_cycle_macro}, micro = {max_cycle_micro}")
        if natorb:
            logger.info("Natural orbital transformation enabled")
        
        # Run CASSCF calculation
        logger.info("Starting CASSCF kernel calculation...")
        try:
            casscf_result = self.mycas.kernel()
        except Exception as e:
            # Use base class error handling
            self._handle_calculation_error(e, 'CASSCF')
        
        # Use base class method to analyze kernel return value
        casscf_energy_float, ci_coefficients, kernel_additional_info = self._analyze_kernel_return_value(casscf_result, 'CASSCF')
        
        # Store CI coefficients for analysis
        if ci_coefficients is not None:
            self.ci_coefficients = ci_coefficients
        elif hasattr(self.mycas, 'ci') and self.mycas.ci is not None:
            self.ci_coefficients = self.mycas.ci
            logger.info("Using CI coefficients from CASSCF object")
        else:
            self.ci_coefficients = None
        
        # Check convergence status
        converged = getattr(self.mycas, 'converged', False)
        if not converged:
            logger.warning("CASSCF calculation did not converge to specified tolerance")
        else:
            logger.info("CASSCF calculation converged successfully")
        
        # Get iteration information
        e_tot_history = getattr(self.mycas, 'e_tot', None)
        if hasattr(e_tot_history, '__len__') and len(e_tot_history) > 1:
            iterations_performed = len(e_tot_history) - 1
        else:
            iterations_performed = getattr(self.mycas, 'itercount', 0)
        
        logger.info(f"CASSCF calculation completed. Energy: {casscf_energy_float:.6f} hartree")
        logger.info(f"Performed {iterations_performed} macro iterations")
        
        # Analyze CASSCF results using base class methods
        casscf_results = self._analyze_casscf_results()
        
        # Enhanced CI coefficients analysis if available
        enhanced_ci_analysis = {}
        if self.ci_coefficients is not None:
            enhanced_ci_analysis = self._analyze_enhanced_ci_coefficients(self.ci_coefficients)
        
        # Combine results with comprehensive information
        results = {
            'scf_energy': scf_energy_float,
            'casscf_energy': casscf_energy_float,
            'correlation_energy': casscf_energy_float - scf_energy_float,
            'converged': converged,
            'macro_iterations': iterations_performed,
            'kernel_return_info': kernel_additional_info,
            'ci_coefficients_available': self.ci_coefficients is not None,
            'enhanced_ci_analysis': enhanced_ci_analysis
        }
        results.update(casscf_results)
        
        logger.info(f"CASSCF calculation completed with enhanced analysis. Energy: {casscf_energy_float:.6f} hartree, Converged: {converged}")
        return results
    
    def _create_scf_method(self, mol):
        """Create SCF method object for CASSCF reference (RHF/UHF)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = scf.RHF(mol)
            logger.info("Using Restricted Hartree-Fock (RHF) for closed-shell CASSCF reference")
        else:
            mf = scf.UHF(mol)
            logger.info("Using Unrestricted Hartree-Fock (UHF) for open-shell CASSCF reference")
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to SCF method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UHF' if spin > 0 else 'RHF'} (CASSCF reference)"
    
    def _requires_geometry_optimization(self) -> bool:
        """CASSCF does not require geometry optimization by default."""
        return False
    
    def _requires_frequency_analysis(self) -> bool:
        """CASSCF does not require frequency analysis by default."""
        return False
    
    # ===== CASSCF-specific Analysis Methods =====
    
    def _analyze_casscf_results(self) -> Dict[str, Any]:
        """
        Perform comprehensive analysis of CASSCF results using base class methods.
        """
        if self.mycas is None:
            raise CalculationError("CASSCF calculation not available for analysis")
        
        logger.info("Starting CASSCF results analysis...")
        
        # Set verbose level for detailed analysis
        original_verbose = self.mycas.verbose
        self.mycas.verbose = 4
        
        # Perform PySCF's built-in analysis
        logger.info("Running PySCF CASSCF analysis...")
        try:
            self.mycas.analyze()
        except Exception as e:
            logger.warning(f"PySCF CASSCF analysis encountered issues: {e}")
        
        # Restore original verbose level
        self.mycas.verbose = original_verbose
        
        # Extract analysis results using base class methods
        analysis_results = {}
        
        # Natural orbital analysis
        try:
            natural_orbital_analysis = self._analyze_natural_orbitals()
            analysis_results['natural_orbital_analysis'] = natural_orbital_analysis
            logger.info("Natural orbital analysis completed successfully")
        except Exception as e:
            logger.error(f"Natural orbital analysis failed: {e}")
            analysis_results['natural_orbital_analysis'] = None
        
        # CI coefficient analysis
        try:
            ci_analysis = self._analyze_ci_coefficients()
            analysis_results['ci_coefficient_analysis'] = ci_analysis
            logger.info("CI coefficient analysis completed successfully")
        except Exception as e:
            logger.error(f"CI coefficient analysis failed: {e}")
            analysis_results['ci_coefficient_analysis'] = None
        
        # Orbital rotation analysis (CASSCF specific)
        try:
            rotation_analysis = self._analyze_orbital_rotation()
            analysis_results['orbital_rotation_analysis'] = rotation_analysis
            logger.info("Orbital rotation analysis completed successfully")
        except Exception as e:
            logger.error(f"Orbital rotation analysis failed: {e}")
            analysis_results['orbital_rotation_analysis'] = None
        
        # Mulliken spin density analysis (for open-shell systems)
        if self.results.get('spin', 0) > 0:
            try:
                spin_analysis = self._calculate_mulliken_spin_density()
                analysis_results['mulliken_spin_analysis'] = spin_analysis
                logger.info("Mulliken spin density analysis completed successfully")
            except Exception as e:
                logger.error(f"Mulliken spin density analysis failed: {e}")
                analysis_results['mulliken_spin_analysis'] = None
        
        # Orbital overlap analysis
        try:
            overlap_analysis = self._calculate_orbital_overlap()
            analysis_results['orbital_overlap_analysis'] = overlap_analysis
            logger.info("Orbital overlap analysis completed successfully")
        except Exception as e:
            logger.error(f"Orbital overlap analysis failed: {e}")
            analysis_results['orbital_overlap_analysis'] = None
        
        logger.info("CASSCF results analysis completed")
        return analysis_results
    

    
    def _analyze_orbital_rotation(self) -> Dict[str, Any]:
        """
        Analyze the orbital rotation that occurred during CASSCF optimization.
        This is specific to CASSCF and not available in CASCI.
        """
        rotation_analysis = {'available': True}
        
        try:
            # Check if we have rotation information
            if hasattr(self.mycas, 'mo_coeff') and hasattr(self.mf, 'mo_coeff'):
                # Calculate orbital rotation matrix
                if hasattr(self.mol, 'get_ovlp'):
                    S = self.mol.get_ovlp()
                else:
                    S = self.mf.get_ovlp()
                
                # Rotation matrix: R = C_initial^T * S * C_final
                rotation_matrix = np.dot(self.mf.mo_coeff.T, np.dot(S, self.mycas.mo_coeff))
                
                # Analyze rotation for different orbital spaces
                ncore = self.mycas.ncore
                ncas = self.results.get('ncas', 4)
                
                # Core orbital rotations
                if ncore > 0:
                    core_rotation = rotation_matrix[:ncore, :ncore]
                    core_rotation_magnitude = np.max(np.abs(core_rotation - np.eye(ncore)))
                else:
                    core_rotation_magnitude = 0.0
                
                # Active orbital rotations
                active_start = ncore
                active_end = ncore + ncas
                active_rotation = rotation_matrix[active_start:active_end, active_start:active_end]
                active_rotation_magnitude = np.max(np.abs(active_rotation - np.eye(ncas)))
                
                # Virtual orbital rotations (if any)
                nvirt = rotation_matrix.shape[0] - active_end
                if nvirt > 0:
                    virtual_rotation = rotation_matrix[active_end:, active_end:]
                    virtual_rotation_magnitude = np.max(np.abs(virtual_rotation - np.eye(nvirt)))
                else:
                    virtual_rotation_magnitude = 0.0
                
                rotation_analysis.update({
                    'core_orbital_rotation_magnitude': float(core_rotation_magnitude),
                    'active_orbital_rotation_magnitude': float(active_rotation_magnitude),
                    'virtual_orbital_rotation_magnitude': float(virtual_rotation_magnitude),
                    'overall_rotation_magnitude': float(np.max([
                        core_rotation_magnitude, 
                        active_rotation_magnitude, 
                        virtual_rotation_magnitude
                    ]))
                })
                
                # Classify rotation extent
                max_rotation = rotation_analysis['overall_rotation_magnitude']
                if max_rotation < 0.1:
                    rotation_extent = "minimal"
                elif max_rotation < 0.3:
                    rotation_extent = "moderate"
                else:
                    rotation_extent = "significant"
                
                rotation_analysis['rotation_extent'] = rotation_extent
                
                logger.info(f"Orbital rotation analysis: core = {core_rotation_magnitude:.3f}, "
                           f"active = {active_rotation_magnitude:.3f}, "
                           f"virtual = {virtual_rotation_magnitude:.3f} ({rotation_extent})")
            else:
                rotation_analysis.update({
                    'available': False,
                    'reason': 'Initial or final orbital coefficients not available'
                })
                
        except Exception as e:
            logger.error(f"Error in orbital rotation analysis: {e}")
            rotation_analysis['error'] = str(e)
        
        return rotation_analysis
