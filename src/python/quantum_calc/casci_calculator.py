"""CASCI calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from pyscf import gto, scf, mcscf

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects

logger = logging.getLogger(__name__)


class CASCICalculator(BaseCalculator):
    """CASCI calculator using PySCF for multi-configurational calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None, optimize_geometry: bool = False):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir, optimize_geometry)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None  # Can be RHF or UHF
        self.mycas: Optional[mcscf.casci.CASCI] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        # Set calculation method for memory management
        self.calculation_method = 'CASCI'
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup CASCI calculation using the base template method."""
        # Call the base template method which handles common setup
        super().setup_calculation(atoms, **kwargs)
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """Validate CASCI-specific parameters."""
        # CASCI-specific parameters
        ncas = kwargs.get('ncas', 4)  # Number of active space orbitals
        nelecas = kwargs.get('nelecas', 4)  # Number of active space electrons
        analyze_nto = kwargs.get('analyze_nto', False)  # Natural transition orbitals
        natorb = kwargs.get('natorb', True)  # Transform to natural orbitals
        
        # Validate and adjust CASCI parameters using base class method
        ncas, nelecas = self._validate_active_space_parameters(ncas, nelecas)
        
        # Store parameters for template method access
        self.ncas = ncas
        self.nelecas = nelecas
        self.analyze_nto = analyze_nto
        self.natorb = natorb
        
        return {
            'ncas': ncas,
            'nelecas': nelecas,
            'analyze_nto': analyze_nto,
            'natorb': natorb,
            'method': 'UHF-CASCI' if kwargs.get('spin', 0) > 0 else 'RHF-CASCI'
        }
    
    def _get_default_memory_mb(self) -> int:
        """Get default memory setting for CASCI calculations."""
        return 6000  # CASCI needs substantial memory (6GB default)
    
    def _get_calculation_method_name(self) -> str:
        """Get the name of the calculation method for logging."""
        return 'CASCI'
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform CASCI-specific calculation after SCF."""
        # Safely convert base SCF energy
        if base_energy is None:
            raise CalculationError("SCF energy is None - calculation may have failed")
        try:
            scf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert SCF energy to float: {e}")
        
        logger.info(f"Reference SCF completed. SCF energy: {scf_energy_float:.6f} hartree")
        
        # CASCI calculation
        logger.info("Starting CASCI calculation...")
        ncas = getattr(self, 'ncas', 4)
        nelecas = getattr(self, 'nelecas', 4)
        natorb = getattr(self, 'natorb', True)
        
        # Create CASCI object
        self.mycas = mcscf.CASCI(self.mf, ncas, nelecas)
        self.mycas.natorb = natorb
        
        # Validate active space parameters against molecular orbital structure
        if hasattr(self.mf, 'mo_energy') and self.mf.mo_energy is not None:
            n_orb = len(self.mf.mo_energy)
            n_occupied = int(np.sum(self.mf.mo_occ > 0)) if hasattr(self.mf, 'mo_occ') else 0
            
            if ncas > n_orb:
                raise InputError(f"Active space size ({ncas}) exceeds total orbitals ({n_orb})")
            if ncas > n_occupied + 10:  # Reasonable check for too many virtual orbitals
                logger.warning(f"Large active space ({ncas} orbitals) includes many virtual orbitals. "
                             f"System has {n_occupied} occupied orbitals out of {n_orb} total.")
        
        logger.info(f"CASCI setup: {ncas} active orbitals, {nelecas} active electrons")
        if natorb:
            logger.info("Natural orbital transformation enabled")
        
        # Run CASCI calculation
        logger.info("Starting CASCI kernel calculation...")
        try:
            casci_result = self.mycas.kernel()
        except Exception as e:
            # Use base class error handling
            self._handle_calculation_error(e, 'CASCI')
        
        # Use base class method to analyze kernel return value
        casci_energy_float, ci_coefficients, kernel_additional_info = self._analyze_kernel_return_value(casci_result, 'CASCI')
        
        # Store CI coefficients for analysis
        if ci_coefficients is not None:
            self.ci_coefficients = ci_coefficients
        elif hasattr(self.mycas, 'ci') and self.mycas.ci is not None:
            self.ci_coefficients = self.mycas.ci
            logger.info("Using CI coefficients from CASCI object")
        else:
            self.ci_coefficients = None
        
        if not self.mycas.converged:
            logger.warning("CASCI calculation did not converge to specified tolerance")
        
        # Analyze CASCI results using base class methods
        casci_results = self._analyze_casci_results()
        
        # Enhanced CI coefficients analysis if available
        enhanced_ci_analysis = {}
        if self.ci_coefficients is not None:
            enhanced_ci_analysis = self._analyze_enhanced_ci_coefficients(self.ci_coefficients)
        
        # Combine results with comprehensive information
        results = {
            'scf_energy': scf_energy_float,
            'casci_energy': casci_energy_float,
            'correlation_energy': casci_energy_float - scf_energy_float,
            'kernel_return_info': kernel_additional_info,
            'ci_coefficients_available': self.ci_coefficients is not None,
            'enhanced_ci_analysis': enhanced_ci_analysis
        }
        results.update(casci_results)
        
        logger.info(f"CASCI calculation completed with enhanced analysis. Energy: {casci_energy_float:.6f} hartree")
        return results
    
    def _create_scf_method(self, mol):
        """Create SCF method object for CASCI reference (RHF/UHF)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = scf.RHF(mol)
            logger.info("Using Restricted Hartree-Fock (RHF) for closed-shell CASCI reference")
        else:
            mf = scf.UHF(mol)
            logger.info("Using Unrestricted Hartree-Fock (UHF) for open-shell CASCI reference")
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to SCF method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UHF' if spin > 0 else 'RHF'} (CASCI reference)"
    
    def _requires_geometry_optimization(self) -> bool:
        """CASCI does not require geometry optimization by default."""
        return False
    
    def _requires_frequency_analysis(self) -> bool:
        """CASCI does not require frequency analysis by default."""
        return False
    
    # ===== CASCI-specific Analysis Methods =====
    
    def _analyze_casci_results(self) -> Dict[str, Any]:
        """
        Perform comprehensive analysis of CASCI results using base class methods.
        """
        if self.mycas is None:
            raise CalculationError("CASCI calculation not available for analysis")
        
        logger.info("Starting CASCI results analysis...")
        
        # Set verbose level for detailed analysis
        original_verbose = self.mycas.verbose
        self.mycas.verbose = 4
        
        # Perform PySCF's built-in analysis
        logger.info("Running PySCF CASCI analysis...")
        try:
            self.mycas.analyze()
        except Exception as e:
            logger.warning(f"PySCF CASCI analysis encountered issues: {e}")
        
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
        
        logger.info("CASCI results analysis completed")
        return analysis_results



