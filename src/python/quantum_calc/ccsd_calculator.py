"""CCSD calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, scf, cc
# geometric_solver is now imported in BaseCalculator

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects
from .config_manager import get_memory_for_method

logger = logging.getLogger(__name__)


class CCSDCalculator(BaseCalculator):
    """CCSD/CCSD(T) calculator using PySCF for structure optimization and CCSD calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None, optimize_geometry: bool = False):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir, optimize_geometry)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None
        self.ccsd: Optional[cc.CCSD] = None
        self.optimized_geometry: Optional[np.ndarray] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        self.calculate_ccsd_t = False  # Flag to indicate if CCSD(T) should be calculated
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup CCSD calculation using the base template method."""
        # Call the base template method which handles common setup
        super().setup_calculation(atoms, **kwargs)
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """Validate CCSD-specific parameters."""
        # CCSD-specific parameters
        frozen_core = kwargs.get('frozen_core', True)  # Enable frozen core by default for CCSD
        ccsd_t = kwargs.get('ccsd_t', False)  # Whether to calculate CCSD(T) correction
        
        # Store CCSD(T) flag for later use
        self.calculate_ccsd_t = ccsd_t
        self.frozen_core = frozen_core
        
        # Set method name based on spin and CCSD(T) flag
        method_name = 'UCCSD' if kwargs.get('spin', 0) > 0 else 'RCCSD'
        if ccsd_t:
            method_name += '(T)'
        
        return {
            'method': method_name,
            'frozen_core': frozen_core,
            'ccsd_t_correction': ccsd_t
        }
    
    def _get_default_memory_mb(self) -> int:
        """Get default memory setting for CCSD calculations from config."""
        return get_memory_for_method('CCSD')
    
    def _get_calculation_method_name(self) -> str:
        """Get the name of the calculation method for logging."""
        return 'CCSD(T)' if self.calculate_ccsd_t else 'CCSD'
    
    # ===== Template Method Pattern Implementation =====

    def _extract_ccsd_diagnostics(self) -> Dict[str, Any]:
        """
        Extract CCSD diagnostic indicators.

        T1 diagnostic: Indicates the reliability of single-reference methods
                      T1 > 0.02 suggests multi-reference character
        D1 diagnostic: Alternative diagnostic based on density matrix
        D2 diagnostic: Diagnostic based on T2 amplitudes
        """
        diagnostics = {}

        if not hasattr(self, 'ccsd') or self.ccsd is None:
            logger.warning("CCSD object not available for diagnostics extraction")
            return diagnostics

        try:
            # T1 diagnostic - most important indicator
            if hasattr(self.ccsd, 'get_t1_diagnostic') and callable(self.ccsd.get_t1_diagnostic):
                t1_diag = self.ccsd.get_t1_diagnostic()
                diagnostics['ccsd_t1_diagnostic'] = float(t1_diag)
                logger.info(f"CCSD T1 diagnostic: {diagnostics['ccsd_t1_diagnostic']:.6f}")

                # Interpret T1 diagnostic
                if t1_diag > 0.02:
                    logger.warning(f"T1 diagnostic ({t1_diag:.6f}) > 0.02: Multi-reference character detected - single-reference CCSD may be unreliable")
                else:
                    logger.info(f"T1 diagnostic ({t1_diag:.6f}) <= 0.02: Single-reference CCSD is appropriate")

            # D1 diagnostic
            if hasattr(self.ccsd, 'get_d1_diagnostic') and callable(self.ccsd.get_d1_diagnostic):
                d1_diag = self.ccsd.get_d1_diagnostic()
                diagnostics['ccsd_d1_diagnostic'] = float(d1_diag)
                logger.info(f"CCSD D1 diagnostic: {diagnostics['ccsd_d1_diagnostic']:.6f}")

            # D2 diagnostic
            if hasattr(self.ccsd, 'get_d2_diagnostic') and callable(self.ccsd.get_d2_diagnostic):
                d2_diag = self.ccsd.get_d2_diagnostic()
                diagnostics['ccsd_d2_diagnostic'] = float(d2_diag)
                logger.info(f"CCSD D2 diagnostic: {diagnostics['ccsd_d2_diagnostic']:.6f}")

            logger.info("CCSD diagnostics extraction completed")

        except Exception as e:
            logger.warning(f"Failed to extract CCSD diagnostics: {e}")

        return diagnostics

    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform CCSD-specific calculation after HF."""
        # 安全にベースHFエネルギーを変換
        if base_energy is None:
            raise CalculationError("HF energy is None - calculation may have failed")
        try:
            hf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert HF energy to float: {e}")

        logger.info(f"HF energy: {hf_energy_float} Hartree")

        # CCSD calculation
        logger.info("Starting CCSD calculation...")

        # Create CCSD object based on HF reference
        self.ccsd = cc.CCSD(self.mf)

        # Apply frozen core if requested
        if getattr(self, 'frozen_core', True):
            # Automatically set frozen core orbitals
            self.ccsd.set_frozen()
            logger.info("Using frozen core approximation")

        # Run CCSD calculation
        ccsd_energy, t1, t2 = self.ccsd.kernel()

        # Check CCSD convergence
        if not self.ccsd.converged:
            raise ConvergenceError("CCSD calculation failed to converge")

        # Get CCSD results
        ccsd_corr_energy = self.ccsd.e_corr
        ccsd_total_energy = self.ccsd.e_tot

        logger.info("CCSD calculation completed successfully.")
        logger.info(f"CCSD correlation energy: {ccsd_corr_energy} Hartree")
        logger.info(f"CCSD total energy: {ccsd_total_energy} Hartree")

        # CCSD(T) calculation if requested
        ccsd_t_correction = 0.0
        ccsd_t_total_energy = ccsd_total_energy

        if self.calculate_ccsd_t:
            logger.info("Starting CCSD(T) perturbative triples correction...")
            ccsd_t_correction = self.ccsd.ccsd_t()
            ccsd_t_total_energy = ccsd_total_energy + ccsd_t_correction

            logger.info("CCSD(T) calculation completed successfully.")
            logger.info(f"CCSD(T) triples correction: {ccsd_t_correction} Hartree")
            logger.info(f"CCSD(T) total energy: {ccsd_t_total_energy} Hartree")

        # 安全にCCSDエネルギーを変換
        if ccsd_corr_energy is None:
            raise CalculationError("CCSD correlation energy is None - calculation may have failed")
        if ccsd_total_energy is None:
            raise CalculationError("CCSD total energy is None - calculation may have failed")

        try:
            ccsd_corr_energy_float = float(ccsd_corr_energy)
            ccsd_total_energy_float = float(ccsd_total_energy)
            ccsd_t_correction_float = float(ccsd_t_correction)
            ccsd_t_total_energy_float = float(ccsd_t_total_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert CCSD energies to float: {e}")

        results = {
            'hf_energy': hf_energy_float,
            'ccsd_correlation_energy': ccsd_corr_energy_float,
            'scf_energy': ccsd_t_total_energy_float if self.calculate_ccsd_t else ccsd_total_energy_float,  # Use CCSD(T) energy if calculated
            'ccsd_total_energy': ccsd_total_energy_float
        }

        # Add CCSD(T) specific results if calculated
        if self.calculate_ccsd_t:
            results.update({
                'ccsd_t_correction': ccsd_t_correction_float,
                'ccsd_t_total_energy': ccsd_t_total_energy_float
            })

        # Add common additional properties (from base calculator)
        common_props = self._extract_common_additional_properties()
        results.update(common_props)

        # Add CCSD diagnostic indicators
        ccsd_diag = self._extract_ccsd_diagnostics()
        results.update(ccsd_diag)

        return results
    
    def _create_scf_method(self, mol):
        """Create HF method object for CCSD reference (RHF/UHF)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = scf.RHF(mol)
            logger.info("Using Restricted Hartree-Fock (RHF) reference for RCCSD")
        else:
            mf = scf.UHF(mol)
            logger.info("Using Unrestricted Hartree-Fock (UHF) reference for UCCSD")
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to HF reference method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UHF' if spin > 0 else 'RHF'} (CCSD reference)"