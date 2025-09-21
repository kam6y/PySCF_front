"""MP2 calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, scf, mp
# geometric_solver is now imported in BaseCalculator

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects
from .config_manager import get_memory_for_method

logger = logging.getLogger(__name__)


class MP2Calculator(BaseCalculator):
    """MP2 calculator using PySCF for structure optimization and MP2 energy calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None, optimize_geometry: bool = True):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir, optimize_geometry)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None
        self.mp2: Optional[mp.MP2] = None
        self.optimized_geometry: Optional[np.ndarray] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup MP2 calculation using the base template method."""
        # Call the base template method which handles common setup
        super().setup_calculation(atoms, **kwargs)
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """Validate MP2-specific parameters."""
        # MP2-specific parameters (MP2 has no method-specific parameters beyond common ones)
        # MP2 uses HF as reference, so method is based on spin
        return {
            'method': 'UMP2' if kwargs.get('spin', 0) > 0 else 'RMP2'
        }
    
    def _get_default_memory_mb(self) -> int:
        """Get default memory setting for MP2 calculations from config."""
        return get_memory_for_method('MP2')
    
    def _get_calculation_method_name(self) -> str:
        """Get the name of the calculation method for logging."""
        return 'MP2'
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform MP2-specific calculation after HF."""
        # 安全にベースHFエネルギーを変換
        if base_energy is None:
            raise CalculationError("HF energy is None - calculation may have failed")
        try:
            hf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert HF energy to float: {e}")
        
        logger.info(f"HF energy: {hf_energy_float} Hartree")
        
        # MP2 calculation
        logger.info("Starting MP2 calculation...")
        
        # Create MP2 object based on HF reference
        self.mp2 = mp.MP2(self.mf)
        
        # Run MP2 calculation
        self.mp2.kernel()
        
        # Get MP2 results
        if not hasattr(self.mp2, 'e_corr') or self.mp2.e_corr is None:
            raise CalculationError("MP2 calculation failed: correlation energy not available")
        
        mp2_corr_energy = self.mp2.e_corr
        mp2_total_energy = self.mp2.e_tot
        
        logger.info("MP2 calculation completed successfully.")
        logger.info(f"MP2 correlation energy: {mp2_corr_energy} Hartree")
        logger.info(f"MP2 total energy: {mp2_total_energy} Hartree")
        
        # 安全にMP2エネルギーを変換
        if mp2_corr_energy is None:
            raise CalculationError("MP2 correlation energy is None - calculation may have failed")
        if mp2_total_energy is None:
            raise CalculationError("MP2 total energy is None - calculation may have failed")
        
        try:
            mp2_corr_energy_float = float(mp2_corr_energy)
            mp2_total_energy_float = float(mp2_total_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert MP2 energies to float: {e}")
        
        return {
            'hf_energy': hf_energy_float,
            'mp2_correlation_energy': mp2_corr_energy_float,
            'scf_energy': mp2_total_energy_float,  # For consistency with other calculators
            'mp2_total_energy': mp2_total_energy_float
        }
    
    def _create_scf_method(self, mol):
        """Create HF method object for MP2 reference (RHF/UHF)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = scf.RHF(mol)
            logger.info("Using Restricted Hartree-Fock (RHF) reference for RMP2")
        else:
            mf = scf.UHF(mol)
            logger.info("Using Unrestricted Hartree-Fock (UHF) reference for UMP2")
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to HF reference method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UHF' if spin > 0 else 'RHF'} (MP2 reference)"
    
    def _perform_geometry_optimization(self) -> None:
        """Perform geometry optimization using MP2 level of theory."""
        from pyscf.geomopt import geometric_solver
        from pyscf import mp
        
        logger.info("Starting MP2 geometry optimization...")
        
        # First, run initial HF calculation to get reference
        logger.info("Running initial HF calculation for MP2 reference...")
        hf_energy = self.mf.kernel()
        
        if not self.mf.converged:
            from .exceptions import ConvergenceError
            raise ConvergenceError("Initial HF calculation failed to converge for MP2 geometry optimization")
        
        logger.info(f"Initial HF energy: {hf_energy} Hartree")
        
        # Create MP2 object for geometry optimization
        logger.info("Creating MP2 object for geometry optimization...")
        mp2_obj = mp.MP2(self.mf)
        
        # Perform MP2 geometry optimization
        logger.info("Performing geometry optimization at MP2 level...")
        optimized_mol = geometric_solver.optimize(mp2_obj)
        self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
        logger.info("MP2 geometry optimization completed")
        
        # Apply coordinate alignment after optimization
        self._align_optimized_geometry()
        
        # Store the MP2 object for later use if needed
        self.mp2 = mp2_obj
    
