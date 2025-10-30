"""DFT calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, dft
# geometric_solver is now imported in BaseCalculator

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects
from .config_manager import get_memory_for_method

logger = logging.getLogger(__name__)



class DFTCalculator(BaseCalculator):
    """DFT calculator using PySCF for structure optimization and orbital analysis."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None, optimize_geometry: bool = True):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir, optimize_geometry)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[dft.RKS] = None
        self.optimized_geometry: Optional[np.ndarray] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup DFT calculation using the base template method."""
        # Call the base template method which handles common setup
        super().setup_calculation(atoms, **kwargs)
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """Validate DFT-specific parameters."""
        # DFT-specific parameters
        xc = kwargs.get('xc', 'B3LYP')  # Exchange-correlation functional
        
        # Store parameters for template method access
        self.xc_functional = xc
        
        return {
            'xc_functional': xc,
            'method': 'UKS' if kwargs.get('spin', 0) > 0 else 'RKS'
        }
    
    def _get_default_memory_mb(self) -> int:
        """Get default memory setting for DFT calculations from config."""
        return get_memory_for_method('DFT')
    
    def _get_calculation_method_name(self) -> str:
        """Get the name of the calculation method for logging."""
        return 'DFT'
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform DFT-specific calculation (just return SCF energy)."""
        # 安全にエネルギーを変換
        if base_energy is None:
            raise CalculationError("SCF energy is None - calculation may have failed")
        try:
            scf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert SCF energy to float: {e}")
        
        return {'scf_energy': scf_energy_float}
    
    def _create_scf_method(self, mol):
        """Create DFT method object (RKS/UKS)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = dft.RKS(mol)
            logger.info("Using Restricted Kohn-Sham (RKS) for closed-shell system")
        else:
            mf = dft.UKS(mol)
            logger.info("Using Unrestricted Kohn-Sham (UKS) for open-shell system")
        
        # Set XC functional
        mf.xc = self.xc_functional
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to DFT method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return 'UKS' if spin > 0 else 'RKS'
    
