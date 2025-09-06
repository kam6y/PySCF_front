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

logger = logging.getLogger(__name__)



class DFTCalculator(BaseCalculator):
    """DFT calculator using PySCF for structure optimization and orbital analysis."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[dft.RKS] = None
        self.optimized_geometry: Optional[np.ndarray] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup DFT calculation with molecular geometry and parameters."""
        try:
            # Extract calculation parameters
            basis = kwargs.get('basis', '6-31G(d)')
            xc = kwargs.get('xc', 'B3LYP')
            charge = kwargs.get('charge', 0)
            spin = kwargs.get('spin', 0)
            max_cycle = kwargs.get('max_cycle', 150)
            solvent_method = kwargs.get('solvent_method', 'none')
            solvent = kwargs.get('solvent', '-')
            memory_mb = kwargs.get('memory_mb', 2000)  # Default 2GB

            # Convert atoms list to PySCF format
            atom_string = self._atoms_to_string(atoms)
            
            # Create molecular object
            self.mol = gto.M(
                atom=atom_string,
                basis=basis,
                charge=charge,
                spin=spin,
                verbose=0
            )
            # 安全なメモリ設定を適用
            if memory_mb and memory_mb > 0:
                self.mol.max_memory = memory_mb
            else:
                self.mol.max_memory = 2000  # デフォルト2GB
            
            # Setup DFT calculation based on spin multiplicity
            # For closed-shell systems (spin=0), use RKS
            # For open-shell systems (spin>0), use UKS
            if spin == 0:
                self.mf = dft.RKS(self.mol)
                logger.info("Using Restricted Kohn-Sham (RKS) for closed-shell system")
            else:
                self.mf = dft.UKS(self.mol)
                logger.info("Using Unrestricted Kohn-Sham (UKS) for open-shell system")
            
            # Apply solvent effects if requested
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.xc = xc
            self.mf.max_cycle = max_cycle
            
            # Store parameters for template method
            self.max_cycle = max_cycle
            self.xc_functional = xc
            self.solvent_method = solvent_method
            self.solvent = solvent
            
            # Store parameters
            self.results.update({
                'basis': basis,
                'xc_functional': xc,
                'charge': charge,
                'spin': spin,
                'max_cycle': max_cycle,
                'solvent_method': solvent_method,
                'solvent': solvent,
                'atom_count': len(atoms),
                'method': 'UKS' if spin > 0 else 'RKS'
            })
            
        except Exception as e:
            raise InputError(f"Failed to setup DFT calculation: {str(e)}")
    
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
    
