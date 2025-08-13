"""DFT calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, dft
from pyscf.geomopt import geometric_solver

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
            
            # Store parameters
            self.results.update({
                'basis': basis,
                'xc_functional': xc,
                'charge': charge,
                'spin_multiplicity': 2 * spin + 1,
                'max_cycle': max_cycle,
                'solvent_method': solvent_method,
                'solvent': solvent,
                'atom_count': len(atoms),
                'method': 'UKS' if spin > 0 else 'RKS'
            })
            
        except Exception as e:
            raise InputError(f"Failed to setup DFT calculation: {str(e)}")
    
    def run_calculation(self) -> Dict[str, Any]:
        """Run structure optimization and orbital analysis."""
        if self.mol is None or self.mf is None:
            raise CalculationError("Calculation not properly setup. Call setup_calculation first.")
        
        try:
            # Step 1: Structure optimization
            logger.info("Starting geometry optimization...")
            optimized_mol = geometric_solver.optimize(self.mf)
            self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
            
            # Step 2: SCF calculation with optimized geometry
            logger.info("Running SCF calculation with optimized geometry...")
            
            # Recreate mf object with optimized geometry while preserving solvent effects
            solvent_method = self.results.get('solvent_method', 'none')
            solvent = self.results.get('solvent', '-')
            spin = (self.results.get('spin_multiplicity', 1) - 1) // 2
            
            if spin == 0:
                self.mf = dft.RKS(optimized_mol)
            else:
                self.mf = dft.UKS(optimized_mol)
            
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.xc = self.results['xc_functional']
            self.mf.max_cycle = self.results['max_cycle']
            
            scf_energy = self.mf.kernel()
            
            # Check SCF convergence
            if not self.mf.converged:
                raise ConvergenceError("SCF calculation failed to converge")
            
            # Verify orbital occupations
            if self.mf.mo_occ is None or len(self.mf.mo_occ) == 0:
                raise CalculationError("SCF calculation failed: mo_occ not properly assigned")
            
            logger.info("SCF calculation completed successfully.")
            logger.info(f"Number of occupied orbitals: {self._count_occupied_orbitals()}")
            
            # Step 3: Orbital analysis
            homo_idx, lumo_idx = self._analyze_orbitals()
            
            # Step 4: Prepare results
            chk_path = self.get_checkpoint_path()
            
            # 安全にエネルギーを変換
            if scf_energy is None:
                raise CalculationError("SCF energy is None - calculation may have failed")
            try:
                scf_energy_float = float(scf_energy)
            except (ValueError, TypeError) as e:
                raise CalculationError(f"Failed to convert SCF energy to float: {e}")
            
            self.results.update({
                'scf_energy': scf_energy_float,
                'converged': True,
                'homo_index': homo_idx,
                'lumo_index': lumo_idx,
                'num_occupied_orbitals': int(self._count_occupied_orbitals()),
                'num_virtual_orbitals': int(self._count_virtual_orbitals()),
                'checkpoint_file': chk_path,
                'checkpoint_exists': os.path.exists(chk_path),
                'working_directory': self.working_dir,
                'optimized_geometry': self._geometry_to_xyz_string()
            })
            
            if self.keep_files:
                self.file_manager.save_calculation_results(self.working_dir, self.results)
                self.file_manager.save_geometry(self.working_dir, self.results['optimized_geometry'])
                logger.info(f"Calculation files saved to: {self.working_dir}")
            
            return self.results
            
        except ConvergenceError:
            raise
        except Exception as e:
            raise CalculationError(f"Calculation failed: {str(e)}")
    
    def _atoms_to_string(self, atoms: List[List]) -> str:
        """Convert atoms list to PySCF atom string format."""
        atom_lines = []
        for atom_data in atoms:
            symbol = atom_data[0]
            coords = atom_data[1]
            if len(coords) != 3:
                raise GeometryError(f"Invalid coordinates for atom {symbol}: expected 3 coordinates, got {len(coords)}")
            atom_lines.append(f"{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        return "\n".join(atom_lines)
    
    def _analyze_orbitals(self) -> tuple[int, int]:
        """Analyze molecular orbitals to find HOMO and LUMO indices."""
        if self.mf.mo_occ is None:
            raise CalculationError("Orbital occupations not available")
        
        # Handle both RKS (1D array) and UKS (2D array) cases
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS case: use alpha orbitals
            mo_occ = mo_occ[0]
        
        # Find HOMO (highest occupied molecular orbital)
        occupied_indices = np.where(mo_occ > 0)[0]
        if len(occupied_indices) == 0:
            raise CalculationError("No occupied orbitals found")
        homo_idx = occupied_indices[-1]
        
        # Find LUMO (lowest unoccupied molecular orbital)
        unoccupied_indices = np.where(mo_occ == 0)[0]
        if len(unoccupied_indices) == 0:
            raise CalculationError("No unoccupied orbitals found")
        lumo_idx = unoccupied_indices[0]
        
        return int(homo_idx), int(lumo_idx)
    
    def _count_occupied_orbitals(self) -> int:
        """Count the number of occupied orbitals."""
        if self.mf.mo_occ is None:
            return 0
        
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS case: count both alpha and beta orbitals
            return int(np.sum(mo_occ > 0))
        else:
            # RKS case: simple sum
            return int(np.sum(mo_occ > 0))
    
    def _count_virtual_orbitals(self) -> int:
        """Count the number of virtual orbitals."""
        if self.mf.mo_occ is None:
            return 0
        
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS case: count both alpha and beta orbitals
            return int(np.sum(mo_occ == 0))
        else:
            # RKS case: simple sum
            return int(np.sum(mo_occ == 0))
    
    def _geometry_to_xyz_string(self) -> str:
        """Convert optimized geometry to XYZ format string."""
        if self.optimized_geometry is None:
            return ""
        
        atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
        lines = [str(self.mol.natm)]
        lines.append("Optimized geometry from PySCF DFT calculation")
        
        for i, (symbol, coords) in enumerate(zip(atom_symbols, self.optimized_geometry)):
            lines.append(f"{symbol:2s} {coords[0]:12.6f} {coords[1]:12.6f} {coords[2]:12.6f}")
        
        return "\n".join(lines)