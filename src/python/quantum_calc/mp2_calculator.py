"""MP2 calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, scf, mp
from pyscf.geomopt import geometric_solver

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects

logger = logging.getLogger(__name__)


class MP2Calculator(BaseCalculator):
    """MP2 calculator using PySCF for structure optimization and MP2 energy calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None
        self.mp2: Optional[mp.MP2] = None
        self.optimized_geometry: Optional[np.ndarray] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup MP2 calculation with molecular geometry and parameters."""
        try:
            # Extract calculation parameters
            basis = kwargs.get('basis', '6-31G(d)')
            charge = kwargs.get('charge', 0)
            spin = kwargs.get('spin', 0)
            max_cycle = kwargs.get('max_cycle', 150)
            solvent_method = kwargs.get('solvent_method', 'none')
            solvent = kwargs.get('solvent', '-')
            cpu_cores = kwargs.get('cpu_cores')
            memory_mb = kwargs.get('memory_mb')
            
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
            
            # Apply resource settings
            self.apply_resource_settings(self.mol, memory_mb, cpu_cores)
            
            # Setup HF calculation first (MP2 requires HF reference)
            # For closed-shell systems (spin=0), use RHF
            # For open-shell systems (spin>0), use UHF
            if spin == 0:
                self.mf = scf.RHF(self.mol)
                logger.info("Using Restricted Hartree-Fock (RHF) reference for RMP2")
            else:
                self.mf = scf.UHF(self.mol)
                logger.info("Using Unrestricted Hartree-Fock (UHF) reference for UMP2")
            
            # Apply solvent effects if requested
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.max_cycle = max_cycle
            
            # Store parameters
            self.results.update({
                'basis': basis,
                'charge': charge,
                'spin_multiplicity': 2 * spin + 1,
                'max_cycle': max_cycle,
                'solvent_method': solvent_method,
                'solvent': solvent,
                'atom_count': len(atoms),
                'method': 'UMP2' if spin > 0 else 'RMP2',
                'cpu_cores': cpu_cores,
                'memory_mb': memory_mb
            })
            
        except Exception as e:
            raise InputError(f"Failed to setup MP2 calculation: {str(e)}")
    
    def run_calculation(self) -> Dict[str, Any]:
        """Run structure optimization and MP2 energy calculation."""
        if self.mol is None or self.mf is None:
            raise CalculationError("Calculation not properly setup. Call setup_calculation first.")
        
        try:
            # Step 1: Structure optimization using HF method
            logger.info("Starting geometry optimization using HF method...")
            optimized_mol = geometric_solver.optimize(self.mf)
            self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
            
            # Step 2: HF calculation with optimized geometry
            logger.info("Running HF calculation with optimized geometry...")
            
            # Recreate mf object with optimized geometry while preserving solvent effects
            solvent_method = self.results.get('solvent_method', 'none')
            solvent = self.results.get('solvent', '-')
            spin = (self.results.get('spin_multiplicity', 1) - 1) // 2
            
            # Apply resource settings to optimized molecule
            memory_mb = self.results.get('memory_mb')
            cpu_cores = self.results.get('cpu_cores')
            self.apply_resource_settings(optimized_mol, memory_mb, cpu_cores)
            
            if spin == 0:
                self.mf = scf.RHF(optimized_mol)
            else:
                self.mf = scf.UHF(optimized_mol)
            
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.max_cycle = self.results['max_cycle']
            
            hf_energy = self.mf.kernel()
            
            # Check HF convergence
            if not self.mf.converged:
                raise ConvergenceError("HF calculation failed to converge")
            
            logger.info("HF calculation completed successfully.")
            logger.info(f"HF energy: {hf_energy} Hartree")
            
            # Step 3: MP2 calculation
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
            
            # Step 4: Orbital analysis using HF orbitals
            homo_idx, lumo_idx = self._analyze_orbitals()
            
            # Step 5: Prepare results
            chk_path = self.get_checkpoint_path()
            self.results.update({
                'hf_energy': float(hf_energy),
                'mp2_correlation_energy': float(mp2_corr_energy),
                'scf_energy': float(mp2_total_energy),  # For consistency with other calculators
                'mp2_total_energy': float(mp2_total_energy),
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
            raise CalculationError(f"MP2 calculation failed: {str(e)}")
    
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
        """Analyze molecular orbitals to find HOMO and LUMO indices using HF orbitals."""
        if self.mf.mo_occ is None:
            raise CalculationError("Orbital occupations not available")
        
        # Handle both RHF (1D array) and UHF (2D array) cases
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UHF case: use alpha orbitals
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
            # UHF case: count both alpha and beta orbitals
            return int(np.sum(mo_occ > 0))
        else:
            # RHF case: simple sum
            return int(np.sum(mo_occ > 0))
    
    def _count_virtual_orbitals(self) -> int:
        """Count the number of virtual orbitals."""
        if self.mf.mo_occ is None:
            return 0
        
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UHF case: count both alpha and beta orbitals
            return int(np.sum(mo_occ == 0))
        else:
            # RHF case: simple sum
            return int(np.sum(mo_occ == 0))
    
    def _geometry_to_xyz_string(self) -> str:
        """Convert optimized geometry to XYZ format string."""
        if self.optimized_geometry is None:
            return ""
        
        atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
        lines = [str(self.mol.natm)]
        lines.append("Optimized geometry from PySCF MP2 calculation")
        
        for i, (symbol, coords) in enumerate(zip(atom_symbols, self.optimized_geometry)):
            lines.append(f"{symbol:2s} {coords[0]:12.6f} {coords[1]:12.6f} {coords[2]:12.6f}")
        
        return "\n".join(lines)