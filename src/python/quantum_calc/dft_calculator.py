"""DFT calculator implementation using PySCF."""

import os
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, dft
from pyscf.geomopt import geometric_solver

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager


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
            
            # Setup DFT calculation
            self.mf = dft.RKS(self.mol)
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
                'atom_count': len(atoms)
            })
            
        except Exception as e:
            raise InputError(f"Failed to setup DFT calculation: {str(e)}")
    
    def run_calculation(self) -> Dict[str, Any]:
        """Run structure optimization and orbital analysis."""
        if self.mol is None or self.mf is None:
            raise CalculationError("Calculation not properly setup. Call setup_calculation first.")
        
        try:
            # Step 1: Structure optimization
            print("Starting geometry optimization...")
            optimized_mol = geometric_solver.optimize(self.mf)
            self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
            
            # Step 2: SCF calculation with optimized geometry
            print("Running SCF calculation with optimized geometry...")
            self.mf = dft.RKS(optimized_mol)
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
            
            print("SCF calculation completed successfully.")
            print(f"Number of occupied orbitals: {sum(self.mf.mo_occ > 0)}")
            
            # Step 3: Orbital analysis
            homo_idx, lumo_idx = self._analyze_orbitals()
            
            # Step 4: Prepare results
            chk_path = self.get_checkpoint_path()
            self.results.update({
                'scf_energy': float(scf_energy),
                'converged': True,
                'homo_index': homo_idx,
                'lumo_index': lumo_idx,
                'num_occupied_orbitals': int(sum(self.mf.mo_occ > 0)),
                'num_virtual_orbitals': int(sum(self.mf.mo_occ == 0)),
                'checkpoint_file': chk_path,
                'checkpoint_exists': os.path.exists(chk_path),
                'working_directory': self.working_dir,
                'optimized_geometry': self._geometry_to_xyz_string()
            })
            
            # Save calculation files and information
            if self.keep_files:
                self.file_manager.save_calculation_info(self.working_dir, self.results)
                self.file_manager.save_geometry(self.working_dir, self.results['optimized_geometry'])
                print(f"Calculation files saved to: {self.working_dir}")
            
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
        
        # Find HOMO (highest occupied molecular orbital)
        occupied_indices = np.where(self.mf.mo_occ > 0)[0]
        if len(occupied_indices) == 0:
            raise CalculationError("No occupied orbitals found")
        homo_idx = occupied_indices[-1]
        
        # Find LUMO (lowest unoccupied molecular orbital)
        unoccupied_indices = np.where(self.mf.mo_occ == 0)[0]
        if len(unoccupied_indices) == 0:
            raise CalculationError("No unoccupied orbitals found")
        lumo_idx = unoccupied_indices[0]
        
        return int(homo_idx), int(lumo_idx)
    
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