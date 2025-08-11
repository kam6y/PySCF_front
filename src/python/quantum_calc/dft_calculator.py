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
            
            # Setup DFT calculation with solvent effects
            if solvent_method == 'none' or solvent == '-':
                self.mf = dft.RKS(self.mol)
            elif solvent_method.lower() in ['pcm', 'ief-pcm', 'c-pcm', 'cosmo']:
                self.mf = dft.RKS(self.mol).PCM()
                self._setup_pcm_solvent(solvent_method, solvent)
            else:
                logger.warning(f"Unsupported solvent method: {solvent_method}, using no solvent")
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
                'solvent_method': solvent_method,
                'solvent': solvent,
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
            logger.info("Starting geometry optimization...")
            optimized_mol = geometric_solver.optimize(self.mf)
            self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
            
            # Step 2: SCF calculation with optimized geometry
            logger.info("Running SCF calculation with optimized geometry...")
            
            # Recreate mf object with optimized geometry while preserving solvent effects
            solvent_method = self.results.get('solvent_method', 'none')
            solvent = self.results.get('solvent', '-')
            
            if solvent_method == 'none' or solvent == '-':
                self.mf = dft.RKS(optimized_mol)
            elif solvent_method.lower() in ['pcm', 'ief-pcm', 'c-pcm', 'cosmo']:
                self.mf = dft.RKS(optimized_mol).PCM()
                self._setup_pcm_solvent(solvent_method, solvent)
            else:
                logger.warning(f"Unsupported solvent method: {solvent_method}, using no solvent")
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
            
            logger.info("SCF calculation completed successfully.")
            logger.info(f"Number of occupied orbitals: {sum(self.mf.mo_occ > 0)}")
            
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
    
    def _setup_pcm_solvent(self, method: str, solvent: str) -> None:
        """Setup PCM solvent effects."""
        # Set PCM method
        method_lower = method.lower()
        if method_lower == 'ief-pcm':
            self.mf.with_solvent.method = 'IEF-PCM'
        elif method_lower == 'c-pcm':
            self.mf.with_solvent.method = 'C-PCM'
        elif method_lower == 'cosmo':
            self.mf.with_solvent.method = 'COSMO'
        else:  # default to IEF-PCM for 'pcm'
            self.mf.with_solvent.method = 'IEF-PCM'
        
        # Set dielectric constant based on solvent
        solvent_dielectric = {
            'water': 78.3553,
            'dimethylsulfoxide': 46.826,
            'dmso': 46.826,  # alias
            'n,n-dimethylformamide': 37.219,
            'dmf': 37.219,  # alias
            'nitromethane': 36.562,
            'methanol': 32.613,
            'ethanol': 24.852,
            'acetone': 20.493,
            'dichloroethane': 10.125,
            'dichloromethane': 8.93,
            'tetrahydrofuran': 7.4297,
            'thf': 7.4297,  # alias
            'chlorobenzene': 5.6968,
            'chloroform': 4.7113,
            'diethylether': 4.2400,
            'toluene': 2.3741,
            'benzene': 2.2706,
            '1,4-dioxane': 2.2099,
            'dioxane': 2.2099,  # alias
            'cyclohexane': 2.0160
        }
        
        if solvent.lower() in solvent_dielectric:
            self.mf.with_solvent.eps = solvent_dielectric[solvent.lower()]
        else:
            # Try to parse as custom dielectric constant
            try:
                eps_value = float(solvent)
                if eps_value > 1.0:
                    self.mf.with_solvent.eps = eps_value
                else:
                    logger.warning(f"Invalid dielectric constant: {solvent}, using water (78.36)")
                    self.mf.with_solvent.eps = 78.3553
            except ValueError:
                logger.warning(f"Unknown solvent: {solvent}, using water dielectric constant")
                self.mf.with_solvent.eps = 78.3553
    
