"""CCSD calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional
from pyscf import gto, scf, cc
from pyscf.geomopt import geometric_solver

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects

logger = logging.getLogger(__name__)


class CCSDCalculator(BaseCalculator):
    """CCSD/CCSD(T) calculator using PySCF for structure optimization and CCSD calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None
        self.ccsd: Optional[cc.CCSD] = None
        self.optimized_geometry: Optional[np.ndarray] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        self.calculate_ccsd_t = False  # Flag to indicate if CCSD(T) should be calculated
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup CCSD calculation with molecular geometry and parameters."""
        try:
            # Extract calculation parameters
            basis = kwargs.get('basis', 'cc-pVDZ')  # CCSD typically uses correlation-consistent basis sets
            charge = kwargs.get('charge', 0)
            spin = kwargs.get('spin', 0)
            max_cycle = kwargs.get('max_cycle', 150)
            solvent_method = kwargs.get('solvent_method', 'none')
            solvent = kwargs.get('solvent', '-')
            memory_mb = kwargs.get('memory_mb', 4000)  # CCSD needs more memory than HF/DFT
            frozen_core = kwargs.get('frozen_core', True)  # Enable frozen core by default for CCSD
            ccsd_t = kwargs.get('ccsd_t', False)  # Whether to calculate CCSD(T) correction

            # Store CCSD(T) flag
            self.calculate_ccsd_t = ccsd_t

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
                self.mol.max_memory = 4000  # CCSDはより多くのメモリが必要
            
            # Setup HF calculation first (CCSD requires HF reference)
            # For closed-shell systems (spin=0), use RHF
            # For open-shell systems (spin>0), use UHF
            if spin == 0:
                self.mf = scf.RHF(self.mol)
                logger.info("Using Restricted Hartree-Fock (RHF) reference for RCCSD")
            else:
                self.mf = scf.UHF(self.mol)
                logger.info("Using Unrestricted Hartree-Fock (UHF) reference for UCCSD")
            
            # Apply solvent effects if requested
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.max_cycle = max_cycle
            
            # Store parameters
            method_name = 'UCCSD' if spin > 0 else 'RCCSD'
            if ccsd_t:
                method_name += '(T)'
            
            self.results.update({
                'basis': basis,
                'charge': charge,
                'spin_multiplicity': 2 * spin + 1,
                'max_cycle': max_cycle,
                'solvent_method': solvent_method,
                'solvent': solvent,
                'atom_count': len(atoms),
                'method': method_name,
                'frozen_core': frozen_core,
                'ccsd_t_correction': ccsd_t
            })
            
        except Exception as e:
            raise InputError(f"Failed to setup CCSD calculation: {str(e)}")
    
    def run_calculation(self) -> Dict[str, Any]:
        """Run structure optimization and CCSD/CCSD(T) energy calculation."""
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
            
            # Step 3: CCSD calculation
            logger.info("Starting CCSD calculation...")
            
            # Create CCSD object based on HF reference
            self.ccsd = cc.CCSD(self.mf)
            
            # Apply frozen core if requested
            if self.results.get('frozen_core', True):
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
            
            # Step 4: CCSD(T) calculation if requested
            ccsd_t_correction = 0.0
            ccsd_t_total_energy = ccsd_total_energy
            
            if self.calculate_ccsd_t:
                logger.info("Starting CCSD(T) perturbative triples correction...")
                ccsd_t_correction = self.ccsd.ccsd_t()
                ccsd_t_total_energy = ccsd_total_energy + ccsd_t_correction
                
                logger.info("CCSD(T) calculation completed successfully.")
                logger.info(f"CCSD(T) triples correction: {ccsd_t_correction} Hartree")
                logger.info(f"CCSD(T) total energy: {ccsd_t_total_energy} Hartree")
            
            # Step 5: Orbital analysis using HF orbitals
            homo_idx, lumo_idx = self._analyze_orbitals()
            
            # Step 6: Prepare results
            chk_path = self.get_checkpoint_path()
            
            # 安全にエネルギーを変換
            if hf_energy is None:
                raise CalculationError("HF energy is None - calculation may have failed")
            if ccsd_corr_energy is None:
                raise CalculationError("CCSD correlation energy is None - calculation may have failed")
            if ccsd_total_energy is None:
                raise CalculationError("CCSD total energy is None - calculation may have failed")
            
            try:
                hf_energy_float = float(hf_energy)
                ccsd_corr_energy_float = float(ccsd_corr_energy)
                ccsd_total_energy_float = float(ccsd_total_energy)
                ccsd_t_correction_float = float(ccsd_t_correction)
                ccsd_t_total_energy_float = float(ccsd_t_total_energy)
            except (ValueError, TypeError) as e:
                raise CalculationError(f"Failed to convert CCSD energies to float: {e}")
            
            # Update results with CCSD data
            self.results.update({
                'hf_energy': hf_energy_float,
                'ccsd_correlation_energy': ccsd_corr_energy_float,
                'scf_energy': ccsd_total_energy_float,  # For consistency with other calculators
                'ccsd_total_energy': ccsd_total_energy_float,
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
            
            # Add CCSD(T) specific results if calculated
            if self.calculate_ccsd_t:
                self.results.update({
                    'ccsd_t_correction': ccsd_t_correction_float,
                    'ccsd_t_total_energy': ccsd_t_total_energy_float
                })
            
            if self.keep_files:
                self.file_manager.save_calculation_results(self.working_dir, self.results)
                self.file_manager.save_geometry(self.working_dir, self.results['optimized_geometry'])
                logger.info(f"Calculation files saved to: {self.working_dir}")
            
            return self.results
            
        except ConvergenceError:
            raise
        except Exception as e:
            raise CalculationError(f"CCSD calculation failed: {str(e)}")