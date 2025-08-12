"""TDDFT calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from pyscf import gto, dft, tddft, tdscf

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects

logger = logging.getLogger(__name__)


class TDDFTCalculator(BaseCalculator):
    """TDDFT calculator using PySCF for excited state calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[dft.RKS] = None
        self.mytd: Optional[tddft.TDDFT] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup TDDFT calculation with molecular geometry and parameters."""
        try:
            # Extract calculation parameters
            basis = kwargs.get('basis', '6-31G(d)')
            xc = kwargs.get('xc', 'B3LYP')
            charge = kwargs.get('charge', 0)
            spin = kwargs.get('spin', 0)
            max_cycle = kwargs.get('max_cycle', 150)
            solvent_method = kwargs.get('solvent_method', 'none')
            solvent = kwargs.get('solvent', '-')
            
            # TDDFT-specific parameters
            nstates = kwargs.get('nstates', 10)
            tddft_method = kwargs.get('tddft_method', 'TDDFT')
            analyze_nto = kwargs.get('analyze_nto', False)
            
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
            
            # Setup DFT calculation (ground state) first
            self.mf = dft.RKS(self.mol)
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
                'tddft_nstates': nstates,
                'tddft_method': tddft_method,
                'tddft_analyze_nto': analyze_nto
            })
            
        except Exception as e:
            raise InputError(f"Failed to setup TDDFT calculation: {str(e)}")
    
    def run_calculation(self) -> Dict[str, Any]:
        """Run ground state DFT followed by TDDFT calculation."""
        if self.mol is None or self.mf is None:
            raise CalculationError("Calculation not properly setup. Call setup_calculation first.")
        
        try:
            # Step 1: Ground state DFT calculation
            logger.info("Starting ground state DFT calculation...")
            scf_energy = self.mf.kernel()
            
            # Check SCF convergence
            if not self.mf.converged:
                raise ConvergenceError("Ground state SCF calculation failed to converge")
                
            logger.info(f"Ground state DFT completed. SCF energy: {scf_energy:.6f} hartree")
            
            # Step 2: TDDFT calculation
            logger.info("Starting TDDFT calculation...")
            nstates = self.results.get('tddft_nstates', 10)
            tddft_method = self.results.get('tddft_method', 'TDDFT')
            
            if tddft_method == 'TDA':
                # Tamm-Dancoff approximation
                self.mytd = tdscf.TDA(self.mf)
            else:
                # Full TDDFT
                self.mytd = tddft.TDDFT(self.mf)
            
            self.mytd.nstates = nstates
            excitation_energies = self.mytd.kernel()
            
            if not hasattr(self.mytd, 'e') or self.mytd.e is None:
                raise CalculationError("TDDFT calculation failed: no excitation energies obtained")
            
            logger.info(f"TDDFT calculation completed. Found {len(self.mytd.e)} excited states.")
            
            # Step 3: Analyze results
            tddft_results = self._analyze_tddft_results()
            
            # Step 4: Orbital analysis for ground state
            homo_idx, lumo_idx = self._analyze_orbitals()
            
            # Step 5: Prepare results
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
                **tddft_results
            })
            
            if self.keep_files:
                self.file_manager.save_calculation_results(self.working_dir, self.results)
                logger.info(f"TDDFT calculation files saved to: {self.working_dir}")
            
            return self.results
            
        except ConvergenceError:
            raise
        except Exception as e:
            logger.error(f"TDDFT calculation failed: {str(e)}")
            raise CalculationError(f"TDDFT calculation failed: {str(e)}")
    
    def _analyze_tddft_results(self) -> Dict[str, Any]:
        """Analyze TDDFT results and extract key information."""
        if self.mytd is None or not hasattr(self.mytd, 'e'):
            raise CalculationError("TDDFT results not available")
        
        # Convert excitation energies from hartree to eV
        ev_to_hartree = 27.2114
        excitation_energies_ev = [float(e * ev_to_hartree) for e in self.mytd.e]
        
        # Convert to wavelengths in nm
        # E(eV) = hc/λ, where hc ≈ 1239.84 eV·nm
        excitation_wavelengths = [1239.84 / e if e > 0 else 0 for e in excitation_energies_ev]
        
        # Get oscillator strengths
        oscillator_strengths = []
        transition_dipoles = []
        
        try:
            # Analyze transitions to get oscillator strengths and dipole moments
            self.mytd.analyze()
            
            if hasattr(self.mytd, 'oscillator_strength'):
                oscillator_strengths = [float(f) for f in self.mytd.oscillator_strength]
            
            if hasattr(self.mytd, 'transition_dipole'):
                for dipole in self.mytd.transition_dipole:
                    if len(dipole) >= 3:
                        transition_dipoles.append({
                            'x': float(dipole[0]),
                            'y': float(dipole[1]),
                            'z': float(dipole[2])
                        })
            
        except Exception as e:
            logger.warning(f"Failed to analyze transitions: {e}")
        
        # Analyze major transitions
        major_transitions = self._analyze_major_transitions(
            excitation_energies_ev, excitation_wavelengths, oscillator_strengths
        )
        
        return {
            'excitation_energies': excitation_energies_ev,
            'excitation_wavelengths': excitation_wavelengths,
            'oscillator_strengths': oscillator_strengths,
            'transition_dipoles': transition_dipoles,
            'major_transitions': major_transitions
        }
    
    def _analyze_major_transitions(
        self, 
        energies: List[float], 
        wavelengths: List[float], 
        oscillator_strengths: List[float]
    ) -> List[Dict[str, Any]]:
        """Analyze and identify major orbital transitions."""
        major_transitions = []
        
        for i, (energy, wavelength) in enumerate(zip(energies, wavelengths)):
            osc_strength = oscillator_strengths[i] if i < len(oscillator_strengths) else 0.0
            
            # Determine dominant transition type based on energy and oscillator strength
            dominant_transition = self._identify_transition_type(energy, osc_strength)
            
            major_transitions.append({
                'state': i + 1,
                'energy': energy,
                'wavelength': wavelength,
                'oscillator_strength': osc_strength,
                'dominant_transition': dominant_transition
            })
        
        return major_transitions
    
    def _identify_transition_type(self, energy_ev: float, osc_strength: float) -> str:
        """Identify the likely transition type based on energy and oscillator strength."""
        if energy_ev > 10.0:
            return "Core/High energy transition"
        elif energy_ev > 6.0:
            if osc_strength > 0.1:
                return "π→π* transition"
            else:
                return "n→π* transition"  
        elif energy_ev > 3.0:
            if osc_strength > 0.05:
                return "HOMO→LUMO (π→π*)"
            else:
                return "HOMO→LUMO (n→π*)"
        else:
            return "Low energy transition"
    
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
    
    def _analyze_orbitals(self) -> Tuple[int, int]:
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