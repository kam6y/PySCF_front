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
            cpu_cores = kwargs.get('cpu_cores')
            memory_mb = kwargs.get('memory_mb')
            
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
            
            # Apply resource settings
            self.apply_resource_settings(self.mol, memory_mb, cpu_cores)
            
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
                'tddft_analyze_nto': analyze_nto,
                'cpu_cores': cpu_cores,
                'memory_mb': memory_mb
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
            
            # Validate nstates setting
            if nstates <= 0:
                raise InputError(f"Invalid number of excited states: {nstates}. Must be positive.")
            
            # Check if requested number of states is reasonable for the system
            n_orb = len(self.mf.mo_energy) if hasattr(self.mf, 'mo_energy') else 0
            if nstates > n_orb // 2:
                logger.warning(f"Requested {nstates} excited states, but system has only {n_orb} orbitals. "
                             f"Consider reducing nstates to avoid convergence issues.")
            
            # Run TDDFT calculation
            try:
                excitation_energies = self.mytd.kernel()
            except Exception as e:
                error_msg = str(e).lower()
                if "singular" in error_msg or "convergence" in error_msg:
                    raise ConvergenceError(f"TDDFT calculation failed to converge: {str(e)}")
                elif "memory" in error_msg:
                    raise CalculationError(f"TDDFT calculation failed due to insufficient memory: {str(e)}")
                else:
                    raise CalculationError(f"TDDFT calculation failed: {str(e)}")
            
            # Validate results
            if not hasattr(self.mytd, 'e') or self.mytd.e is None:
                raise CalculationError("TDDFT calculation failed: no excitation energies obtained")
                
            if len(self.mytd.e) == 0:
                raise CalculationError("TDDFT calculation failed: no excited states found")
                
            if len(self.mytd.e) < nstates:
                logger.warning(f"Only {len(self.mytd.e)} excited states found, but {nstates} were requested. "
                             f"This may indicate convergence issues.")
            
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
        """
        Analyze TDDFT results and extract key spectroscopic information.
        
        This method processes the TDDFT calculation results to compute:
        - Excitation energies in eV
        - Wavelengths in nm using E = hc/λ relationship  
        - Oscillator strengths using length gauge
        - Transition dipole moments
        - Major orbital transitions with spectroscopic classification
        - Natural Transition Orbital (NTO) analysis if requested
        
        Returns:
            Dict containing all calculated spectroscopic properties
            
        Note:
            Uses PySCF's standard methods without analyze() as per PySCF documentation.
            Oscillator strengths calculated in length gauge for better accuracy.
        """
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
        
        # Calculate oscillator strengths using length gauge (standard approach)
        try:
            osc_strengths_result = self.mytd.oscillator_strength(gauge='length')
            if osc_strengths_result is not None:
                oscillator_strengths = [float(f) for f in osc_strengths_result]
                logger.info(f"Retrieved {len(oscillator_strengths)} oscillator strengths (length gauge)")
        except Exception as osc_error:
            logger.warning(f"Failed to get oscillator strengths: {osc_error}")
        
        # Calculate transition dipole moments
        try:
            dipole_result = self.mytd.transition_dipole()
            if dipole_result is not None:
                for dipole in dipole_result:
                    if hasattr(dipole, '__len__') and len(dipole) >= 3:
                        transition_dipoles.append({
                            'x': float(dipole[0]),
                            'y': float(dipole[1]),
                            'z': float(dipole[2])
                        })
                logger.info(f"Retrieved {len(transition_dipoles)} transition dipole moments")
        except Exception as dipole_error:
            logger.warning(f"Failed to get transition dipole moments: {dipole_error}")
        
        # Analyze major transitions
        major_transitions = self._analyze_major_transitions(
            excitation_energies_ev, excitation_wavelengths, oscillator_strengths
        )
        
        # Perform NTO analysis if requested
        nto_analysis = None
        analyze_nto = self.results.get('tddft_analyze_nto', False)
        logger.info(f"NTO analysis requested: {analyze_nto}")
        
        if analyze_nto:
            logger.info("Starting Natural Transition Orbital (NTO) analysis...")
            try:
                nto_analysis = self._analyze_natural_transition_orbitals(excitation_energies_ev)
                if nto_analysis:
                    logger.info(f"NTO analysis completed successfully for {len(nto_analysis)} excited states")
                    # Log NTO pairs count for each state
                    for i, state_data in enumerate(nto_analysis):
                        pairs_count = len(state_data.get('nto_pairs', []))
                        logger.info(f"State S{state_data.get('state', i+1)}: {pairs_count} NTO pairs found")
                else:
                    logger.warning("NTO analysis returned no results")
            except Exception as e:
                logger.error(f"NTO analysis failed with exception: {e}", exc_info=True)
        else:
            logger.info("NTO analysis skipped (not requested)")
        
        return {
            'excitation_energies': excitation_energies_ev,
            'excitation_wavelengths': excitation_wavelengths,
            'oscillator_strengths': oscillator_strengths,
            'transition_dipoles': transition_dipoles,
            'major_transitions': major_transitions,
            'nto_analysis': nto_analysis
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
        """
        Identify the likely transition type based on wavelength and oscillator strength.
        
        Classification based on literature:
        - π→π* transitions: Allowed, high intensity (f > 0.1), typically 165-400 nm
        - n→π* transitions: Forbidden, low intensity (f < 0.01), typically 275-350 nm  
        - σ→σ* transitions: Very high energy, UV-C region < 200 nm
        
        References:
        - Pavia et al., Introduction to Spectroscopy
        - Skoog et al., Principles of Instrumental Analysis
        """
        # Convert energy to wavelength (nm)
        wavelength_nm = 1239.84 / energy_ev if energy_ev > 0 else float('inf')
        
        # High energy transitions (UV-C region)
        if wavelength_nm < 200:
            return "σ→σ*/n→σ* transition"
        
        # UV-B and UV-A region analysis
        elif wavelength_nm < 280:
            # UV-C to UV-B boundary - likely σ→π* or very high energy π→π*
            if osc_strength > 0.1:
                return "π→π* transition (high energy)"
            else:
                return "n→σ*/σ→π* transition"
        
        elif wavelength_nm < 350:
            # UV-B region - characteristic of n→π* and some π→π* transitions
            if osc_strength > 0.01:
                return "π→π* transition"
            else:
                return "n→π* transition (forbidden)"
        
        elif wavelength_nm < 400:
            # UV-A region - primarily π→π* transitions in conjugated systems
            if osc_strength > 0.01:
                return "π→π* transition (conjugated)"
            else:
                return "n→π* transition (weak)"
        
        elif wavelength_nm < 700:
            # Visible region - extended conjugation π→π*
            return "π→π* transition (extended conjugation)"
        
        else:
            # Near-IR region - very low energy transitions
            return "Low energy transition (>700 nm)"
    
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
    
    def _analyze_natural_transition_orbitals(self, excitation_energies: List[float]) -> List[Dict[str, Any]]:
        """
        Perform Natural Transition Orbital (NTO) analysis for excited states.
        
        NTO analysis provides a compact description of electronic excitations by 
        performing Singular Value Decomposition (SVD) of the transition density matrix:
        T = U * diag(λ) * V†
        
        Where:
        - U contains hole orbitals (where electron is promoted from)
        - V contains particle orbitals (where electron is promoted to)  
        - λ are singular values (weights) representing pair contributions
        
        Args:
            excitation_energies: List of excitation energies in eV
            
        Returns:
            List of dictionaries, one per excited state, containing:
            - state: Excited state number (1-indexed)
            - energy: Excitation energy in eV
            - nto_pairs: List of hole-particle orbital pairs with weights
            - total_nto_pairs: Total number of pairs analyzed
            
        References:
            Martin, R. L., J. Chem. Phys. 118, 4775-4777 (2003)
        """
        if self.mytd is None:
            raise CalculationError("TDDFT calculation not available for NTO analysis")
        
        logger.info(f"Starting NTO analysis for {len(excitation_energies)} excited states")
        nto_results = []
        
        # Check if TDDFT object has get_nto method
        if not hasattr(self.mytd, 'get_nto'):
            logger.error("TDDFT object does not have get_nto method")
            raise CalculationError("TDDFT object does not support NTO analysis (get_nto method not available)")
        
        # Analyze NTOs for each excited state
        for state_idx in range(len(excitation_energies)):
            state_number = state_idx + 1  # PySCF uses 1-indexed states
            logger.info(f"Analyzing NTOs for excited state S{state_number} (index {state_idx})")
            
            try:
                # Get NTO weights and coefficients for this state
                # The get_nto method returns (weights, nto_coefficients)
                # weights: singular values (λ) from SVD decomposition  
                # nto_coefficients: natural transition orbitals in AO basis
                logger.debug(f"Calling mytd.get_nto(state={state_number})")
                nto_result = self.mytd.get_nto(state=state_number)
                
                if nto_result is None or len(nto_result) != 2:
                    logger.warning(f"Invalid NTO result for state S{state_number}: {nto_result}")
                    continue
                    
                weights, nto_coeff = nto_result
                logger.info(f"Successfully retrieved NTO data for state S{state_number}")
                logger.debug(f"Weights type: {type(weights)}, shape: {getattr(weights, 'shape', 'unknown')}")
                logger.debug(f"NTO coeff type: {type(nto_coeff)}, shape: {getattr(nto_coeff, 'shape', 'unknown')}")
                
                # Process NTO data
                nto_pairs = self._process_nto_data(weights, nto_coeff, state_number)
                
                state_result = {
                    'state': state_number,
                    'energy': float(excitation_energies[state_idx]),
                    'nto_pairs': nto_pairs,
                    'total_nto_pairs': len(nto_pairs)
                }
                
                nto_results.append(state_result)
                
            except Exception as e:
                logger.error(f"Failed to analyze NTO for state {state_number}: {type(e).__name__}: {e}", exc_info=True)
                # Add empty result for this state to maintain consistency
                nto_results.append({
                    'state': state_number,
                    'energy': float(excitation_energies[state_idx]),
                    'nto_pairs': [],
                    'total_nto_pairs': 0
                })
        
        return nto_results
    
    def _process_nto_data(self, weights: np.ndarray, nto_coeff: np.ndarray, state_number: int) -> List[Dict[str, Any]]:
        """Process NTO weights and coefficients to extract orbital pair information."""
        logger.info(f"Processing NTO data for state S{state_number}")
        nto_pairs = []
        
        # Validate input data
        if weights is None or nto_coeff is None:
            logger.error(f"Invalid NTO data: weights={weights}, nto_coeff={nto_coeff}")
            return nto_pairs
            
        if not isinstance(weights, np.ndarray) or not isinstance(nto_coeff, np.ndarray):
            logger.error(f"NTO data not numpy arrays: weights type={type(weights)}, nto_coeff type={type(nto_coeff)}")
            return nto_pairs
            
        logger.debug(f"Weights array: shape={weights.shape}, min={np.min(weights)}, max={np.max(weights)}")
        logger.debug(f"NTO coeff array: shape={nto_coeff.shape}")
        
        # Find HOMO and LUMO indices for reference
        if self.mf.mo_occ is None:
            logger.warning("Orbital occupations not available for NTO analysis")
            return nto_pairs
        
        occupied_indices = np.where(self.mf.mo_occ > 0)[0]
        virtual_indices = np.where(self.mf.mo_occ == 0)[0]
        
        if len(occupied_indices) == 0 or len(virtual_indices) == 0:
            logger.warning("Cannot determine HOMO/LUMO for NTO analysis")
            return nto_pairs
        
        homo_idx = occupied_indices[-1]
        lumo_idx = virtual_indices[0]
        logger.debug(f"Reference orbitals: HOMO idx={homo_idx}, LUMO idx={lumo_idx}")
        
        # NTO weights interpretation based on SVD theory:
        # In SVD of transition density matrix T = U * diag(weights) * V†
        # weights are singular values representing the amplitude of each hole-particle pair
        # Squares of weights give the relative contribution of each pair
        
        # Calculate total contribution for normalization
        # The weights should already be normalized by PySCF, but we verify
        total_contribution = np.sum(weights**2)
        logger.debug(f"Total NTO contribution (should be ~1.0): {total_contribution:.6f}")
        
        # Sort weights in descending order to get most important transitions first
        sorted_indices = np.argsort(weights)[::-1]
        significant_weights = weights[sorted_indices]
        logger.debug(f"Top 5 NTO weights: {significant_weights[:5]}")
        
        # Process NTO pairs starting from most significant
        # Typically only 1-3 pairs contribute significantly (>5% each)
        max_pairs = min(len(weights), 10)  # Reasonable upper limit
        logger.info(f"Processing up to {max_pairs} NTO pairs out of {len(weights)} total")
        
        for i in range(max_pairs):
            idx = sorted_indices[i]
            weight = float(weights[idx])
            
            # Calculate contribution as percentage of total transition
            # Using weight² as this represents the actual contribution amplitude
            contribution_percent = (weight**2 / total_contribution) * 100
            
            # Skip pairs with very small contributions (< 1% of total transition)
            if contribution_percent < 1.0:
                logger.debug(f"Skipping NTO pair {i+1} with contribution {contribution_percent:.2f}%")
                break
            
            # Determine hole and particle orbital descriptions
            hole_desc, hole_idx_desc = self._get_orbital_description(idx, homo_idx, True)
            particle_desc, particle_idx_desc = self._get_orbital_description(idx, lumo_idx, False)
            
            nto_pair = {
                'hole_orbital': hole_desc,
                'particle_orbital': particle_desc,
                'weight': weight,
                'contribution': contribution_percent,  # Now properly calculated as percentage
                'hole_orbital_index': hole_idx_desc,
                'particle_orbital_index': particle_idx_desc
            }
            
            nto_pairs.append(nto_pair)
        
        return nto_pairs
    
    def _get_orbital_description(self, nto_index: int, reference_idx: int, is_hole: bool) -> Tuple[str, int]:
        """Generate human-readable orbital descriptions for NTOs."""
        # This is a simplified version - in reality, the mapping between
        # NTO indices and canonical orbital indices is more complex
        
        if is_hole:
            # For hole orbitals (occupied)
            if nto_index == 0:
                return "HOMO", reference_idx
            elif nto_index == 1:
                return "HOMO-1", reference_idx - 1
            elif nto_index == 2:
                return "HOMO-2", reference_idx - 2
            else:
                return f"HOMO-{nto_index}", reference_idx - nto_index
        else:
            # For particle orbitals (virtual)
            if nto_index == 0:
                return "LUMO", reference_idx
            elif nto_index == 1:
                return "LUMO+1", reference_idx + 1
            elif nto_index == 2:
                return "LUMO+2", reference_idx + 2
            else:
                return f"LUMO+{nto_index}", reference_idx + nto_index