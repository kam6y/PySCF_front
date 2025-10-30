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
from .config_manager import get_memory_for_method

logger = logging.getLogger(__name__)


class TDDFTCalculator(BaseCalculator):
    """TDDFT calculator using PySCF for excited state calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None, optimize_geometry: bool = False):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir, optimize_geometry)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[dft.RKS] = None  # Can be either RKS or UKS
        self.mytd: Optional[tddft.TDDFT] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup TDDFT calculation using the base template method."""
        # Call the base template method which handles common setup
        super().setup_calculation(atoms, **kwargs)
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """Validate TDDFT-specific parameters."""
        # DFT-related parameters
        xc = kwargs.get('xc', 'B3LYP')  # Exchange-correlation functional
        
        # TDDFT-specific parameters
        nstates = kwargs.get('nstates', 10)
        tddft_method = kwargs.get('tddft_method', 'TDDFT')
        analyze_nto = kwargs.get('analyze_nto', False)
        
        # Store parameters for template method access
        self.xc_functional = xc
        self.tddft_nstates = nstates
        self.tddft_method = tddft_method
        self.analyze_nto = analyze_nto
        
        return {
            'xc_functional': xc,
            'method': 'UKS' if kwargs.get('spin', 0) > 0 else 'RKS',
            'tddft_nstates': nstates,
            'tddft_method': tddft_method,
            'tddft_analyze_nto': analyze_nto
        }
    
    def _get_default_memory_mb(self) -> int:
        """Get default memory setting for TDDFT calculations from config."""
        return get_memory_for_method('TDDFT')
    
    def _get_calculation_method_name(self) -> str:
        """Get the name of the calculation method for logging."""
        return 'TDDFT'
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform TDDFT-specific calculation after ground state DFT."""
        # 安全にベースDFTエネルギーを変換
        if base_energy is None:
            raise CalculationError("SCF energy is None - calculation may have failed")
        try:
            scf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert SCF energy to float: {e}")
        
        logger.info(f"Ground state DFT completed. SCF energy: {scf_energy_float:.6f} hartree")
        
        # TDDFT calculation
        logger.info("Starting TDDFT calculation...")
        nstates = getattr(self, 'tddft_nstates', 10)
        tddft_method = getattr(self, 'tddft_method', 'TDDFT')
        
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
        n_occupied = int(np.sum(self.mf.mo_occ > 0)) if hasattr(self.mf, 'mo_occ') else 0
        n_virtual = n_orb - n_occupied
        
        # Calculate a reasonable maximum number of excited states
        # Use the smaller of: half total orbitals or occupied * virtual / 4
        max_reasonable_states = min(n_orb // 2, max(1, (n_occupied * n_virtual) // 4))
        
        if n_orb < 4:
            # Very small systems - limit to 1-2 states
            max_reasonable_states = min(1, max_reasonable_states)
            logger.warning(f"Very small molecular system detected ({n_orb} orbitals, {n_occupied} occupied). "
                           f"TDDFT may not be suitable for such small systems.")
        
        if nstates > max_reasonable_states:
            original_nstates = nstates
            nstates = max(1, max_reasonable_states)
            logger.warning(f"Requested {original_nstates} excited states, but system has only {n_orb} orbitals "
                           f"({n_occupied} occupied, {n_virtual} virtual). Automatically reducing to {nstates} states "
                           f"to avoid convergence issues.")
            # Update the stored parameter
            self.results['tddft_nstates'] = nstates
            self.mytd.nstates = nstates
        
        logger.info(f"TDDFT calculation setup: {nstates} excited states for system with "
                   f"{n_orb} total orbitals ({n_occupied} occupied, {n_virtual} virtual)")
        
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
        
        # Analyze TDDFT results
        tddft_results = self._analyze_tddft_results()
        
        # Combine results
        results = {'scf_energy': scf_energy_float}
        results.update(tddft_results)
        
        return results
    
    def _create_scf_method(self, mol):
        """Create DFT method object for TDDFT ground state (RKS/UKS)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = dft.RKS(mol)
            logger.info("Using Restricted Kohn-Sham (RKS) for closed-shell TDDFT ground state")
        else:
            mf = dft.UKS(mol)
            logger.info("Using Unrestricted Kohn-Sham (UKS) for open-shell TDDFT ground state")
        
        # Set XC functional
        mf.xc = self.xc_functional
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to DFT method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UKS' if spin > 0 else 'RKS'} (TDDFT ground state)"
    
    def _requires_geometry_optimization(self) -> bool:
        """TDDFT does not require geometry optimization."""
        return False
    
    def _requires_frequency_analysis(self) -> bool:
        """TDDFT does not require frequency analysis."""
        return False

    # ===== 修正済みの完全な関数 =====
    def _analyze_tddft_results(self) -> Dict[str, Any]:
        """
        Analyze TDDFT results and extract key spectroscopic information.
        """
        if self.mytd is None or not hasattr(self.mytd, 'e'):
            raise CalculationError("TDDFT results not available")
        
        # Convert excitation energies from hartree to eV
        ev_to_hartree = 27.2114
        
        # 堅牢なエネルギーフィルタリング - None、NaN、無限大、非数値をチェック
        logger.debug(f"Raw excitation energies type: {type(self.mytd.e)}")
        logger.debug(f"Raw excitation energies shape: {getattr(self.mytd.e, 'shape', 'no shape')}")
        logger.debug(f"Raw excitation energies contents: {self.mytd.e}")
        
        valid_energies_hartree = []
        invalid_count = 0
        
        # 安全にイテレートしてフィルタリング
        try:
            for i, energy in enumerate(self.mytd.e):
                # 詳細な値チェック
                if energy is None:
                    logger.warning(f"Excitation energy {i} is None")
                    invalid_count += 1
                    continue
                
                # numpy値の場合の処理
                if hasattr(energy, 'item'):
                    try:
                        energy_val = float(energy.item())
                    except (ValueError, TypeError) as e:
                        logger.warning(f"Excitation energy {i} conversion failed: {e}")
                        invalid_count += 1
                        continue
                else:
                    try:
                        energy_val = float(energy)
                    except (ValueError, TypeError) as e:
                        logger.warning(f"Excitation energy {i} conversion failed: {e}")
                        invalid_count += 1
                        continue
                
                # NaN、無限大チェック
                if not np.isfinite(energy_val):
                    logger.warning(f"Excitation energy {i} is not finite: {energy_val}")
                    invalid_count += 1
                    continue
                
                # 物理的に意味のある値かチェック（励起エネルギーは正値）
                if energy_val <= 0:
                    logger.warning(f"Excitation energy {i} is non-positive: {energy_val}")
                    invalid_count += 1
                    continue
                
                valid_energies_hartree.append(energy_val)
                
        except Exception as filter_error:
            logger.error(f"Error during energy filtering: {filter_error}")
            raise CalculationError(f"Failed to filter excitation energies: {str(filter_error)}")
        
        logger.info(f"Energy filtering results: {len(valid_energies_hartree)} valid, {invalid_count} invalid")
        
        if not valid_energies_hartree:
            # 有効なエネルギーが一つもない場合はエラーとする
            raise CalculationError(f"TDDFT calculation failed: no valid excitation energies obtained. {invalid_count} invalid values found.")

        # 安全にeVに変換
        try:
            excitation_energies_ev = [float(e * ev_to_hartree) for e in valid_energies_hartree]
            logger.debug(f"Converted energies to eV: {excitation_energies_ev[:5]}...")
        except Exception as conversion_error:
            logger.error(f"Error converting energies to eV: {conversion_error}")
            raise CalculationError(f"Failed to convert excitation energies to eV: {str(conversion_error)}")
        
        # Convert to wavelengths in nm
        # E(eV) = hc/λ, where hc ≈ 1239.84 eV·nm
        excitation_wavelengths = [1239.84 / e if e > 0 else 0 for e in excitation_energies_ev]
        
        # Get oscillator strengths
        oscillator_strengths = []
        transition_dipoles = []
        
        # Calculate oscillator strengths using length gauge (standard approach)
        try:
            # 計算された状態の数だけスライスする
            osc_strengths_result = self.mytd.oscillator_strength(gauge='length')[:len(valid_energies_hartree)]
            if osc_strengths_result is not None:
                # 安全にfloatに変換
                for i, f in enumerate(osc_strengths_result):
                    try:
                        if f is not None and np.isfinite(f):
                            oscillator_strengths.append(float(f))
                        else:
                            logger.warning(f"Invalid oscillator strength at index {i}: {f}")
                            oscillator_strengths.append(0.0)
                    except (ValueError, TypeError) as e:
                        logger.warning(f"Failed to convert oscillator strength {i}: {e}")
                        oscillator_strengths.append(0.0)
                logger.info(f"Retrieved {len(oscillator_strengths)} oscillator strengths (length gauge)")
        except Exception as osc_error:
            logger.warning(f"Failed to get oscillator strengths: {osc_error}")
        
        # Calculate transition dipole moments
        try:
            # 計算された状態の数だけスライスする
            dipole_result = self.mytd.transition_dipole()[:len(valid_energies_hartree)]
            if dipole_result is not None:
                for i, dipole in enumerate(dipole_result):
                    try:
                        if hasattr(dipole, '__len__') and len(dipole) >= 3:
                            # 各成分を安全に変換
                            x = float(dipole[0]) if dipole[0] is not None and np.isfinite(dipole[0]) else 0.0
                            y = float(dipole[1]) if dipole[1] is not None and np.isfinite(dipole[1]) else 0.0
                            z = float(dipole[2]) if dipole[2] is not None and np.isfinite(dipole[2]) else 0.0
                            
                            transition_dipoles.append({
                                'x': x,
                                'y': y,
                                'z': z
                            })
                        else:
                            logger.warning(f"Invalid dipole moment at index {i}: {dipole}")
                            transition_dipoles.append({'x': 0.0, 'y': 0.0, 'z': 0.0})
                    except (ValueError, TypeError, IndexError) as e:
                        logger.warning(f"Failed to process dipole moment {i}: {e}")
                        transition_dipoles.append({'x': 0.0, 'y': 0.0, 'z': 0.0})
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
    
    
    def _analyze_natural_transition_orbitals(self, excitation_energies: List[float]) -> List[Dict[str, Any]]:
        """
        Perform Natural Transition Orbital (NTO) analysis for excited states.
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
        
        # Calculate total contribution for normalization
        total_contribution = np.sum(weights**2)
        logger.debug(f"Total NTO contribution (should be ~1.0): {total_contribution:.6f}")
        
        # Sort weights in descending order to get most important transitions first
        sorted_indices = np.argsort(weights)[::-1]
        significant_weights = weights[sorted_indices]
        logger.debug(f"Top 5 NTO weights: {significant_weights[:5]}")
        
        # Process NTO pairs starting from most significant
        max_pairs = min(len(weights), 10)  # Reasonable upper limit
        logger.info(f"Processing up to {max_pairs} NTO pairs out of {len(weights)} total")
        
        for i in range(max_pairs):
            idx = sorted_indices[i]
            weight = float(weights[idx])
            
            # Calculate contribution as percentage of total transition
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
                'contribution': contribution_percent,
                'hole_orbital_index': hole_idx_desc,
                'particle_orbital_index': particle_idx_desc
            }
            
            nto_pairs.append(nto_pair)
        
        return nto_pairs
    
    def _get_orbital_description(self, nto_index: int, reference_idx: int, is_hole: bool) -> Tuple[str, int]:
        """Generate human-readable orbital descriptions for NTOs."""
        
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