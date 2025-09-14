"""CASSCF calculator implementation using PySCF."""

import os
import logging
import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from pyscf import gto, scf, mcscf

from .base_calculator import BaseCalculator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError
from .file_manager import CalculationFileManager
from .solvent_effects import setup_solvent_effects

logger = logging.getLogger(__name__)


class CASSCFCalculator(BaseCalculator):
    """CASSCF calculator using PySCF for multi-configurational self-consistent field calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None  # Can be RHF or UHF
        self.mycas: Optional[mcscf.casscf.CASSCF] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        # Set calculation method for memory management
        self.calculation_method = 'CASSCF'
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup CASSCF calculation with molecular geometry and parameters."""
        logger.info("Starting CASSCF calculation setup...")
        logger.info(f"Received parameters: {list(kwargs.keys())}")
        try:
            # Extract basic calculation parameters
            basis = kwargs.get('basis', '6-31G(d)')
            charge = kwargs.get('charge', 0)
            spin = kwargs.get('spin', 0)
            max_cycle = kwargs.get('max_cycle', 150)
            solvent_method = kwargs.get('solvent_method', 'none')
            solvent = kwargs.get('solvent', '-')
            memory_mb = kwargs.get('memory_mb', 6000)  # Default 6GB for CASSCF
            
            # CASSCF-specific parameters
            ncas = kwargs.get('ncas', 4)  # Number of active space orbitals
            nelecas = kwargs.get('nelecas', 4)  # Number of active space electrons
            max_cycle_macro = kwargs.get('max_cycle_macro', 50)  # CASSCF macro iterations
            max_cycle_micro = kwargs.get('max_cycle_micro', 3)  # CI solver micro iterations
            analyze_nto = kwargs.get('analyze_nto', False)  # Natural transition orbitals
            natorb = kwargs.get('natorb', True)  # Transform to natural orbitals
            conv_tol = kwargs.get('conv_tol', 1e-6)  # Convergence tolerance
            conv_tol_grad = kwargs.get('conv_tol_grad', 1e-4)  # Gradient tolerance
            
            # AH solver parameters for improved convergence
            ah_conv_tol = kwargs.get('ah_conv_tol', 1e-12)  # AH solver convergence tolerance
            ah_max_cycle = kwargs.get('ah_max_cycle', 30)  # AH solver max iterations
            ah_lindep = kwargs.get('ah_lindep', 1e-14)  # AH linear dependence threshold
            
            # Validate and adjust CASSCF parameters with improved defaults
            if ncas <= 0:
                logger.warning(f"Invalid ncas={ncas}, adjusting to default ncas=4")
                ncas = 4
            if nelecas <= 0:
                logger.warning(f"Invalid nelecas={nelecas}, adjusting to default nelecas=4")
                nelecas = 4
            if nelecas > 2 * ncas:
                logger.warning(f"Too many electrons ({nelecas}) for active space size ({ncas} orbitals)")
                # Adjust nelecas to maximum possible for the given ncas
                nelecas = 2 * ncas
                logger.info(f"Adjusted nelecas to maximum possible: {nelecas}")
            if max_cycle_macro <= 0:
                logger.warning(f"Invalid max_cycle_macro={max_cycle_macro}, adjusting to default 50")
                max_cycle_macro = 50

            # Additional sanity checks with warnings instead of errors
            if ncas > 20:
                logger.warning(f"Large active space (ncas={ncas}) may require substantial computational resources")
            if nelecas > 20:
                logger.warning(f"Many active electrons (nelecas={nelecas}) may require substantial computational resources")
            if max_cycle_macro > 100:
                logger.warning(f"Many CASSCF macro cycles ({max_cycle_macro}) may require long computation time")
            
            # Convert atoms list to PySCF format
            atom_string = self._atoms_to_string(atoms)
            
            # Create molecular object
            logger.info(f"Creating PySCF molecular object with {len(atoms)} atoms, basis={basis}, charge={charge}, spin={spin}")
            self.mol = gto.M(
                atom=atom_string,
                basis=basis,
                charge=charge,
                spin=spin,
                verbose=0
            )
            logger.info("PySCF molecular object created successfully")
            
            # Apply memory settings
            if memory_mb and memory_mb > 0:
                self.mol.max_memory = memory_mb
            else:
                self.mol.max_memory = 6000  # CASSCF needs substantial memory (6GB default)
            
            # Setup SCF calculation (prerequisite for CASSCF)
            # For closed-shell systems (spin=0), use RHF
            # For open-shell systems (spin>0), use UHF
            if spin == 0:
                self.mf = scf.RHF(self.mol)
                logger.info("Using Restricted Hartree-Fock (RHF) for closed-shell CASSCF reference")
            else:
                self.mf = scf.UHF(self.mol)
                logger.info("Using Unrestricted Hartree-Fock (UHF) for open-shell CASSCF reference")
            
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.max_cycle = max_cycle
            
            # Store parameters for template method (use adjusted values)
            self.max_cycle = max_cycle
            self.solvent_method = solvent_method
            self.solvent = solvent
            self.ncas = ncas  # This now contains the adjusted value
            self.nelecas = nelecas  # This now contains the adjusted value
            self.max_cycle_macro = max_cycle_macro  # This now contains the adjusted value
            self.max_cycle_micro = max_cycle_micro
            self.analyze_nto = analyze_nto
            self.natorb = natorb
            self.conv_tol = conv_tol
            self.conv_tol_grad = conv_tol_grad
            
            # Store parameters in results for reference (use adjusted values)
            self.results.update({
                'basis': basis,
                'charge': charge,
                'spin': spin,
                'max_cycle': max_cycle,
                'solvent_method': solvent_method,
                'solvent': solvent,
                'atom_count': len(atoms),
                'ncas': ncas,  # Adjusted value
                'nelecas': nelecas,  # Adjusted value
                'max_cycle_macro': max_cycle_macro,  # Adjusted value
                'max_cycle_micro': max_cycle_micro,
                'analyze_nto': analyze_nto,
                'natorb': natorb,
                'conv_tol': conv_tol,
                'conv_tol_grad': conv_tol_grad,
                'ah_conv_tol': ah_conv_tol,
                'ah_max_cycle': ah_max_cycle,
                'ah_lindep': ah_lindep,
                'method': 'UHF-CASSCF' if spin > 0 else 'RHF-CASSCF'
            })
            
            logger.info(f"CASSCF setup completed: {ncas} orbitals, {nelecas} electrons in active space")
            
        except Exception as e:
            raise InputError(f"Failed to setup CASSCF calculation: {str(e)}")
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform CASSCF-specific calculation after SCF."""
        # Safely convert base SCF energy
        if base_energy is None:
            raise CalculationError("SCF energy is None - calculation may have failed")
        try:
            scf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert SCF energy to float: {e}")
        
        logger.info(f"Reference SCF completed. SCF energy: {scf_energy_float:.6f} hartree")
        
        # CASSCF calculation
        logger.info("Starting CASSCF calculation...")
        ncas = getattr(self, 'ncas', 4)
        nelecas = getattr(self, 'nelecas', 4)
        natorb = getattr(self, 'natorb', True)
        max_cycle_macro = getattr(self, 'max_cycle_macro', 50)
        max_cycle_micro = getattr(self, 'max_cycle_micro', 3)
        conv_tol = getattr(self, 'conv_tol', 1e-6)
        conv_tol_grad = getattr(self, 'conv_tol_grad', 1e-4)
        ah_conv_tol = getattr(self, 'ah_conv_tol', 1e-12)
        ah_max_cycle = getattr(self, 'ah_max_cycle', 30)
        ah_lindep = getattr(self, 'ah_lindep', 1e-14)
        
        # Create CASSCF object
        self.mycas = mcscf.CASSCF(self.mf, ncas, nelecas)
        
        # Set CASSCF parameters
        self.mycas.natorb = natorb
        self.mycas.max_cycle_macro = max_cycle_macro
        self.mycas.max_cycle_micro = max_cycle_micro
        self.mycas.conv_tol = conv_tol
        self.mycas.conv_tol_grad = conv_tol_grad
        
        # Set AH solver parameters for improved convergence
        self.mycas.ah_conv_tol = ah_conv_tol
        self.mycas.ah_max_cycle = ah_max_cycle
        self.mycas.ah_lindep = ah_lindep
        
        # Validate active space parameters against molecular orbital structure
        if hasattr(self.mf, 'mo_energy') and self.mf.mo_energy is not None:
            n_orb = len(self.mf.mo_energy)
            n_occupied = int(np.sum(self.mf.mo_occ > 0)) if hasattr(self.mf, 'mo_occ') else 0
            
            if ncas > n_orb:
                raise InputError(f"Active space size ({ncas}) exceeds total orbitals ({n_orb})")
            if ncas > n_occupied + 10:  # Reasonable check for too many virtual orbitals
                logger.warning(f"Large active space ({ncas} orbitals) includes many virtual orbitals. "
                             f"System has {n_occupied} occupied orbitals out of {n_orb} total.")
        
        logger.info(f"CASSCF setup: {ncas} active orbitals, {nelecas} active electrons")
        logger.info(f"Convergence criteria: energy = {conv_tol:.2e}, gradient = {conv_tol_grad:.2e}")
        logger.info(f"Maximum cycles: macro = {max_cycle_macro}, micro = {max_cycle_micro}")
        if natorb:
            logger.info("Natural orbital transformation enabled")
        
        # Run CASSCF calculation
        logger.info("Starting CASSCF kernel calculation...")
        try:
            casscf_result = self.mycas.kernel()
        except Exception as e:
            import traceback
            logger.error(f"CASSCF kernel() failed with exception type: {type(e).__name__}")
            logger.error(f"CASSCF kernel() exception message: '{str(e)}'")
            logger.error(f"CASSCF kernel() full traceback:\n{traceback.format_exc()}")
            
            error_msg = str(e)
            error_msg_lower = error_msg.lower()
            
            # Handle empty error messages
            if not error_msg.strip():
                error_msg = f"CASSCF kernel() failed with {type(e).__name__} (empty error message)"
                logger.error(f"Empty error message detected in CASSCF kernel(), using: {error_msg}")
            
            if "singular" in error_msg_lower or "convergence" in error_msg_lower:
                raise ConvergenceError(f"CASSCF calculation failed to converge: {error_msg}")
            elif "memory" in error_msg_lower:
                raise CalculationError(f"CASSCF calculation failed due to insufficient memory: {error_msg}")
            elif "maximum" in error_msg_lower and "cycle" in error_msg_lower:
                raise ConvergenceError(f"CASSCF reached maximum cycles without convergence: {error_msg}")
            elif isinstance(e, AssertionError):
                # AssertionError in PySCF often indicates active space parameter or convergence issues
                detailed_msg = f"CASSCF calculation failed with AssertionError: {error_msg}"
                suggestions = []
                suggestions.append(f"Current active space: ncas={ncas}, nelecas={nelecas}")

                # Check if parameters are reasonable
                if ncas > 8:
                    suggestions.append("Large active space may be too demanding - try reducing ncas")
                if nelecas > ncas * 2:
                    suggestions.append("Too many electrons for active space size - check nelecas")
                if max_cycle_macro < 30:
                    suggestions.append("Low macro cycle limit may prevent convergence - try increasing max_cycle_macro")
                if hasattr(self.mf, 'mo_energy') and self.mf.mo_energy is not None:
                    n_orb = len(self.mf.mo_energy)
                    if ncas > n_orb // 2:
                        suggestions.append("Active space too large relative to basis set - reduce ncas")

                if suggestions:
                    detailed_msg += f". Suggestions: {'; '.join(suggestions)}"

                logger.error(f"AssertionError diagnosis: {detailed_msg}")
                raise CalculationError(detailed_msg)
            else:
                raise CalculationError(f"CASSCF calculation failed: {error_msg}")
        
        # Validate results
        if casscf_result is None:
            raise CalculationError("CASSCF calculation failed: no result obtained")
        
        # Process kernel() return value comprehensively
        casscf_energy = None
        ci_coefficients = None
        kernel_additional_info = {}
        
        if isinstance(casscf_result, tuple):
            logger.info(f"CASSCF kernel returned tuple with {len(casscf_result)} elements")
            
            # Extract energy (always first element)
            if len(casscf_result) >= 1:
                casscf_energy = casscf_result[0]
                logger.info("Extracted energy from CASSCF kernel result tuple")
            
            # Extract CI coefficients (typically second element)
            if len(casscf_result) >= 2:
                ci_coefficients = casscf_result[1]
                logger.info("Extracted CI coefficients from CASSCF kernel result")
                kernel_additional_info['ci_coefficients_shape'] = str(getattr(ci_coefficients, 'shape', 'no_shape'))
            
            # Store additional information for analysis
            kernel_additional_info['tuple_length'] = len(casscf_result)
            for i, item in enumerate(casscf_result[2:], start=2):
                try:
                    if hasattr(item, 'shape'):
                        kernel_additional_info[f'element_{i}_shape'] = str(item.shape)
                    elif isinstance(item, (int, float, complex)):
                        kernel_additional_info[f'element_{i}_value'] = float(item)
                    else:
                        kernel_additional_info[f'element_{i}_type'] = type(item).__name__
                except Exception:
                    kernel_additional_info[f'element_{i}_type'] = 'unknown'
        else:
            # Single value returned
            casscf_energy = casscf_result
            logger.info("CASSCF kernel returned single value (energy)")
        
        # Validate energy
        if casscf_energy is None:
            raise CalculationError("CASSCF calculation failed: no energy obtained")
        
        try:
            casscf_energy_float = float(casscf_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert CASSCF energy to float: {e}")
        
        # Store CI coefficients for analysis
        if ci_coefficients is not None:
            self.ci_coefficients = ci_coefficients
        elif hasattr(self.mycas, 'ci') and self.mycas.ci is not None:
            self.ci_coefficients = self.mycas.ci
            logger.info("Using CI coefficients from CASSCF object")
        else:
            self.ci_coefficients = None
        
        logger.info(f"CASSCF kernel result processing complete. Energy: {casscf_energy_float:.6f} hartree")
        
        # Check convergence status
        converged = getattr(self.mycas, 'converged', False)
        if not converged:
            logger.warning("CASSCF calculation did not converge to specified tolerance")
        else:
            logger.info("CASSCF calculation converged successfully")
        
        # Get iteration information
        e_tot_history = getattr(self.mycas, 'e_tot', None)
        if hasattr(e_tot_history, '__len__') and len(e_tot_history) > 1:
            iterations_performed = len(e_tot_history) - 1
        else:
            iterations_performed = getattr(self.mycas, 'itercount', 0)
        
        logger.info(f"CASSCF calculation completed. Energy: {casscf_energy_float:.6f} hartree")
        logger.info(f"Performed {iterations_performed} macro iterations")
        
        # Analyze CASSCF results (including CI coefficients)
        casscf_results = self._analyze_casscf_results()
        
        # Enhanced CI coefficients analysis if available
        enhanced_ci_analysis = {}
        if self.ci_coefficients is not None:
            enhanced_ci_analysis = self._analyze_kernel_ci_coefficients(self.ci_coefficients)
        
        # Combine results with comprehensive information
        results = {
            'scf_energy': scf_energy_float,
            'casscf_energy': casscf_energy_float,
            'correlation_energy': casscf_energy_float - scf_energy_float,
            'converged': converged,
            'macro_iterations': iterations_performed,
            'kernel_return_info': kernel_additional_info,
            'ci_coefficients_available': self.ci_coefficients is not None,
            'enhanced_ci_analysis': enhanced_ci_analysis
        }
        results.update(casscf_results)
        
        logger.info(f"CASSCF calculation completed with enhanced analysis. Energy: {casscf_energy_float:.6f} hartree, Converged: {converged}")
        return results
    
    def _create_scf_method(self, mol):
        """Create SCF method object for CASSCF reference (RHF/UHF)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = scf.RHF(mol)
            logger.info("Using Restricted Hartree-Fock (RHF) for closed-shell CASSCF reference")
        else:
            mf = scf.UHF(mol)
            logger.info("Using Unrestricted Hartree-Fock (UHF) for open-shell CASSCF reference")
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to SCF method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UHF' if spin > 0 else 'RHF'} (CASSCF reference)"
    
    def _requires_geometry_optimization(self) -> bool:
        """CASSCF does not require geometry optimization by default."""
        return False
    
    def _requires_frequency_analysis(self) -> bool:
        """CASSCF does not require frequency analysis by default."""
        return False
    
    # ===== CASSCF-specific Analysis Methods =====
    
    def _analyze_casscf_results(self) -> Dict[str, Any]:
        """
        Perform comprehensive analysis of CASSCF results.
        """
        if self.mycas is None:
            raise CalculationError("CASSCF calculation not available for analysis")
        
        logger.info("Starting CASSCF results analysis...")
        
        # Set verbose level for detailed analysis
        original_verbose = self.mycas.verbose
        self.mycas.verbose = 4
        
        # Perform PySCF's built-in analysis
        logger.info("Running PySCF CASSCF analysis...")
        try:
            self.mycas.analyze()
        except Exception as e:
            logger.warning(f"PySCF CASSCF analysis encountered issues: {e}")
        
        # Restore original verbose level
        self.mycas.verbose = original_verbose
        
        # Extract analysis results
        analysis_results = {}
        
        # Natural orbital analysis
        try:
            natural_orbital_analysis = self._analyze_natural_orbitals()
            analysis_results['natural_orbital_analysis'] = natural_orbital_analysis
            logger.info("Natural orbital analysis completed successfully")
        except Exception as e:
            logger.error(f"Natural orbital analysis failed: {e}")
            analysis_results['natural_orbital_analysis'] = None
        
        # CI coefficient analysis
        try:
            ci_analysis = self._analyze_ci_coefficients()
            analysis_results['ci_coefficient_analysis'] = ci_analysis
            logger.info("CI coefficient analysis completed successfully")
        except Exception as e:
            logger.error(f"CI coefficient analysis failed: {e}")
            analysis_results['ci_coefficient_analysis'] = None
        
        # Orbital rotation analysis (CASSCF specific)
        try:
            rotation_analysis = self._analyze_orbital_rotation()
            analysis_results['orbital_rotation_analysis'] = rotation_analysis
            logger.info("Orbital rotation analysis completed successfully")
        except Exception as e:
            logger.error(f"Orbital rotation analysis failed: {e}")
            analysis_results['orbital_rotation_analysis'] = None
        
        # Mulliken spin density analysis (for open-shell systems)
        if self.results.get('spin', 0) > 0:
            try:
                spin_analysis = self._calculate_mulliken_spin_density()
                analysis_results['mulliken_spin_analysis'] = spin_analysis
                logger.info("Mulliken spin density analysis completed successfully")
            except Exception as e:
                logger.error(f"Mulliken spin density analysis failed: {e}")
                analysis_results['mulliken_spin_analysis'] = None
        
        # Orbital overlap analysis
        try:
            overlap_analysis = self._calculate_orbital_overlap()
            analysis_results['orbital_overlap_analysis'] = overlap_analysis
            logger.info("Orbital overlap analysis completed successfully")
        except Exception as e:
            logger.error(f"Orbital overlap analysis failed: {e}")
            analysis_results['orbital_overlap_analysis'] = None
        
        logger.info("CASSCF results analysis completed")
        return analysis_results
    
    def _analyze_natural_orbitals(self) -> Dict[str, Any]:
        """
        Analyze natural orbitals and their occupation numbers.
        (Identical to CASCI implementation but may have different orbital coefficients due to optimization)
        """
        if not hasattr(self.mycas, 'natorb') or not self.mycas.natorb:
            logger.warning("Natural orbital analysis requested but natorb is not enabled")
            return {'enabled': False, 'reason': 'Natural orbital transformation not enabled'}
        
        natural_orbital_data = {'enabled': True}
        
        # Extract natural orbital occupation numbers
        if hasattr(self.mycas, 'ci') and self.mycas.ci is not None:
            # Get natural orbitals and occupation numbers
            try:
                # Use PySCF's cas_natorb function to get natural orbitals
                natorbs, ci_coeffs, occs = self.mycas.cas_natorb()
                
                # Process occupation numbers
                if occs is not None and len(occs) > 0:
                    # Convert to Python floats for JSON serialization
                    natural_occupations = [float(occ) for occ in occs]
                    
                    # Analyze occupation pattern
                    strongly_occupied = [occ for occ in natural_occupations if occ > 1.5]
                    weakly_occupied = [occ for occ in natural_occupations if 0.1 < occ <= 1.5]
                    virtual = [occ for occ in natural_occupations if occ <= 0.1]
                    
                    natural_orbital_data.update({
                        'occupation_numbers': natural_occupations,
                        'strongly_occupied_count': len(strongly_occupied),
                        'weakly_occupied_count': len(weakly_occupied),
                        'virtual_count': len(virtual),
                        'total_orbitals': len(natural_occupations)
                    })
                    
                    # Calculate effective number of electron pairs and unpaired electrons
                    total_electrons = sum(natural_occupations)
                    effective_pairs = sum([min(occ, 2.0) for occ in natural_occupations]) / 2.0
                    unpaired_electrons = total_electrons - 2.0 * effective_pairs
                    
                    natural_orbital_data.update({
                        'total_active_electrons': float(total_electrons),
                        'effective_electron_pairs': float(effective_pairs),
                        'effective_unpaired_electrons': float(unpaired_electrons)
                    })
                    
                    logger.info(f"Natural orbital analysis: {len(natural_occupations)} orbitals, "
                               f"{total_electrons:.2f} electrons, "
                               f"{len(strongly_occupied)} strongly occupied, "
                               f"{len(weakly_occupied)} weakly occupied")
                else:
                    natural_orbital_data['occupation_numbers'] = None
                    logger.warning("No natural orbital occupation numbers found")
                    
            except Exception as e:
                logger.error(f"Error extracting natural orbital data: {e}")
                natural_orbital_data['error'] = str(e)
        else:
            natural_orbital_data['occupation_numbers'] = None
            logger.warning("No CI coefficients available for natural orbital analysis")
        
        return natural_orbital_data
    
    def _analyze_ci_coefficients(self) -> Dict[str, Any]:
        """
        Analyze CI coefficients to identify major configurations.
        (Identical to CASCI implementation)
        """
        ci_analysis = {}
        
        if not hasattr(self.mycas, 'ci') or self.mycas.ci is None:
            logger.warning("No CI coefficients available for analysis")
            return {'available': False, 'reason': 'No CI coefficients found'}
        
        ci_analysis['available'] = True
        
        try:
            # For CASSCF, ci is typically a numpy array
            ci_coeffs = self.mycas.ci
            
            if hasattr(ci_coeffs, 'flatten'):
                # Multi-dimensional CI vector - flatten it
                ci_flat = ci_coeffs.flatten()
            else:
                ci_flat = np.array(ci_coeffs)
            
            # Find configurations with significant contributions
            ci_squared = ci_flat**2
            total_norm = np.sum(ci_squared)
            
            # Sort by magnitude
            sorted_indices = np.argsort(np.abs(ci_flat))[::-1]
            
            # Extract major configurations (contributions > 1%)
            major_configurations = []
            cumulative_weight = 0.0
            threshold = 0.01  # 1% threshold
            
            for i, idx in enumerate(sorted_indices):
                coeff = float(ci_flat[idx])
                contribution = float(ci_squared[idx] / total_norm)
                cumulative_weight += contribution
                
                if contribution >= threshold and len(major_configurations) < 20:  # Limit to top 20
                    major_configurations.append({
                        'configuration_index': int(idx),
                        'coefficient': coeff,
                        'contribution_percent': contribution * 100.0,
                        'cumulative_percent': cumulative_weight * 100.0
                    })
                
                # Stop when we've captured 98% of the wavefunction
                if cumulative_weight >= 0.98:
                    break
            
            ci_analysis.update({
                'total_configurations': len(ci_flat),
                'major_configurations': major_configurations,
                'leading_coefficient': float(ci_flat[sorted_indices[0]]),
                'leading_contribution_percent': float(ci_squared[sorted_indices[0]] / total_norm * 100),
                'multiconfigurational_character': 100.0 - float(ci_squared[sorted_indices[0]] / total_norm * 100)
            })
            
            logger.info(f"CI analysis: {len(major_configurations)} major configurations identified, "
                       f"leading contribution: {ci_analysis['leading_contribution_percent']:.1f}%")
            
        except Exception as e:
            logger.error(f"Error in CI coefficient analysis: {e}")
            ci_analysis['error'] = str(e)
        
        return ci_analysis
    
    def _analyze_orbital_rotation(self) -> Dict[str, Any]:
        """
        Analyze the orbital rotation that occurred during CASSCF optimization.
        This is specific to CASSCF and not available in CASCI.
        """
        rotation_analysis = {'available': True}
        
        try:
            # Check if we have rotation information
            if hasattr(self.mycas, 'mo_coeff') and hasattr(self.mf, 'mo_coeff'):
                # Calculate orbital rotation matrix
                if hasattr(self.mol, 'get_ovlp'):
                    S = self.mol.get_ovlp()
                else:
                    S = self.mf.get_ovlp()
                
                # Rotation matrix: R = C_initial^T * S * C_final
                rotation_matrix = np.dot(self.mf.mo_coeff.T, np.dot(S, self.mycas.mo_coeff))
                
                # Analyze rotation for different orbital spaces
                ncore = self.mycas.ncore
                ncas = self.results.get('ncas', 4)
                
                # Core orbital rotations
                if ncore > 0:
                    core_rotation = rotation_matrix[:ncore, :ncore]
                    core_rotation_magnitude = np.max(np.abs(core_rotation - np.eye(ncore)))
                else:
                    core_rotation_magnitude = 0.0
                
                # Active orbital rotations
                active_start = ncore
                active_end = ncore + ncas
                active_rotation = rotation_matrix[active_start:active_end, active_start:active_end]
                active_rotation_magnitude = np.max(np.abs(active_rotation - np.eye(ncas)))
                
                # Virtual orbital rotations (if any)
                nvirt = rotation_matrix.shape[0] - active_end
                if nvirt > 0:
                    virtual_rotation = rotation_matrix[active_end:, active_end:]
                    virtual_rotation_magnitude = np.max(np.abs(virtual_rotation - np.eye(nvirt)))
                else:
                    virtual_rotation_magnitude = 0.0
                
                rotation_analysis.update({
                    'core_orbital_rotation_magnitude': float(core_rotation_magnitude),
                    'active_orbital_rotation_magnitude': float(active_rotation_magnitude),
                    'virtual_orbital_rotation_magnitude': float(virtual_rotation_magnitude),
                    'overall_rotation_magnitude': float(np.max([
                        core_rotation_magnitude, 
                        active_rotation_magnitude, 
                        virtual_rotation_magnitude
                    ]))
                })
                
                # Classify rotation extent
                max_rotation = rotation_analysis['overall_rotation_magnitude']
                if max_rotation < 0.1:
                    rotation_extent = "minimal"
                elif max_rotation < 0.3:
                    rotation_extent = "moderate"
                else:
                    rotation_extent = "significant"
                
                rotation_analysis['rotation_extent'] = rotation_extent
                
                logger.info(f"Orbital rotation analysis: core = {core_rotation_magnitude:.3f}, "
                           f"active = {active_rotation_magnitude:.3f}, "
                           f"virtual = {virtual_rotation_magnitude:.3f} ({rotation_extent})")
            else:
                rotation_analysis.update({
                    'available': False,
                    'reason': 'Initial or final orbital coefficients not available'
                })
                
        except Exception as e:
            logger.error(f"Error in orbital rotation analysis: {e}")
            rotation_analysis['error'] = str(e)
        
        return rotation_analysis
    
    def _calculate_mulliken_spin_density(self) -> Dict[str, Any]:
        """
        Calculate Mulliken atomic spin densities for open-shell systems.
        (Identical to CASCI implementation)
        """
        if not hasattr(self.mycas, 'make_rdm1s') or not callable(self.mycas.make_rdm1s):
            logger.warning("CASSCF object does not support spin density calculation")
            return {'available': False, 'reason': 'Spin density not supported for this CASSCF method'}
        
        spin_analysis = {'available': True}
        
        try:
            # Get spin density matrices
            dm1a, dm1b = self.mycas.make_rdm1s()  # Alpha and beta density matrices
            
            if dm1a is None or dm1b is None:
                logger.warning("Could not obtain alpha/beta density matrices")
                return {'available': False, 'reason': 'Density matrices not available'}
            
            # Calculate spin density matrix (alpha - beta)
            spin_dm = dm1a - dm1b
            
            # Perform Mulliken population analysis on spin density
            if hasattr(self.mol, 'get_ovlp'):
                ovlp = self.mol.get_ovlp()
            else:
                ovlp = self.mycas._scf.get_ovlp()
            
            # Mulliken spin populations
            spin_pop = np.einsum('ij,ji->i', spin_dm, ovlp)
            
            # Group by atoms
            atomic_spin_densities = []
            ao_offset = 0
            
            for iatom in range(self.mol.natm):
                atom_symbol = self.mol.atom_symbol(iatom)
                nao = self.mol.aoslice_by_atom()[iatom][3] - self.mol.aoslice_by_atom()[iatom][2]
                
                # Sum spin populations for this atom
                atom_spin = float(np.sum(spin_pop[ao_offset:ao_offset + nao]))
                
                atomic_spin_densities.append({
                    'atom_index': iatom,
                    'element': atom_symbol,
                    'spin_density': atom_spin,
                    'abs_spin_density': abs(atom_spin)
                })
                
                ao_offset += nao
            
            # Calculate total spin
            total_spin_density = sum([atom['spin_density'] for atom in atomic_spin_densities])
            total_abs_spin = sum([atom['abs_spin_density'] for atom in atomic_spin_densities])
            
            spin_analysis.update({
                'atomic_spin_densities': atomic_spin_densities,
                'total_spin_density': float(total_spin_density),
                'total_absolute_spin_density': float(total_abs_spin),
                'expected_spin': float(self.mol.spin)
            })
            
            logger.info(f"Spin density analysis: total = {total_spin_density:.3f}, "
                       f"expected = {self.mol.spin}, atoms analyzed = {len(atomic_spin_densities)}")
            
        except Exception as e:
            logger.error(f"Error in Mulliken spin density calculation: {e}")
            spin_analysis['error'] = str(e)
        
        return spin_analysis
    
    def _calculate_orbital_overlap(self) -> Dict[str, Any]:
        """
        Calculate overlap between CASSCF orbitals and reference SCF orbitals.
        (Identical to CASCI implementation but with optimized orbitals)
        """
        overlap_analysis = {}
        
        if not hasattr(self.mycas, 'mo_coeff') or self.mycas.mo_coeff is None:
            logger.warning("CASSCF orbitals not available for overlap analysis")
            return {'available': False, 'reason': 'CASSCF orbitals not found'}
        
        if not hasattr(self.mf, 'mo_coeff') or self.mf.mo_coeff is None:
            logger.warning("SCF orbitals not available for overlap analysis")
            return {'available': False, 'reason': 'SCF reference orbitals not found'}
        
        overlap_analysis['available'] = True
        
        try:
            # Get overlap matrix
            if hasattr(self.mol, 'get_ovlp'):
                S = self.mol.get_ovlp()
            else:
                S = self.mf.get_ovlp()
            
            casscf_orbs = self.mycas.mo_coeff
            scf_orbs = self.mf.mo_coeff
            
            # Calculate orbital overlap matrix: <CASSCF_i|SCF_j>
            # Overlap = C_CASSCF^T * S * C_SCF
            overlap_matrix = np.dot(casscf_orbs.T, np.dot(S, scf_orbs))
            
            # Analyze overlap for active space orbitals
            ncas = self.results.get('ncas', 4)
            active_start = self.mycas.ncore  # Core orbitals
            active_end = active_start + ncas
            
            # Extract overlaps for active space orbitals
            active_overlaps = overlap_matrix[active_start:active_end, :]
            
            # Find dominant SCF character for each active orbital
            active_orbital_analysis = []
            for i, active_orb_overlaps in enumerate(active_overlaps):
                abs_overlaps = np.abs(active_orb_overlaps)
                max_overlap_idx = np.argmax(abs_overlaps)
                max_overlap_val = float(abs_overlaps[max_overlap_idx])
                
                # Determine orbital type in SCF reference
                if hasattr(self.mf, 'mo_occ'):
                    if max_overlap_idx < len(self.mf.mo_occ):
                        scf_occ = float(self.mf.mo_occ[max_overlap_idx])
                        if scf_occ > 1.5:
                            orbital_type = "occupied"
                        elif scf_occ > 0.1:
                            orbital_type = "partially_occupied" 
                        else:
                            orbital_type = "virtual"
                    else:
                        orbital_type = "unknown"
                else:
                    orbital_type = "unknown"
                
                active_orbital_analysis.append({
                    'active_orbital_index': i,
                    'dominant_scf_orbital': int(max_overlap_idx),
                    'max_overlap': max_overlap_val,
                    'scf_orbital_type': orbital_type
                })
            
            # Calculate overall transformation character
            avg_max_overlap = np.mean([orb['max_overlap'] for orb in active_orbital_analysis])
            
            overlap_analysis.update({
                'active_space_orbitals': ncas,
                'active_orbital_analysis': active_orbital_analysis,
                'average_max_overlap': float(avg_max_overlap),
                'orbital_transformation_character': "minimal" if avg_max_overlap > 0.9 else "significant"
            })
            
            logger.info(f"Orbital overlap analysis: {ncas} active orbitals, "
                       f"average max overlap = {avg_max_overlap:.3f}")
            
        except Exception as e:
            logger.error(f"Error in orbital overlap calculation: {e}")
            overlap_analysis['error'] = str(e)
        
        return overlap_analysis
    
    def _analyze_kernel_ci_coefficients(self, ci_coefficients) -> Dict[str, Any]:
        """
        Enhanced analysis of CI coefficients directly from kernel() return value.
        This provides more detailed analysis than the existing CI coefficient methods.
        """
        ci_analysis = {
            'source': 'kernel_return',
            'available': True
        }
        
        if ci_coefficients is None:
            return {'available': False, 'reason': 'No CI coefficients from kernel'}
        
        try:
            # Handle different CI coefficient formats
            if hasattr(ci_coefficients, 'flatten'):
                ci_flat = ci_coefficients.flatten()
            elif isinstance(ci_coefficients, (list, tuple)):
                ci_flat = np.array(ci_coefficients).flatten()
            else:
                ci_flat = np.array([ci_coefficients]).flatten()
            
            # Basic statistics
            ci_analysis['total_coefficients'] = len(ci_flat)
            ci_analysis['ci_vector_shape'] = str(getattr(ci_coefficients, 'shape', 'scalar'))
            
            # Compute squared coefficients for probability analysis
            ci_squared = np.abs(ci_flat)**2
            total_norm = np.sum(ci_squared)
            
            if total_norm > 0:
                # Normalize if not already normalized
                ci_squared_norm = ci_squared / total_norm
                
                # Find dominant configurations
                sorted_indices = np.argsort(ci_squared)[::-1]
                
                # Extract major configurations (>0.5% contribution)
                major_configs = []
                cumulative_weight = 0.0
                threshold = 0.005  # 0.5% threshold
                
                for i, idx in enumerate(sorted_indices[:50]):  # Limit to top 50
                    coeff = float(ci_flat[idx])
                    contribution = float(ci_squared_norm[idx])
                    cumulative_weight += contribution
                    
                    if contribution >= threshold:
                        major_configs.append({
                            'configuration_index': int(idx),
                            'coefficient': coeff,
                            'contribution_percent': contribution * 100.0,
                            'cumulative_percent': cumulative_weight * 100.0
                        })
                    
                    # Stop when we've captured 99% of the wavefunction
                    if cumulative_weight >= 0.99:
                        break
                
                ci_analysis.update({
                    'major_configurations': major_configs,
                    'leading_coefficient': float(ci_flat[sorted_indices[0]]),
                    'leading_contribution_percent': float(ci_squared_norm[sorted_indices[0]] * 100),
                    'multiconfigurational_character': 100.0 - float(ci_squared_norm[sorted_indices[0]] * 100),
                    'effective_configurations': len(major_configs),
                    'wavefunction_entropy': float(-np.sum(ci_squared_norm * np.log(ci_squared_norm + 1e-12))),
                    'normalization': float(total_norm)
                })
                
                # Configuration type analysis
                if len(major_configs) == 1:
                    dominant_character = "single_configuration"
                elif len(major_configs) <= 5:
                    dominant_character = "few_configuration"
                else:
                    dominant_character = "multiconfigurational"
                
                ci_analysis['wavefunction_character'] = dominant_character
                
                logger.info(f"Enhanced CI analysis: {len(major_configs)} major configurations, "
                           f"leading: {ci_analysis['leading_contribution_percent']:.1f}%, "
                           f"character: {dominant_character}")
            else:
                logger.warning("CI coefficients have zero norm")
                ci_analysis['error'] = 'zero_norm'
                
        except Exception as e:
            logger.error(f"Error in enhanced CI coefficient analysis: {e}")
            ci_analysis['error'] = str(e)
        
        return ci_analysis