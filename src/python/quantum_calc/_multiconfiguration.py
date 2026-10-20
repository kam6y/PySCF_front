from __future__ import annotations

import logging
from typing import Any, Dict, Tuple

import numpy as np

from .exceptions import CalculationError

logger = logging.getLogger(__name__)

class MulticonfigurationAnalysisMixin:
    """Mixin for shared CASCI/CASSCF multiconfiguration analyses."""
    def _validate_active_space_parameters(self, ncas: int, nelecas: int) -> Tuple[int, int]:
        """
        Validate and adjust CASCI/CASSCF active space parameters.
        
        Args:
            ncas: Number of active space orbitals
            nelecas: Number of active space electrons
            
        Returns:
            Tuple of (adjusted_ncas, adjusted_nelecas)
        """
        # Validate and adjust ncas
        if ncas <= 0:
            logger.warning(f"Invalid ncas={ncas}, adjusting to default ncas=4")
            ncas = 4
            
        # Validate and adjust nelecas
        if nelecas <= 0:
            logger.warning(f"Invalid nelecas={nelecas}, adjusting to default nelecas=4")
            nelecas = 4
            
        # Check electron count vs orbital count
        if nelecas > 2 * ncas:
            logger.warning(f"Too many electrons ({nelecas}) for active space size ({ncas} orbitals)")
            # Adjust nelecas to maximum possible for the given ncas
            nelecas = 2 * ncas
            logger.info(f"Adjusted nelecas to maximum possible: {nelecas}")
            
        # Additional sanity checks with warnings
        if ncas > 20:
            logger.warning(f"Large active space (ncas={ncas}) may require substantial computational resources")
        if nelecas > 20:
            logger.warning(f"Many active electrons (nelecas={nelecas}) may require substantial computational resources")
            
        return ncas, nelecas

    def _handle_calculation_error(self, e: Exception, calculation_type: str) -> None:
        """
        Common error handling for calculation failures.
        
        Args:
            e: The exception that occurred
            calculation_type: Type of calculation (e.g., 'CASCI', 'CASSCF', 'DFT')
        """
        import traceback
        from .exceptions import CalculationError, ConvergenceError
        
        # Enhanced error logging
        logger.error(f"Calculation failed in {calculation_type} run_calculation:")
        logger.error(f"Exception type: {type(e).__name__}")
        logger.error(f"Exception message: '{str(e)}'")
        logger.error(f"Full traceback:\n{traceback.format_exc()}")
        
        # Handle empty error messages
        error_message = str(e)
        if not error_message.strip():
            error_message = f"Unknown error in {calculation_type} calculation (empty error message)"
            logger.error(f"Empty error message detected, using: {error_message}")
        
        error_msg_lower = error_message.lower()
        
        # Classify error types
        if "singular" in error_msg_lower or "convergence" in error_msg_lower:
            raise ConvergenceError(f"{calculation_type} calculation failed to converge: {error_message}")
        elif "memory" in error_msg_lower:
            raise CalculationError(f"{calculation_type} calculation failed due to insufficient memory: {error_message}")
        elif "maximum" in error_msg_lower and "cycle" in error_msg_lower:
            raise ConvergenceError(f"{calculation_type} reached maximum cycles without convergence: {error_message}")
        elif isinstance(e, (CalculationError, ConvergenceError)):
            raise
        else:
            raise CalculationError(f"{calculation_type} calculation failed: {error_message}")

    def _analyze_kernel_return_value(self, kernel_result, calculation_type: str) -> Tuple[float, Any, Dict[str, Any]]:
        """
        Common analysis of PySCF kernel() return values.
        
        Args:
            kernel_result: The result returned by PySCF kernel() method
            calculation_type: Type of calculation for error messages
            
        Returns:
            Tuple of (energy, ci_coefficients, additional_info)
        """
        energy = None
        ci_coefficients = None
        additional_info = {}
        
        if kernel_result is None:
            raise CalculationError(f"{calculation_type} calculation failed: no result obtained")
        
        if isinstance(kernel_result, tuple):
            logger.info(f"{calculation_type} kernel returned tuple with {len(kernel_result)} elements")
            
            # Extract energy (always first element)
            if len(kernel_result) >= 1:
                energy = kernel_result[0]
                logger.info(f"Extracted energy from {calculation_type} kernel result tuple")
            
            # Extract CI coefficients (typically second element)
            if len(kernel_result) >= 2:
                ci_coefficients = kernel_result[1]
                logger.info(f"Extracted CI coefficients from {calculation_type} kernel result")
                additional_info['ci_coefficients_shape'] = str(getattr(ci_coefficients, 'shape', 'no_shape'))
            
            # Store additional information for analysis
            additional_info['tuple_length'] = len(kernel_result)
            for i, item in enumerate(kernel_result[2:], start=2):
                try:
                    if hasattr(item, 'shape'):
                        additional_info[f'element_{i}_shape'] = str(item.shape)
                    elif isinstance(item, (int, float, complex)):
                        additional_info[f'element_{i}_value'] = float(item)
                    else:
                        additional_info[f'element_{i}_type'] = type(item).__name__
                except Exception:
                    additional_info[f'element_{i}_type'] = 'unknown'
        else:
            # Single value returned
            energy = kernel_result
            logger.info(f"{calculation_type} kernel returned single value (energy)")
        
        # Validate energy
        if energy is None:
            raise CalculationError(f"{calculation_type} calculation failed: no energy obtained")
        
        try:
            energy_float = float(energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert {calculation_type} energy to float: {e}")
        
        logger.info(f"{calculation_type} kernel result processing complete. Energy: {energy_float:.6f} hartree")
        
        return energy_float, ci_coefficients, additional_info

    def _analyze_natural_orbitals(self) -> Dict[str, Any]:
        """
        Analyze natural orbitals and their occupation numbers.
        Common method for CASCI and CASSCF calculations.
        """
        if not hasattr(self, 'mycas') or self.mycas is None:
            logger.warning("CASCI/CASSCF object not available for natural orbital analysis")
            return {'enabled': False, 'reason': 'CASCI/CASSCF object not available'}
        
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
        Common method for CASCI and CASSCF calculations.
        """
        ci_analysis = {}
        
        if not hasattr(self, 'mycas') or self.mycas is None:
            logger.warning("CASCI/CASSCF object not available for CI coefficient analysis")
            return {'available': False, 'reason': 'CASCI/CASSCF object not available'}
        
        if not hasattr(self.mycas, 'ci') or self.mycas.ci is None:
            logger.warning("No CI coefficients available for analysis")
            return {'available': False, 'reason': 'No CI coefficients found'}
        
        ci_analysis['available'] = True
        
        try:
            # For CASCI/CASSCF, ci is typically a numpy array
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

    def _calculate_mulliken_spin_density(self) -> Dict[str, Any]:
        """
        Calculate Mulliken atomic spin densities for open-shell systems.
        Common method for CASCI and CASSCF calculations.
        """
        if not hasattr(self, 'mycas') or self.mycas is None:
            logger.warning("CASCI/CASSCF object not available for spin density calculation")
            return {'available': False, 'reason': 'CASCI/CASSCF object not available'}
        
        if not hasattr(self.mycas, 'make_rdm1s') or not callable(self.mycas.make_rdm1s):
            logger.warning("CASCI/CASSCF object does not support spin density calculation")
            return {'available': False, 'reason': 'Spin density not supported for this method'}
        
        spin_analysis = {'available': True}
        
        try:
            # Get spin density matrices
            dm1a, dm1b = self.mycas.make_rdm1s()  # Alpha and beta density matrices
            
            if dm1a is None or dm1b is None:
                logger.warning("Could not obtain alpha/beta density matrices")
                return {'available': False, 'reason': 'Density matrices not available'}
            
            # Calculate spin density matrix (alpha - beta)
            spin_dm = dm1a - dm1b

            # Use self.mf.mol to ensure consistency with the converged calculation
            mol = self.mf.mol if hasattr(self, 'mf') and self.mf is not None else self.mol

            # Perform Mulliken population analysis on spin density
            if hasattr(mol, 'get_ovlp'):
                ovlp = mol.get_ovlp()
            else:
                ovlp = self.mycas._scf.get_ovlp()

            # Mulliken spin populations
            spin_pop = np.einsum('ij,ji->i', spin_dm, ovlp)

            # Group by atoms
            atomic_spin_densities = []
            ao_offset = 0

            for iatom in range(mol.natm):
                atom_symbol = mol.atom_symbol(iatom)
                nao = mol.aoslice_by_atom()[iatom][3] - mol.aoslice_by_atom()[iatom][2]

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
                'expected_spin': float(mol.spin)
            })

            logger.info(f"Spin density analysis: total = {total_spin_density:.3f}, "
                       f"expected = {mol.spin}, atoms analyzed = {len(atomic_spin_densities)}")
            
        except Exception as e:
            logger.error(f"Error in Mulliken spin density calculation: {e}")
            spin_analysis['error'] = str(e)
        
        return spin_analysis

    def _calculate_orbital_overlap(self) -> Dict[str, Any]:
        """
        Calculate overlap between CASCI/CASSCF orbitals and reference SCF orbitals.
        Common method for CASCI and CASSCF calculations.
        """
        overlap_analysis = {}
        
        if not hasattr(self, 'mycas') or self.mycas is None:
            logger.warning("CASCI/CASSCF object not available for overlap analysis")
            return {'available': False, 'reason': 'CASCI/CASSCF object not available'}
        
        if not hasattr(self.mycas, 'mo_coeff') or self.mycas.mo_coeff is None:
            logger.warning("CASCI/CASSCF orbitals not available for overlap analysis")
            return {'available': False, 'reason': 'CASCI/CASSCF orbitals not found'}
        
        if not hasattr(self.mf, 'mo_coeff') or self.mf.mo_coeff is None:
            logger.warning("SCF orbitals not available for overlap analysis")
            return {'available': False, 'reason': 'SCF reference orbitals not found'}
        
        overlap_analysis['available'] = True

        # Use self.mf.mol to ensure consistency with the converged calculation
        mol = self.mf.mol if hasattr(self, 'mf') and self.mf is not None else self.mol

        try:
            # Get overlap matrix
            S = mol.get_ovlp() if hasattr(mol, 'get_ovlp') else self.mf.get_ovlp()
            
            cas_orbs = self.mycas.mo_coeff
            scf_orbs = self.mf.mo_coeff
            
            # Calculate orbital overlap matrix: <CAS_i|SCF_j>
            # Overlap = C_CAS^T * S * C_SCF
            overlap_matrix = np.dot(cas_orbs.T, np.dot(S, scf_orbs))
            
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

    def _analyze_enhanced_ci_coefficients(self, ci_coefficients) -> Dict[str, Any]:
        """
        Enhanced analysis of CI coefficients directly from kernel() return value.
        This provides more detailed analysis than the existing CI coefficient methods.
        Common method for CASCI and CASSCF calculations.
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
