"""CASCI calculator implementation using PySCF."""

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


class CASCICalculator(BaseCalculator):
    """CASCI calculator using PySCF for multi-configurational calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, keep_files: bool = False, molecule_name: Optional[str] = None):
        # Use file manager for better organization
        self.file_manager = CalculationFileManager()
        if working_dir is None:
            working_dir = self.file_manager.create_calculation_dir(molecule_name)
        super().__init__(working_dir)
        self.mol: Optional[gto.Mole] = None
        self.mf: Optional[scf.hf.SCF] = None  # Can be RHF or UHF
        self.mycas: Optional[mcscf.casci.CASCI] = None
        self.keep_files = keep_files
        self.molecule_name = molecule_name
        
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup CASCI calculation with molecular geometry and parameters."""
        try:
            # Extract basic calculation parameters
            basis = kwargs.get('basis', '6-31G(d)')
            charge = kwargs.get('charge', 0)
            spin = kwargs.get('spin', 0)
            max_cycle = kwargs.get('max_cycle', 150)
            solvent_method = kwargs.get('solvent_method', 'none')
            solvent = kwargs.get('solvent', '-')
            memory_mb = kwargs.get('memory_mb', 4000)  # Default 4GB
            
            # CASCI-specific parameters
            ncas = kwargs.get('ncas', 6)  # Number of active space orbitals
            nelecas = kwargs.get('nelecas', 8)  # Number of active space electrons
            analyze_nto = kwargs.get('analyze_nto', False)  # Natural transition orbitals
            natorb = kwargs.get('natorb', True)  # Transform to natural orbitals
            
            # Validate CASCI parameters
            if ncas <= 0:
                raise InputError(f"Number of active space orbitals (ncas={ncas}) must be positive")
            if nelecas <= 0:
                raise InputError(f"Number of active space electrons (nelecas={nelecas}) must be positive")
            if nelecas > 2 * ncas:
                raise InputError(f"Too many electrons ({nelecas}) for active space size ({ncas} orbitals)")
            
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
            
            # Apply memory settings
            if memory_mb and memory_mb > 0:
                self.mol.max_memory = memory_mb
            else:
                self.mol.max_memory = 4000  # CASCI needs substantial memory
            
            # Setup SCF calculation (prerequisite for CASCI)
            # For closed-shell systems (spin=0), use RHF
            # For open-shell systems (spin>0), use UHF
            if spin == 0:
                self.mf = scf.RHF(self.mol)
                logger.info("Using Restricted Hartree-Fock (RHF) for closed-shell CASCI reference")
            else:
                self.mf = scf.UHF(self.mol)
                logger.info("Using Unrestricted Hartree-Fock (UHF) for open-shell CASCI reference")
            
            self.mf = setup_solvent_effects(self.mf, solvent_method, solvent)
            
            self.mf.chkfile = self.get_checkpoint_path()
            self.mf.max_cycle = max_cycle
            
            # Store parameters for template method
            self.max_cycle = max_cycle
            self.solvent_method = solvent_method
            self.solvent = solvent
            self.ncas = ncas
            self.nelecas = nelecas
            self.analyze_nto = analyze_nto
            self.natorb = natorb
            
            # Store parameters in results for reference
            self.results.update({
                'basis': basis,
                'charge': charge,
                'spin': spin,
                'max_cycle': max_cycle,
                'solvent_method': solvent_method,
                'solvent': solvent,
                'atom_count': len(atoms),
                'ncas': ncas,
                'nelecas': nelecas,
                'analyze_nto': analyze_nto,
                'natorb': natorb,
                'method': 'UHF-CASCI' if spin > 0 else 'RHF-CASCI'
            })
            
            logger.info(f"CASCI setup completed: {ncas} orbitals, {nelecas} electrons in active space")
            
        except Exception as e:
            raise InputError(f"Failed to setup CASCI calculation: {str(e)}")
    
    # ===== Template Method Pattern Implementation =====
    
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform CASCI-specific calculation after SCF."""
        # Safely convert base SCF energy
        if base_energy is None:
            raise CalculationError("SCF energy is None - calculation may have failed")
        try:
            scf_energy_float = float(base_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert SCF energy to float: {e}")
        
        logger.info(f"Reference SCF completed. SCF energy: {scf_energy_float:.6f} hartree")
        
        # CASCI calculation
        logger.info("Starting CASCI calculation...")
        ncas = getattr(self, 'ncas', 6)
        nelecas = getattr(self, 'nelecas', 8)
        natorb = getattr(self, 'natorb', True)
        
        # Create CASCI object
        self.mycas = mcscf.CASCI(self.mf, ncas, nelecas)
        self.mycas.natorb = natorb
        
        # Validate active space parameters against molecular orbital structure
        if hasattr(self.mf, 'mo_energy') and self.mf.mo_energy is not None:
            n_orb = len(self.mf.mo_energy)
            n_occupied = int(np.sum(self.mf.mo_occ > 0)) if hasattr(self.mf, 'mo_occ') else 0
            
            if ncas > n_orb:
                raise InputError(f"Active space size ({ncas}) exceeds total orbitals ({n_orb})")
            if ncas > n_occupied + 10:  # Reasonable check for too many virtual orbitals
                logger.warning(f"Large active space ({ncas} orbitals) includes many virtual orbitals. "
                             f"System has {n_occupied} occupied orbitals out of {n_orb} total.")
        
        logger.info(f"CASCI setup: {ncas} active orbitals, {nelecas} active electrons")
        if natorb:
            logger.info("Natural orbital transformation enabled")
        
        # Run CASCI calculation
        logger.info("Starting CASCI kernel calculation...")
        try:
            casci_result = self.mycas.kernel()
        except Exception as e:
            error_msg = str(e).lower()
            if "singular" in error_msg or "convergence" in error_msg:
                raise ConvergenceError(f"CASCI calculation failed to converge: {str(e)}")
            elif "memory" in error_msg:
                raise CalculationError(f"CASCI calculation failed due to insufficient memory: {str(e)}")
            else:
                raise CalculationError(f"CASCI calculation failed: {str(e)}")
        
        # Validate results
        if casci_result is None:
            raise CalculationError("CASCI calculation failed: no result obtained")
        
        # Process kernel() return value comprehensively
        casci_energy = None
        ci_coefficients = None
        kernel_additional_info = {}
        
        if isinstance(casci_result, tuple):
            logger.info(f"CASCI kernel returned tuple with {len(casci_result)} elements")
            
            # Extract energy (always first element)
            if len(casci_result) >= 1:
                casci_energy = casci_result[0]
                logger.info("Extracted energy from kernel result tuple")
            
            # Extract CI coefficients (typically second element)
            if len(casci_result) >= 2:
                ci_coefficients = casci_result[1]
                logger.info("Extracted CI coefficients from kernel result")
                kernel_additional_info['ci_coefficients_shape'] = str(getattr(ci_coefficients, 'shape', 'no_shape'))
            
            # Store additional information for analysis
            kernel_additional_info['tuple_length'] = len(casci_result)
            for i, item in enumerate(casci_result[2:], start=2):
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
            casci_energy = casci_result
            logger.info("CASCI kernel returned single value (energy)")
        
        # Validate energy
        if casci_energy is None:
            raise CalculationError("CASCI calculation failed: no energy obtained")
        
        try:
            casci_energy_float = float(casci_energy)
        except (ValueError, TypeError) as e:
            raise CalculationError(f"Failed to convert CASCI energy to float: {e}")
        
        # Store CI coefficients for analysis
        if ci_coefficients is not None:
            self.ci_coefficients = ci_coefficients
        elif hasattr(self.mycas, 'ci') and self.mycas.ci is not None:
            self.ci_coefficients = self.mycas.ci
            logger.info("Using CI coefficients from CASCI object")
        else:
            self.ci_coefficients = None
        
        logger.info(f"CASCI kernel result processing complete. Energy: {casci_energy_float:.6f} hartree")
        
        if not self.mycas.converged:
            logger.warning("CASCI calculation did not converge to specified tolerance")
        
        # Analyze CASCI results (including CI coefficients)
        casci_results = self._analyze_casci_results()
        
        # Enhanced CI coefficients analysis if available
        enhanced_ci_analysis = {}
        if self.ci_coefficients is not None:
            enhanced_ci_analysis = self._analyze_kernel_ci_coefficients(self.ci_coefficients)
        
        # Combine results with comprehensive information
        results = {
            'scf_energy': scf_energy_float,
            'casci_energy': casci_energy_float,
            'correlation_energy': casci_energy_float - scf_energy_float,
            'kernel_return_info': kernel_additional_info,
            'ci_coefficients_available': self.ci_coefficients is not None,
            'enhanced_ci_analysis': enhanced_ci_analysis
        }
        results.update(casci_results)
        
        logger.info(f"CASCI calculation completed with enhanced analysis. Energy: {casci_energy_float:.6f} hartree")
        return results
    
    def _create_scf_method(self, mol):
        """Create SCF method object for CASCI reference (RHF/UHF)."""
        spin = self.results.get('spin', 0)
        
        if spin == 0:
            mf = scf.RHF(mol)
            logger.info("Using Restricted Hartree-Fock (RHF) for closed-shell CASCI reference")
        else:
            mf = scf.UHF(mol)
            logger.info("Using Unrestricted Hartree-Fock (UHF) for open-shell CASCI reference")
        
        return mf
    
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to SCF method."""
        return setup_solvent_effects(mf, self.solvent_method, self.solvent)
    
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        spin = self.results.get('spin', 0)
        return f"{'UHF' if spin > 0 else 'RHF'} (CASCI reference)"
    
    def _requires_geometry_optimization(self) -> bool:
        """CASCI does not require geometry optimization by default."""
        return False
    
    def _requires_frequency_analysis(self) -> bool:
        """CASCI does not require frequency analysis by default."""
        return False
    
    # ===== CASCI-specific Analysis Methods =====
    
    def _analyze_casci_results(self) -> Dict[str, Any]:
        """
        Perform comprehensive analysis of CASCI results.
        """
        if self.mycas is None:
            raise CalculationError("CASCI calculation not available for analysis")
        
        logger.info("Starting CASCI results analysis...")
        
        # Set verbose level for detailed analysis
        original_verbose = self.mycas.verbose
        self.mycas.verbose = 4
        
        # Perform PySCF's built-in analysis
        logger.info("Running PySCF CASCI analysis...")
        try:
            self.mycas.analyze()
        except Exception as e:
            logger.warning(f"PySCF CASCI analysis encountered issues: {e}")
        
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
        
        logger.info("CASCI results analysis completed")
        return analysis_results
    
    def _analyze_natural_orbitals(self) -> Dict[str, Any]:
        """
        Analyze natural orbitals and their occupation numbers.
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
        """
        ci_analysis = {}
        
        if not hasattr(self.mycas, 'ci') or self.mycas.ci is None:
            logger.warning("No CI coefficients available for analysis")
            return {'available': False, 'reason': 'No CI coefficients found'}
        
        ci_analysis['available'] = True
        
        try:
            # For CASCI, ci is typically a numpy array
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
        """
        if not hasattr(self.mycas, 'make_rdm1s') or not callable(self.mycas.make_rdm1s):
            logger.warning("CASCI object does not support spin density calculation")
            return {'available': False, 'reason': 'Spin density not supported for this CASCI method'}
        
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
        Calculate overlap between CASCI orbitals and reference SCF orbitals.
        """
        overlap_analysis = {}
        
        if not hasattr(self.mycas, 'mo_coeff') or self.mycas.mo_coeff is None:
            logger.warning("CASCI orbitals not available for overlap analysis")
            return {'available': False, 'reason': 'CASCI orbitals not found'}
        
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
            
            casci_orbs = self.mycas.mo_coeff
            scf_orbs = self.mf.mo_coeff
            
            # Calculate orbital overlap matrix: <CASCI_i|SCF_j>
            # Overlap = C_CASCI^T * S * C_SCF
            overlap_matrix = np.dot(casci_orbs.T, np.dot(S, scf_orbs))
            
            # Analyze overlap for active space orbitals
            ncas = self.results.get('ncas', 6)
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