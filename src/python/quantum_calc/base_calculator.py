"""Base calculator class for quantum chemistry calculations."""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Tuple
import os
import tempfile
import logging
from contextlib import contextmanager
import numpy as np

logger = logging.getLogger(__name__)


class BaseCalculator(ABC):
    """Abstract base class for quantum chemistry calculations."""
    
    def __init__(self, working_dir: Optional[str] = None):
        """Initialize calculator with optional working directory."""
        self.working_dir = working_dir or tempfile.mkdtemp(prefix="pyscf_calc_")
        self.results: Dict[str, Any] = {}
        
    def parse_xyz(self, xyz_string: str) -> List[List]:
        """Parse XYZ format string into atom list."""
        lines = xyz_string.strip().split('\n')
        if len(lines) < 3:
            raise ValueError("Invalid XYZ format: insufficient lines")
        
        try:
            atom_count = int(lines[0])
        except ValueError:
            raise ValueError("Invalid XYZ format: first line must be atom count")
        
        if len(lines) < atom_count + 2:
            raise ValueError(f"Invalid XYZ format: expected {atom_count + 2} lines, got {len(lines)}")
        
        atoms = []
        for i in range(2, atom_count + 2):
            parts = lines[i].split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ format at line {i + 1}: insufficient columns")
            
            symbol = parts[0]
            try:
                coords = [float(parts[j]) for j in range(1, 4)]
            except ValueError:
                raise ValueError(f"Invalid XYZ format at line {i + 1}: invalid coordinates")
            
            atoms.append([symbol, coords])
        
        return atoms
    
    def get_checkpoint_path(self) -> str:
        """Get the path to the checkpoint file."""
        return os.path.join(self.working_dir, "calculation.chk")
    
    def apply_resource_settings(self, mol, memory_mb: Optional[int] = None, cpu_cores: Optional[int] = None) -> None:
        """Apply resource settings to PySCF molecule object."""
        if memory_mb is not None and memory_mb > 0:
            # PySCF expects memory in MB
            mol.max_memory = int(memory_mb)
            print(f"Set PySCF max_memory to {memory_mb} MB")
        else:
            # Use calculation-specific default memory settings
            calculation_method = getattr(self, 'calculation_method', 'DFT')
            if calculation_method in ['CASCI', 'CASSCF']:
                mol.max_memory = 6000  # 6GB for CASCI/CASSCF
                print("Using default PySCF max_memory: 6000 MB (CASCI/CASSCF)")
            elif calculation_method in ['CCSD', 'CCSD_T']:
                mol.max_memory = 4000  # 4GB for coupled cluster
                print("Using default PySCF max_memory: 4000 MB (CCSD)")
            else:
                mol.max_memory = 2000  # 2GB for other methods
                print("Using default PySCF max_memory: 2000 MB")
        
        # CPU cores are now configured at the process level in process_manager.py
        # This avoids conflicts and ensures proper timing of environment variable setup
    
    @contextmanager
    def controlled_threading(self, cpu_cores: Optional[int] = None):
        """
        Context manager to control BLAS/LAPACK threading for specific operations.
        
        Args:
            cpu_cores: Number of threads to use for BLAS/LAPACK operations
        """
        try:
            from threadpoolctl import threadpool_limits
            
            if cpu_cores is not None and cpu_cores > 0:
                print(f"Applying threadpool_limits(limits={cpu_cores}, user_api='blas')")
                with threadpool_limits(limits=int(cpu_cores), user_api='blas'):
                    yield
            else:
                # No thread control - proceed normally
                yield
        except ImportError:
            print("threadpoolctl not available, proceeding without thread control")
            yield
        except Exception as e:
            print(f"Error in thread control: {e}, proceeding without thread control")
            yield
    
    @abstractmethod
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup the quantum chemistry calculation."""
        pass
    
    def run_calculation(self) -> Dict[str, Any]:
        """Template method for running quantum chemistry calculations."""
        try:
            self._pre_calculation_check()
            
            if self._requires_geometry_optimization():
                self._perform_geometry_optimization()
            
            self._setup_final_calculation()
            base_energy = self._run_base_scf_calculation()
            self._verify_scf_convergence()
            
            specific_results = self._perform_specific_calculation(base_energy)
            
            if self._requires_orbital_analysis():
                self._perform_orbital_analysis()
            
            if self._requires_frequency_analysis():
                self._perform_frequency_analysis()
            
            return self._prepare_final_results(specific_results)
            
        except Exception as e:
            import traceback
            from .exceptions import CalculationError
            
            # Enhanced error logging for CASCI/CASSCF
            calculation_method = getattr(self, 'calculation_method', 'Unknown')
            logger.error(f"Calculation failed in {calculation_method} run_calculation:")
            logger.error(f"Exception type: {type(e).__name__}")
            logger.error(f"Exception message: '{str(e)}'")
            logger.error(f"Full traceback:\n{traceback.format_exc()}")
            
            # Handle empty error messages specifically for CASCI/CASSCF
            error_message = str(e)
            if not error_message.strip():
                error_message = f"Unknown error in {calculation_method} calculation (empty error message)"
                logger.error(f"Empty error message detected, using: {error_message}")
            
            if isinstance(e, (CalculationError,)):
                raise
            raise CalculationError(f"{calculation_method} calculation failed: {error_message}")
    
    def cleanup(self, keep_files: bool = False) -> None:
        """Clean up temporary files."""
        import shutil
        if not keep_files and os.path.exists(self.working_dir) and "pyscf_calc_" in self.working_dir:
            shutil.rmtree(self.working_dir)
        elif keep_files:
            print(f"Calculation files preserved in: {self.working_dir}")
    
    def _atoms_to_string(self, atoms: List[List]) -> str:
        """Convert atoms list to PySCF atom string format."""
        from .exceptions import GeometryError
        
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
        from .exceptions import CalculationError
        
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mo_occ is None:
            raise CalculationError("Orbital occupations not available")
        
        # Handle both RKS/RHF (1D array) and UKS/UHF (2D array) cases
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS/UHF case: use alpha orbitals
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
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mo_occ is None:
            return 0
        
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS/UHF case: count both alpha and beta orbitals
            return int(np.sum(mo_occ > 0))
        else:
            # RKS/RHF case: simple sum
            return int(np.sum(mo_occ > 0))
    
    def _count_virtual_orbitals(self) -> int:
        """Count the number of virtual orbitals."""
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mo_occ is None:
            return 0
        
        mo_occ = self.mf.mo_occ
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS/UHF case: count both alpha and beta orbitals
            return int(np.sum(mo_occ == 0))
        else:
            # RKS/RHF case: simple sum
            return int(np.sum(mo_occ == 0))
    
    def _geometry_to_xyz_string(self) -> str:
        """Convert optimized geometry to XYZ format string."""
        if not hasattr(self, 'optimized_geometry') or self.optimized_geometry is None:
            return ""
        if not hasattr(self, 'mol') or self.mol is None:
            return ""
        
        atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
        lines = [str(self.mol.natm)]
        lines.append("Optimized geometry from PySCF calculation")
        
        for i, (symbol, coords) in enumerate(zip(atom_symbols, self.optimized_geometry)):
            lines.append(f"{symbol:2s} {coords[0]:12.6f} {coords[1]:12.6f} {coords[2]:12.6f}")
        
        return "\n".join(lines)
    
    def _calculate_mulliken_charges(self) -> Optional[List[Dict[str, Any]]]:
        """Calculate Mulliken population analysis charges for each atom."""
        if not hasattr(self, 'mf') or self.mf is None:
            return None
        if not hasattr(self, 'mol') or self.mol is None:
            return None
        
        try:
            # Perform Mulliken population analysis
            # This returns (pop, charges) where pop are populations and charges are atomic charges
            pop, charges = self.mf.mulliken_pop()
            
            # Extract charges for each atom
            mulliken_charges = []
            for i in range(self.mol.natm):
                atom_symbol = self.mol.atom_symbol(i)
                # Convert numpy float to Python float for JSON serialization
                charge = float(charges[i])
                
                mulliken_charges.append({
                    'atom_index': i,
                    'element': atom_symbol,
                    'charge': charge
                })
            
            # Verify total charge conservation (should equal molecular charge)
            total_charge = sum(item['charge'] for item in mulliken_charges)
            expected_charge = float(self.mol.charge)
            logger.info(f"Mulliken analysis: calculated total charge = {total_charge:.6f}, expected = {expected_charge:.6f}")
            
            if abs(total_charge - expected_charge) > 0.01:
                logger.warning(f"Mulliken total charge ({total_charge:.6f}) differs from molecular charge ({expected_charge:.6f}) by more than 0.01")
            
            return mulliken_charges
            
        except Exception as e:
            logger.warning(f"Mulliken population analysis failed: {str(e)}")
            return None
    
    def _safe_float_conversion(self, value: Any, param_name: str) -> Optional[float]:
        """
        Safely convert PySCF thermochemistry values to float.
        
        Args:
            value: Value from PySCF thermo function (could be scalar, tuple, array)
            param_name: Parameter name for logging purposes
            
        Returns:
            Float value or None if conversion fails
        """
        try:
            # Handle different data types returned by PySCF
            if isinstance(value, (list, tuple)):
                # Take the first element if it's a sequence
                if len(value) > 0:
                    return float(value[0])
                else:
                    logger.warning(f"Empty sequence for {param_name}")
                    return None
            elif hasattr(value, '__iter__') and not isinstance(value, str):
                # Handle numpy arrays and other iterables
                try:
                    import numpy as np
                    if isinstance(value, np.ndarray):
                        return float(value.flat[0])  # Get first element safely
                    else:
                        # Generic iterable
                        first_val = next(iter(value))
                        return float(first_val)
                except (StopIteration, IndexError):
                    logger.warning(f"Empty iterable for {param_name}")
                    return None
            else:
                # Direct scalar conversion
                return float(value)
        except (ValueError, TypeError) as e:
            logger.warning(f"Failed to convert {param_name} (type: {type(value)}, value: {value}): {e}")
            return None
    
    def _calculate_frequencies(self) -> Optional[Dict[str, Any]]:
        """Calculate vibrational frequencies and thermochemical properties."""
        if not hasattr(self, 'mf') or self.mf is None:
            return None
        if not hasattr(self, 'mol') or self.mol is None:
            return None
        
        # Initialize partial result structure
        frequency_results = {
            'frequency_analysis_performed': False,
            'vibrational_frequencies': None,
            'imaginary_frequencies_count': None,
            'zero_point_energy': None,
            'thermal_energy_298K': None,
            'entropy_298K': None,
            'gibbs_free_energy_298K': None,
            'heat_capacity_298K': None
        }
        
        try:
            logger.info("Starting vibrational frequency analysis...")
            
            # Import thermochemistry module early to catch import errors
            from pyscf.hessian import thermo
            
            # Step 1: Calculate Hessian matrix
            try:
                logger.info("Computing Hessian matrix...")
                hessian_calc = self.mf.Hessian()
                hessian = hessian_calc.kernel()
                logger.info("Hessian calculation completed successfully")
            except Exception as e:
                logger.error(f"Hessian calculation failed: {str(e)}")
                return frequency_results
            
            # Step 2: Perform harmonic analysis
            try:
                logger.info("Performing harmonic analysis...")
                freq_info = thermo.harmonic_analysis(self.mol, hessian)
                
                # Extract frequencies in cm^-1
                frequencies_au = freq_info['freq_au']
                frequencies_cm = freq_info['freq_wavenumber']
                
                # Filter out low frequencies (below 80 cm^-1 threshold)
                # and count imaginary frequencies
                imaginary_count = 0
                positive_frequencies = []
                
                for freq_cm in frequencies_cm:
                    if freq_cm < 0:
                        # Imaginary frequency (negative eigenvalue)
                        imaginary_count += 1
                    elif freq_cm >= 80.0:  # PySCF default threshold
                        # Real, significant frequency
                        positive_frequencies.append(float(freq_cm))
                
                logger.info(f"Found {len(positive_frequencies)} positive frequencies (≥80 cm⁻¹)")
                logger.info(f"Found {imaginary_count} imaginary frequencies")
                
                # Update results with frequency information
                frequency_results.update({
                    'frequency_analysis_performed': True,
                    'vibrational_frequencies': positive_frequencies,
                    'imaginary_frequencies_count': imaginary_count
                })
                
                # Log optimization quality assessment
                if imaginary_count == 0:
                    logger.info("Geometry optimization successful: no imaginary frequencies detected")
                elif imaginary_count == 1:
                    logger.warning("One imaginary frequency detected: possible transition state")
                else:
                    logger.warning(f"Multiple imaginary frequencies ({imaginary_count}) detected: poor optimization")
                    
            except Exception as e:
                logger.error(f"Harmonic analysis failed: {str(e)}")
                return frequency_results
            
            # Step 3: Calculate thermochemical properties (optional, non-critical)
            try:
                temperature = 298.15  # K
                pressure = 101325     # Pa
                
                logger.info("Computing thermochemical properties at 298.15 K...")
                thermo_info = thermo.thermo(self.mf, frequencies_au, temperature, pressure)
                
                # Safely extract thermochemical properties
                logger.debug(f"Thermo info keys: {list(thermo_info.keys())}")
                
                zero_point_energy = self._safe_float_conversion(thermo_info.get('ZPE'), 'ZPE')
                thermal_energy = self._safe_float_conversion(thermo_info.get('E_tot'), 'E_tot')
                entropy = self._safe_float_conversion(thermo_info.get('S_tot'), 'S_tot')
                gibbs_free_energy = self._safe_float_conversion(thermo_info.get('G_tot'), 'G_tot')
                heat_capacity = self._safe_float_conversion(thermo_info.get('Cv_tot'), 'Cv_tot')
                
                # Log successful conversions
                successful_conversions = []
                if zero_point_energy is not None:
                    successful_conversions.append('ZPE')
                if thermal_energy is not None:
                    successful_conversions.append('E_tot')
                if entropy is not None:
                    successful_conversions.append('S_tot')
                if gibbs_free_energy is not None:
                    successful_conversions.append('G_tot')
                if heat_capacity is not None:
                    successful_conversions.append('Cv_tot')
                
                logger.info(f"Successfully converted thermochemical properties: {successful_conversions}")
                
                # Update results with thermochemical data
                frequency_results.update({
                    'zero_point_energy': zero_point_energy,
                    'thermal_energy_298K': thermal_energy,
                    'entropy_298K': entropy,
                    'gibbs_free_energy_298K': gibbs_free_energy,
                    'heat_capacity_298K': heat_capacity
                })
                
                logger.info("Thermochemical analysis completed successfully")
                
            except Exception as e:
                logger.warning(f"Thermochemical analysis failed (frequency data preserved): {str(e)}")
                # Continue without thermochemical data - frequency data is still valuable
            
            logger.info("Vibrational frequency analysis completed successfully")
            return frequency_results
            
        except ImportError as e:
            logger.warning(f"Frequency analysis skipped: missing PySCF hessian module ({str(e)})")
            return None
        except Exception as e:
            logger.error(f"Vibrational frequency analysis failed completely: {str(e)}")
            # Return partial results if any were obtained
            if frequency_results['frequency_analysis_performed']:
                logger.info("Returning partial frequency analysis results")
                return frequency_results
            return None
    
    # ===== Template Method Pattern Implementation =====
    
    def _pre_calculation_check(self) -> None:
        """Check prerequisites for calculation."""
        if not hasattr(self, 'mol') or self.mol is None:
            raise CalculationError("Molecular object not setup. Call setup_calculation first.")
        if not hasattr(self, 'mf') or self.mf is None:
            raise CalculationError("Mean field object not setup. Call setup_calculation first.")
    
    def _perform_geometry_optimization(self) -> None:
        """Perform geometry optimization using geometric_solver."""
        from pyscf.geomopt import geometric_solver
        
        logger.info("Starting geometry optimization...")
        optimized_mol = geometric_solver.optimize(self.mf)
        self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
        logger.info("Geometry optimization completed")
        
        # Apply coordinate alignment after optimization
        self._align_optimized_geometry()
    
    def _align_optimized_geometry(self) -> None:
        """
        Align optimized geometry using ASE with correct rotation matrix approach:
        1. Move center of mass to origin
        2. Calculate principal axes using inertia tensor
        3. Create rotation matrix to align principal axes with coordinate axes
        4. Apply rotation matrix to all atomic coordinates
        """
        if not hasattr(self, 'optimized_geometry') or self.optimized_geometry is None:
            logger.warning("No optimized geometry available for alignment")
            return
        
        if not hasattr(self, 'mol') or self.mol is None:
            logger.warning("Molecular object not available for alignment")
            return
        
        try:
            # Import required libraries
            from ase import Atoms
            
            # Get atom symbols from the molecule
            atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
            positions = self.optimized_geometry.copy()
            
            logger.info("Starting molecular alignment with ASE...")
            logger.info(f"Original coordinates (first few atoms): {positions[:min(3, len(positions))]}")
            
            # Create ASE Atoms object
            atoms = Atoms(symbols=atom_symbols, positions=positions)
            
            # Move center of mass to origin
            atoms.center()
            centered_positions = atoms.get_positions()
            logger.info("Moved center of mass to origin")
            
            # Calculate moments of inertia and get principal axes
            inertia_moments, principal_axes = atoms.get_moments_of_inertia(vectors=True)
            logger.info(f"Original inertia moments: {inertia_moments}")
            logger.info(f"Principal axes shape: {principal_axes.shape}")
            
            # Check for special cases
            if len(atoms) == 1:
                logger.info("Single atom molecule - only centering applied")
                aligned_positions = centered_positions
                rotation_applied = False
            elif np.any(inertia_moments < 1e-10):
                logger.info("Linear molecule detected - limited alignment applied")
                # For linear molecules, align the molecular axis with Z-axis
                aligned_positions, rotation_applied = self._align_linear_molecule(centered_positions, principal_axes)
            else:
                # Apply rotation to align principal axes with coordinate axes
                logger.info("Applying rotation to align principal axes with coordinate axes")
                aligned_positions, rotation_applied = self._apply_rotation_matrix(centered_positions, principal_axes)
            
            # Update the optimized geometry with aligned coordinates
            self.optimized_geometry = aligned_positions
            
            # Verify alignment by checking the new inertia tensor
            atoms_aligned = Atoms(symbols=atom_symbols, positions=aligned_positions)
            new_inertia_moments = atoms_aligned.get_moments_of_inertia()
            
            if rotation_applied:
                # Calculate rotation angle for logging
                rotation_angle = self._calculate_rotation_angle(inertia_moments, new_inertia_moments)
                logger.info(f"Molecular alignment completed successfully (rotation angle: {rotation_angle:.2f}°)")
                logger.info(f"Aligned coordinates (first few atoms): {aligned_positions[:min(3, len(aligned_positions))]}")
                logger.info(f"Original inertia moments: {inertia_moments}")
                logger.info(f"Aligned inertia moments: {new_inertia_moments}")
                
                # Verify that rotation was meaningful
                if not self._validate_rotation(inertia_moments, new_inertia_moments):
                    logger.warning("Rotation validation failed - inertia moments unchanged")
            else:
                logger.info("No rotation applied - molecule already aligned or special case")
                logger.info(f"Final coordinates (first few atoms): {aligned_positions[:min(3, len(aligned_positions))]}")
            
        except ImportError as e:
            logger.error(f"ASE library not available for molecular alignment: {str(e)}")
            logger.info("Proceeding without molecular alignment")
        except Exception as e:
            logger.error(f"Molecular alignment failed: {str(e)}")
            logger.info("Proceeding with unaligned geometry")
    
    def _apply_rotation_matrix(self, positions: np.ndarray, principal_axes: np.ndarray) -> tuple[np.ndarray, bool]:
        """
        Apply rotation matrix to align principal axes with coordinate axes.
        
        Args:
            positions: Atomic positions (N x 3 array)
            principal_axes: Principal axes from inertia tensor (3 x 3 matrix)
            
        Returns:
            Tuple of (rotated_positions, rotation_applied_flag)
        """
        try:
            # The principal_axes matrix contains eigenvectors as columns
            # We want to rotate so that these axes align with [1,0,0], [0,1,0], [0,0,1]
            # This means we want: principal_axes @ rotation_matrix = identity
            # So: rotation_matrix = principal_axes^T
            rotation_matrix = principal_axes.T
            
            # Apply rotation to all atomic positions
            rotated_positions = positions @ rotation_matrix
            
            logger.info(f"Applied rotation matrix with determinant: {np.linalg.det(rotation_matrix):.6f}")
            
            return rotated_positions, True
            
        except Exception as e:
            logger.error(f"Failed to apply rotation matrix: {str(e)}")
            return positions, False
    
    def _align_linear_molecule(self, positions: np.ndarray, principal_axes: np.ndarray) -> tuple[np.ndarray, bool]:
        """
        Align linear molecule with Z-axis.
        
        Args:
            positions: Atomic positions (N x 3 array)  
            principal_axes: Principal axes from inertia tensor (3 x 3 matrix)
            
        Returns:
            Tuple of (aligned_positions, rotation_applied_flag)
        """
        try:
            # For linear molecules, find the axis with the largest moment of inertia
            # and align it with the Z-axis
            molecular_axis = principal_axes[:, -1]  # Last eigenvector (largest moment)
            z_axis = np.array([0, 0, 1])
            
            # Calculate rotation axis and angle
            rotation_axis = np.cross(molecular_axis, z_axis)
            rotation_angle = np.arccos(np.clip(np.dot(molecular_axis, z_axis), -1, 1))
            
            if rotation_angle < 1e-6:  # Already aligned
                logger.info("Linear molecule already aligned with Z-axis")
                return positions, False
            
            # Create rotation matrix using Rodrigues' formula
            rotation_axis_normalized = rotation_axis / np.linalg.norm(rotation_axis)
            cos_angle = np.cos(rotation_angle)
            sin_angle = np.sin(rotation_angle)
            
            # Rodrigues' rotation matrix formula
            K = np.array([[0, -rotation_axis_normalized[2], rotation_axis_normalized[1]],
                         [rotation_axis_normalized[2], 0, -rotation_axis_normalized[0]],
                         [-rotation_axis_normalized[1], rotation_axis_normalized[0], 0]])
            
            rotation_matrix = np.eye(3) + sin_angle * K + (1 - cos_angle) * K @ K
            
            # Apply rotation
            aligned_positions = positions @ rotation_matrix.T
            
            logger.info(f"Linear molecule aligned with Z-axis (rotation angle: {np.degrees(rotation_angle):.2f}°)")
            
            return aligned_positions, True
            
        except Exception as e:
            logger.error(f"Failed to align linear molecule: {str(e)}")
            return positions, False
    
    def _calculate_rotation_angle(self, original_moments: np.ndarray, aligned_moments: np.ndarray) -> float:
        """Calculate approximate rotation angle from inertia moment changes."""
        try:
            # Use the change in the smallest moment as an indicator
            original_sorted = np.sort(original_moments)
            aligned_sorted = np.sort(aligned_moments)
            
            # Calculate relative change
            relative_change = np.abs(original_sorted - aligned_sorted) / (original_sorted + 1e-10)
            max_change = np.max(relative_change)
            
            # Rough estimate: larger changes suggest larger rotations
            estimated_angle = np.degrees(np.arctan(max_change * 10))  # Heuristic scaling
            
            return min(estimated_angle, 180.0)  # Cap at 180 degrees
            
        except Exception:
            return 0.0
    
    def _validate_rotation(self, original_moments: np.ndarray, aligned_moments: np.ndarray) -> bool:
        """Validate that rotation was meaningful by checking inertia moment changes."""
        try:
            # Check if any moment changed significantly (more than 0.1%)
            relative_changes = np.abs(original_moments - aligned_moments) / (original_moments + 1e-10)
            max_change = np.max(relative_changes)
            
            # If max change is less than 0.001 (0.1%), rotation was likely ineffective
            if max_change < 0.001:
                return False
                
            # Check if moments are properly ordered (sorted)
            aligned_sorted = np.sort(aligned_moments)
            is_properly_ordered = np.allclose(aligned_moments, aligned_sorted, rtol=1e-6)
            
            return is_properly_ordered
            
        except Exception:
            return False
    
    def _setup_final_calculation(self) -> None:
        """Setup calculation with optimized geometry."""
        if not hasattr(self, 'optimized_geometry') or self.optimized_geometry is None:
            logger.info("No geometry optimization performed, using original geometry")
            return
        
        # Recreate molecular object with optimized geometry
        optimized_mol = self._create_optimized_molecule()
        
        # Recreate mean field object with optimized geometry
        self.mf = self._create_scf_method(optimized_mol)
        self.mf = self._apply_solvent_effects(self.mf)
        self._apply_calculation_settings()
        
        logger.info("Final calculation setup completed with optimized geometry")
    
    def _create_optimized_molecule(self):
        """Create molecular object with optimized geometry."""
        from pyscf import gto
        
        # Get atom symbols from original molecule
        atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
        
        # Create atom string with optimized coordinates
        atom_lines = []
        for symbol, coords in zip(atom_symbols, self.optimized_geometry):
            atom_lines.append(f"{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        atom_string = "\n".join(atom_lines)
        
        # Create new molecular object
        optimized_mol = gto.M(
            atom=atom_string,
            basis=self.mol.basis,
            charge=self.mol.charge,
            spin=self.mol.spin,
            verbose=0
        )
        optimized_mol.max_memory = self.mol.max_memory
        
        return optimized_mol
    
    def _apply_calculation_settings(self) -> None:
        """Apply common calculation settings."""
        self.mf.chkfile = self.get_checkpoint_path()
        if hasattr(self, 'max_cycle'):
            self.mf.max_cycle = self.max_cycle
        elif 'max_cycle' in self.results:
            self.mf.max_cycle = self.results['max_cycle']
    
    def _run_base_scf_calculation(self) -> float:
        """Run base SCF calculation and return energy."""
        logger.info(f"Running {self._get_base_method_description()} calculation...")
        energy = self.mf.kernel()
        logger.info(f"{self._get_base_method_description()} calculation completed")
        return energy
    
    def _verify_scf_convergence(self) -> None:
        """Verify SCF convergence and orbital data."""
        from .exceptions import ConvergenceError
        
        if not self.mf.converged:
            raise ConvergenceError(f"{self._get_base_method_description()} calculation failed to converge")
        
        if self.mf.mo_occ is None or len(self.mf.mo_occ) == 0:
            raise CalculationError(f"{self._get_base_method_description()} calculation failed: mo_occ not properly assigned")
        
        logger.info(f"Number of occupied orbitals: {self._count_occupied_orbitals()}")
    
    def _perform_orbital_analysis(self) -> None:
        """Perform orbital analysis and store results."""
        logger.info("Performing orbital analysis...")
        homo_idx, lumo_idx = self._analyze_orbitals()
        
        self.results.update({
            'homo_index': homo_idx,
            'lumo_index': lumo_idx,
            'num_occupied_orbitals': int(self._count_occupied_orbitals()),
            'num_virtual_orbitals': int(self._count_virtual_orbitals())
        })
        logger.info("Orbital analysis completed")
    
    def _perform_mulliken_analysis(self) -> None:
        """Perform Mulliken population analysis."""
        logger.info("Performing Mulliken population analysis...")
        mulliken_charges = self._calculate_mulliken_charges()
        
        self.results['mulliken_charges'] = mulliken_charges
        
        if mulliken_charges is not None:
            logger.info(f"Calculated Mulliken charges for {len(mulliken_charges)} atoms")
        else:
            logger.warning("Mulliken population analysis failed or was skipped")
    
    def _perform_frequency_analysis(self) -> None:
        """Perform vibrational frequency analysis."""
        logger.info("Performing vibrational frequency analysis...")
        frequency_results = self._calculate_frequencies()
        
        if frequency_results is not None:
            self.results.update(frequency_results)
            logger.info("Vibrational frequency analysis completed successfully")
        else:
            logger.warning("Vibrational frequency analysis failed or was skipped")
    
    def _prepare_final_results(self, specific_results: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare final results dictionary."""
        # Update with calculation-specific results
        self.results.update(specific_results)
        
        # Add common final results
        chk_path = self.get_checkpoint_path()
        self.results.update({
            'converged': True,
            'checkpoint_file': chk_path,
            'checkpoint_exists': os.path.exists(chk_path),
            'working_directory': self.working_dir,
            'optimized_geometry': self._geometry_to_xyz_string()
        })
        
        # Save files if requested
        if hasattr(self, 'keep_files') and self.keep_files:
            if hasattr(self, 'file_manager'):
                self.file_manager.save_calculation_results(self.working_dir, self.results)
                self.file_manager.save_geometry(self.working_dir, self.results['optimized_geometry'])
                logger.info(f"Calculation files saved to: {self.working_dir}")
        
        return self.results
    
    # ===== Abstract Methods for Subclasses =====
    
    @abstractmethod
    def _perform_specific_calculation(self, base_energy: float) -> Dict[str, Any]:
        """Perform calculation-specific computations."""
        pass
    
    @abstractmethod
    def _create_scf_method(self, mol):
        """Create appropriate SCF method object (RKS/UKS, RHF/UHF, etc.)."""
        pass
    
    @abstractmethod
    def _apply_solvent_effects(self, mf):
        """Apply solvent effects to mean field object."""
        pass
    
    @abstractmethod
    def _get_base_method_description(self) -> str:
        """Get description of base method for logging."""
        pass
    
    # ===== Default Implementations for Optional Methods =====
    
    def _requires_geometry_optimization(self) -> bool:
        """Whether this calculation requires geometry optimization."""
        return True
    
    def _requires_orbital_analysis(self) -> bool:
        """Whether this calculation requires orbital analysis."""
        return True
    
    def _requires_frequency_analysis(self) -> bool:
        """Whether this calculation requires frequency analysis."""
        return True