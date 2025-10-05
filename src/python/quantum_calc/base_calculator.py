"""Base calculator class for quantum chemistry calculations."""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Tuple
import os
import tempfile
import logging
from contextlib import contextmanager
import numpy as np
from .config_manager import get_memory_for_method, get_max_cycle

logger = logging.getLogger(__name__)


class BaseCalculator(ABC):
    """Abstract base class for quantum chemistry calculations."""
    
    def __init__(self, working_dir: Optional[str] = None, optimize_geometry: bool = True):
        """Initialize calculator with optional working directory and geometry optimization flag."""
        self.working_dir = working_dir or tempfile.mkdtemp(prefix="pyscf_calc_")
        self.results: Dict[str, Any] = {}
        self.optimize_geometry = optimize_geometry
        
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
            # Use calculation-specific default memory settings from config
            calculation_method = getattr(self, 'calculation_method', 'DFT')
            default_memory = get_memory_for_method(calculation_method)
            mol.max_memory = default_memory
            print(f"Using config-based PySCF max_memory: {default_memory} MB ({calculation_method})")
        
        # CPU cores are now configured at the process level in process_manager.py
        # This avoids conflicts and ensures proper timing of environment variable setup

    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """
        Template method for setting up quantum chemistry calculations.
        
        This method implements a common setup workflow:
        1. Extract and validate common parameters
        2. Validate calculation-specific parameters  
        3. Create molecular object
        4. Setup SCF method
        5. Store calculation parameters
        """
        try:
            logger.info(f"Starting {self._get_calculation_method_name()} calculation setup...")
            logger.info(f"Received parameters: {list(kwargs.keys())}")
            
            # Step 1: Extract and validate common parameters
            common_params = self._extract_common_parameters(**kwargs)
            
            # Step 2: Validate calculation-specific parameters
            specific_params = self._validate_specific_parameters(**kwargs)
            
            # Step 3: Create molecular object
            self._create_molecule_object(atoms, common_params)
            
            # Step 4: Setup SCF method
            self._setup_scf_method(common_params)
            
            # Step 5: Store calculation parameters
            self._store_calculation_parameters(common_params, specific_params, len(atoms))
            
            logger.info(f"{self._get_calculation_method_name()} setup completed successfully")
            
        except Exception as e:
            from .exceptions import InputError
            raise InputError(f"Failed to setup {self._get_calculation_method_name()} calculation: {str(e)}")
    
    def _extract_common_parameters(self, **kwargs) -> Dict[str, Any]:
        """Extract and validate common calculation parameters."""
        return {
            'basis': kwargs.get('basis', '6-31G(d)'),
            'charge': kwargs.get('charge', 0),
            'spin': kwargs.get('spin', 0),
            'max_cycle': kwargs.get('max_cycle', get_max_cycle()),
            'solvent_method': kwargs.get('solvent_method', 'none'),
            'solvent': kwargs.get('solvent', '-'),
            'memory_mb': kwargs.get('memory_mb', self._get_default_memory_mb())
        }
    
    def _validate_specific_parameters(self, **kwargs) -> Dict[str, Any]:
        """
        Validate calculation-specific parameters. 
        Subclasses should override this method for their specific validation needs.
        """
        return {}
    
    def _create_molecule_object(self, atoms: List[List], common_params: Dict[str, Any]) -> None:
        """Create PySCF molecular object with common parameters."""
        from pyscf import gto
        
        # Convert atoms list to PySCF format
        atom_string = self._atoms_to_string(atoms)
        
        # Create molecular object
        logger.info(f"Creating PySCF molecular object with {len(atoms)} atoms, "
                   f"basis={common_params['basis']}, charge={common_params['charge']}, "
                   f"spin={common_params['spin']}")
        
        self.mol = gto.M(
            atom=atom_string,
            basis=common_params['basis'],
            charge=common_params['charge'],
            spin=common_params['spin'],
            verbose=0
        )
        
        # Apply memory settings
        memory_mb = common_params['memory_mb']
        if memory_mb and memory_mb > 0:
            self.mol.max_memory = memory_mb
            logger.info(f"Set PySCF max_memory to {memory_mb} MB")
        else:
            default_memory = self._get_default_memory_mb()
            self.mol.max_memory = default_memory
            logger.info(f"Using default PySCF max_memory: {default_memory} MB")
        
        logger.info("PySCF molecular object created successfully")
    
    def _setup_scf_method(self, common_params: Dict[str, Any]) -> None:
        """Setup SCF method with common parameters."""
        # Store common parameters for template method access (before using them)
        self.max_cycle = common_params['max_cycle']
        self.solvent_method = common_params['solvent_method']
        self.solvent = common_params['solvent']
        
        # Store spin in results so _create_scf_method can access it
        self.results['spin'] = common_params['spin']
        self.results['charge'] = common_params['charge']
        self.results['basis'] = common_params['basis']

        # Create SCF method object (RHF/UHF, RKS/UKS, etc.)
        self.mf = self._create_scf_method(self.mol)

        # Apply solvent effects (now that solvent parameters are available)
        self.mf = self._apply_solvent_effects(self.mf)

        # Set common SCF parameters
        self.mf.chkfile = self.get_checkpoint_path()
        self.mf.max_cycle = common_params['max_cycle']

        logger.info(f"SCF method setup completed: {self._get_base_method_description()}")
    
    def _store_calculation_parameters(self, common_params: Dict[str, Any], 
                                    specific_params: Dict[str, Any], atom_count: int) -> None:
        """Store calculation parameters in results dictionary."""
        # Store common parameters
        self.results.update({
            'basis': common_params['basis'],
            'charge': common_params['charge'],
            'spin': common_params['spin'],
            'max_cycle': common_params['max_cycle'],
            'solvent_method': common_params['solvent_method'],
            'solvent': common_params['solvent'],
            'atom_count': atom_count,
            'method': self._get_method_description()
        })
        
        # Store calculation-specific parameters
        self.results.update(specific_params)
    
    def _get_default_memory_mb(self) -> int:
        """
        Get default memory setting for this calculation method.
        Uses configuration file or falls back to hardcoded values.
        """
        calculation_method = getattr(self, 'calculation_method', 'default')
        return get_memory_for_method(calculation_method)
    
    def _get_calculation_method_name(self) -> str:
        """
        Get the name of the calculation method for logging.
        Subclasses should override this method.
        """
        return getattr(self, 'calculation_method', 'Unknown')
    
    def _get_method_description(self) -> str:
        """
        Get the method description for results storage.
        Subclasses should override this method.
        """
        spin = self.results.get('spin', 0)
        base_method = self._get_base_method_description()
        return f"{base_method} ({self._get_calculation_method_name()})" if base_method else self._get_calculation_method_name()
    
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
        return self.optimize_geometry
    
    def _requires_orbital_analysis(self) -> bool:
        """Whether this calculation requires orbital analysis."""
        return True
    
    def _requires_frequency_analysis(self) -> bool:
        """Whether this calculation requires frequency analysis."""
        return True
    
    # ===== Common Validation and Helper Methods =====
    
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
    
    # ===== Common CASCI/CASSCF Analysis Methods =====
    
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
        
        try:
            # Get overlap matrix
            if hasattr(self.mol, 'get_ovlp'):
                S = self.mol.get_ovlp()
            else:
                S = self.mf.get_ovlp()
            
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