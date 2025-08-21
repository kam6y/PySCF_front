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
            # デフォルトメモリ設定を適用
            mol.max_memory = 2000
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
    
    @abstractmethod
    def run_calculation(self) -> Dict[str, Any]:
        """Run the quantum chemistry calculation and return results."""
        pass
    
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