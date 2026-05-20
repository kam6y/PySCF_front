"""Base calculator class for quantum chemistry calculations."""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
import os
import tempfile
import logging
import sys
import re
import shutil
import subprocess
from importlib import util
import numpy as np
from .config_manager import get_memory_for_method, get_max_cycle
from .exceptions import CalculationError, PauseRequestedException
from ._checkpoint_resume import CheckpointResumeMixin
from ._frequency_analysis import FrequencyAnalysisMixin
from ._geometry_optimization import GeometryOptimizationMixin
from ._multiconfiguration import MulticonfigurationAnalysisMixin
from ._common_properties import CommonPropertiesMixin

logger = logging.getLogger(__name__)


class BaseCalculator(
    GeometryOptimizationMixin,
    FrequencyAnalysisMixin,
    CheckpointResumeMixin,
    MulticonfigurationAnalysisMixin,
    CommonPropertiesMixin,
    ABC
):
    """Abstract base class for quantum chemistry calculations."""
    
    def __init__(
        self,
        working_dir: Optional[str] = None,
        optimize_geometry: bool = True,
        geomopt_maxsteps: Optional[int] = None,
        geomopt_conv_energy: Optional[float] = None
    ):
        """Initialize calculator with optional working directory and geometry optimization flag."""
        self.working_dir = working_dir or tempfile.mkdtemp(prefix="pyscf_calc_")
        self.results: Dict[str, Any] = {}
        self.optimize_geometry = optimize_geometry
        self.geomopt_maxsteps = geomopt_maxsteps
        self.geomopt_conv_energy = geomopt_conv_energy
        self.gpu_enabled = False
        self._gpu4pyscf_available: Optional[bool] = None
        self._cuda_supported: Optional[bool] = None

    def _is_gpu_acceleration_enabled(self) -> bool:
        """Check whether GPU acceleration is enabled in app settings."""
        try:
            from .settings_manager import get_current_settings
            settings = get_current_settings()
            return bool(getattr(settings, "gpu_acceleration_enabled", False))
        except Exception:
            return False

    def _is_gpu4pyscf_available(self) -> bool:
        """Check whether GPU4PySCF is available on this system."""
        if not self._is_gpu_acceleration_enabled():
            return False
        if self._gpu4pyscf_available is None:
            if not sys.platform.startswith("linux"):
                self._gpu4pyscf_available = False
            else:
                self._gpu4pyscf_available = (
                    self._is_cuda_supported() and util.find_spec("gpu4pyscf") is not None
                )
        return self._gpu4pyscf_available

    def _is_cuda_supported(self) -> bool:
        """Check whether a supported CUDA Toolkit is available (via nvcc)."""
        if self._cuda_supported is None:
            if not sys.platform.startswith("linux"):
                self._cuda_supported = False
                return self._cuda_supported
            if shutil.which("nvcc") is None:
                self._cuda_supported = False
                return self._cuda_supported
            try:
                result = subprocess.run(
                    ["nvcc", "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
            except subprocess.TimeoutExpired:
                self._cuda_supported = False
                return self._cuda_supported
            if result.returncode != 0:
                self._cuda_supported = False
                return self._cuda_supported
            output = (result.stdout or "") + "\n" + (result.stderr or "")
            match = re.search(r"release\s+(\d+)\.(\d+)", output) or re.search(r"V(\d+)\.(\d+)", output)
            if not match:
                self._cuda_supported = False
                return self._cuda_supported
            major = int(match.group(1))
            self._cuda_supported = major in {11, 12, 13}
        return self._cuda_supported

    def _to_numpy(self, value: Any) -> Any:
        """Convert cupy arrays to numpy arrays if needed."""
        if value is None:
            return None
        try:
            import cupy as cp
        except Exception:
            return value
        if isinstance(value, cp.ndarray):
            return cp.asnumpy(value)
        if isinstance(value, (list, tuple)):
            return type(value)(self._to_numpy(item) for item in value)
        return value

    def _as_numpy_array(self, value: Any) -> Any:
        """Normalize array-like values to numpy arrays when possible."""
        value = self._to_numpy(value)
        if isinstance(value, (list, tuple)):
            return np.asarray(value)
        return value
        
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
            'memory_mb': kwargs.get('memory_mb', self._get_default_memory_mb()),
            'density_fitting': kwargs.get('density_fitting', False),
            'auxiliary_basis': kwargs.get('auxiliary_basis', None),
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

        # Apply density fitting if enabled
        self.density_fitting = common_params.get('density_fitting', False)
        self.auxiliary_basis = common_params.get('auxiliary_basis')
        if self.density_fitting:
            auxbasis = self.auxiliary_basis if self.auxiliary_basis else None
            self.mf = self.mf.density_fit(auxbasis=auxbasis)
            logger.info(f"Density fitting enabled with auxiliary basis: {auxbasis or 'auto'}")

        # Apply solvent effects (now that solvent parameters are available)
        self.mf = self._apply_solvent_effects(self.mf)

        # Set common SCF parameters
        self.mf.chkfile = self.get_checkpoint_path()
        self.mf.max_cycle = common_params['max_cycle']

        # Integrate SCF callback for pause support
        self.mf.callback = self._scf_callback

        logger.info(f"SCF method setup completed with pause callback: {self._get_base_method_description()}")
    
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
            'method': self._get_method_description(),
            'density_fitting': common_params.get('density_fitting', False),
            'auxiliary_basis': common_params.get('auxiliary_basis'),
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

            if self._requires_mulliken_analysis():
                self._perform_mulliken_analysis()

            if self._requires_frequency_analysis():
                self._perform_frequency_analysis()

            return self._prepare_final_results(specific_results)

        except PauseRequestedException:
            # Re-raise pause exception without modification
            # This will be caught by ProcessManager and handled appropriately
            logger.info("Calculation paused by user request")
            raise

        except Exception as e:
            import traceback

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
        if not keep_files and os.path.exists(self.working_dir) and "pyscf_calc_" in self.working_dir:
            shutil.rmtree(self.working_dir)
        elif keep_files:
            print(f"Calculation files preserved in: {self.working_dir}")

    # ===== Template Method Pattern Implementation =====
    
    def _pre_calculation_check(self) -> None:
        """Check prerequisites for calculation."""
        if not hasattr(self, 'mol') or self.mol is None:
            raise CalculationError("Molecular object not setup. Call setup_calculation first.")
        if not hasattr(self, 'mf') or self.mf is None:
            raise CalculationError("Mean field object not setup. Call setup_calculation first.")

    def _setup_final_calculation(self) -> None:
        """Setup calculation with optimized geometry."""
        if not hasattr(self, 'optimized_geometry') or self.optimized_geometry is None:
            logger.info("No geometry optimization performed, using original geometry")
            return
        
        # Recreate molecular object with optimized geometry
        optimized_mol = self._create_optimized_molecule()
        
        # Recreate mean field object with optimized geometry
        self.mf = self._create_scf_method(optimized_mol)
        if getattr(self, 'density_fitting', False):
            auxbasis = getattr(self, 'auxiliary_basis', None) or None
            self.mf = self.mf.density_fit(auxbasis=auxbasis)
        self.mf = self._apply_solvent_effects(self.mf)
        self._apply_calculation_settings()
        
        logger.info("Final calculation setup completed with optimized geometry")
    
    
    def _apply_calculation_settings(self) -> None:
        """Apply common calculation settings."""
        self.mf.chkfile = self.get_checkpoint_path()
        if hasattr(self, 'max_cycle'):
            self.mf.max_cycle = self.max_cycle
        elif 'max_cycle' in self.results:
            self.mf.max_cycle = self.results['max_cycle']

        # Integrate SCF callback for pause support
        self.mf.callback = self._scf_callback

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

    def _prepare_final_results(self, specific_results: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare final results dictionary."""
        # Update with calculation-specific results
        self.results.update(specific_results)

        resolved_auxiliary_basis = self._resolve_actual_auxiliary_basis()
        if resolved_auxiliary_basis is not None:
            self.results['resolved_auxiliary_basis'] = resolved_auxiliary_basis
        
        # Add common final results
        chk_path = self.get_checkpoint_path()
        self.results.update({
            'converged': True,
            'checkpoint_file': chk_path,
            'checkpoint_exists': os.path.exists(chk_path),
            'working_directory': self.working_dir,
            'optimized_geometry': self._geometry_to_xyz_string(),
            'gpu_enabled': bool(self.gpu_enabled)
        })
        
        # Save files if requested
        if hasattr(self, 'keep_files') and self.keep_files and hasattr(self, 'file_manager'):
            self.file_manager.save_calculation_results(self.working_dir, self.results)
            self.file_manager.save_geometry(self.working_dir, self.results['optimized_geometry'])
            logger.info(f"Calculation files saved to: {self.working_dir}")
        
        return self.results

    def _resolve_actual_auxiliary_basis(self) -> Optional[str]:
        """Resolve the auxiliary basis actually used by PySCF for density fitting."""
        if not getattr(self, 'density_fitting', False):
            return None

        mf = getattr(self, 'mf', None)
        if mf is None:
            return None

        with_df = getattr(mf, 'with_df', None)
        if with_df is None:
            return None

        auxbasis = getattr(with_df, 'auxbasis', None)
        if isinstance(auxbasis, str) and auxbasis.strip():
            return auxbasis

        auxmol = getattr(with_df, 'auxmol', None)
        auxmol_basis = getattr(auxmol, 'basis', None) if auxmol is not None else None
        return self._format_auxiliary_basis_value(auxmol_basis)

    def _format_auxiliary_basis_value(self, basis_value: Any) -> Optional[str]:
        """Format PySCF auxiliary basis metadata for result display."""
        if basis_value is None:
            return None

        if isinstance(basis_value, str):
            return basis_value

        if isinstance(basis_value, dict):
            formatted_entries = {
                str(atom): formatted
                for atom, value in basis_value.items()
                if (formatted := self._format_auxiliary_basis_value(value))
            }
            if not formatted_entries:
                return None

            unique_values = sorted(set(formatted_entries.values()))
            if len(unique_values) == 1:
                return unique_values[0]

            return ', '.join(
                f'{atom}: {formatted_entries[atom]}'
                for atom in sorted(formatted_entries)
            )

        return None
    
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

    def _requires_mulliken_analysis(self) -> bool:
        """Whether this calculation requires Mulliken population analysis."""
        return True
