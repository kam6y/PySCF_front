"""Worker runtime functions for process-based quantum chemistry calculations."""

import os
import sys
import logging

from .config_manager import get_memory_for_method
from .pause_manager import pause_manager
from .exceptions import PauseRequestedException

logger = logging.getLogger(__name__)


# ========== Worker Process Initialization ==========

def _worker_initializer():
    """
    Initialize worker process environment before any imports.

    This function is called by ProcessPoolExecutor when a new worker process starts,
    BEFORE the calculation_worker function is executed. It sets thread control
    environment variables to ensure that BLAS/LAPACK libraries initialize with
    the default thread count of 1.

    The actual thread count for each calculation is set later in _setup_worker_environment()
    based on the user's cpu_cores parameter.
    """
    import os

    thread_vars = [
        'OMP_NUM_THREADS',
        'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS',
        'BLIS_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS',
        'NUMEXPR_NUM_THREADS',
    ]

    for var in thread_vars:
        os.environ[var] = '1'


# ========== Private Helper Functions for calculation_worker ==========

def _setup_worker_environment(parameters: dict, process_logger) -> tuple:
    """
    Set up worker process environment variables and memory configuration.

    Thread Control Strategy (Unified Approach):
    This is the ONLY place where thread counts are configured for calculations.
    The configuration uses a three-layered approach:

    1. Environment variables (set here): Control BLAS/LAPACK libraries at the OS level
       - OMP_NUM_THREADS, MKL_NUM_THREADS, etc.
       - These must be set BEFORE importing PySCF to take effect

    2. PySCF lib.num_threads() (set in calculation_worker): Controls PySCF's internal threading
       - Set in calculation_worker function after PySCF import

    3. threadpoolctl context (used in calculation_worker): Runtime control for specific operations
       - Applied as a context manager around calculation execution
       - Provides additional safety layer

    Returns (cpu_cores, memory_mb).
    """
    # Get user-specified CPU cores or default to 1
    cpu_cores = parameters.get('cpu_cores') or 1
    cpu_cores_str = str(int(cpu_cores))

    # Set all parallel processing environment variables to control CPU usage
    # These environment variables must be set before any BLAS/LAPACK library is loaded
    os.environ['OMP_NUM_THREADS'] = cpu_cores_str
    os.environ['MKL_NUM_THREADS'] = cpu_cores_str
    os.environ['OPENBLAS_NUM_THREADS'] = cpu_cores_str
    os.environ['BLIS_NUM_THREADS'] = cpu_cores_str
    os.environ['VECLIB_MAXIMUM_THREADS'] = cpu_cores_str
    os.environ['NUMEXPR_NUM_THREADS'] = cpu_cores_str

    # Set appropriate memory defaults based on calculation method
    calculation_method = parameters.get('calculation_method', 'DFT')
    default_memory = get_memory_for_method(calculation_method)
    memory_mb = parameters.get('memory_mb') or default_memory

    # Log memory allocation
    if parameters.get('memory_mb'):
        process_logger.info(f"Using user-specified memory: {memory_mb} MB for {calculation_method}")
    else:
        process_logger.info(f"Using default memory: {memory_mb} MB for {calculation_method}")

    return cpu_cores, memory_mb


def _check_casci_dependencies(calculation_method: str, process_logger) -> None:
    """Check PySCF dependencies for CASCI/CASSCF calculations."""
    if calculation_method not in ['CASCI', 'CASSCF']:
        return

    process_logger.info("Performing PySCF dependency checks for CASCI/CASSCF...")
    try:
        import pyscf
        process_logger.info(f"PySCF version: {pyscf.__version__}")

        from pyscf import mcscf
        process_logger.info("PySCF mcscf module loaded successfully")

        from pyscf import gto
        test_mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', verbose=0)
        process_logger.info("PySCF basic functionality test passed")
    except ImportError as e:
        process_logger.error(f"PySCF dependency check failed: {e}")
        process_logger.error("CASCI/CASSCF calculations will likely fail")
    except Exception as e:
        process_logger.warning(f"PySCF functionality test encountered issues: {e}")
        process_logger.warning("CASCI/CASSCF calculations may have issues")


def _import_calculator_classes(process_logger):
    """
    Import calculator classes and return them as a dict.
    Returns (calculators_dict, exception_classes_tuple).
    """
    from quantum_calc import DFTCalculator, HFCalculator, MP2Calculator, CCSDCalculator, TDDFTCalculator
    from quantum_calc import CalculationError, ConvergenceError, InputError, PauseRequestedException

    calculators = {
        'DFT': DFTCalculator,
        'HF': HFCalculator,
        'MP2': MP2Calculator,
        'CCSD': CCSDCalculator,
        'CCSD_T': CCSDCalculator,
        'TDDFT': TDDFTCalculator,
        'CASCI': None,
        'CASSCF': None
    }

    # Try importing CASCI/CASSCF
    try:
        from quantum_calc import CASCICalculator, CASSCFCalculator
        calculators['CASCI'] = CASCICalculator
        calculators['CASSCF'] = CASSCFCalculator
        process_logger.info("Successfully imported CASCI/CASSCF calculators")
    except ImportError as e:
        process_logger.error(f"Failed to import CASCI/CASSCF calculators: {e}")
        process_logger.error("CASCI and CASSCF calculations will not be available")
    except Exception as e:
        process_logger.error(f"Unexpected error importing CASCI/CASSCF calculators: {e}")
        process_logger.error("CASCI and CASSCF calculations will not be available")

    return calculators, (CalculationError, ConvergenceError, InputError, PauseRequestedException)


def _create_calculator_instance(calculation_method: str, parameters: dict,
                                calc_dir: str, calculator_classes: dict, process_logger):
    """Create and return appropriate calculator instance."""
    optimize_geometry = parameters.get('optimize_geometry', True)
    molecule_name = parameters['name']

    calculator_class = calculator_classes.get(calculation_method)

    if calculator_class is None:
        if calculation_method in ['CASCI', 'CASSCF']:
            raise ImportError(
                f"{calculation_method} calculator is not available. "
                "Please check PySCF mcscf module installation."
            )
        # Default to DFT if unknown method
        calculator_class = calculator_classes['DFT']
        process_logger.warning(f"Unknown calculation method '{calculation_method}', defaulting to DFT")

    return calculator_class(
        working_dir=calc_dir,
        keep_files=True,
        molecule_name=molecule_name,
        optimize_geometry=optimize_geometry
    )


def _prepare_setup_parameters(parameters: dict, memory_mb: int) -> dict:
    """Prepare setup parameters dict for calculator."""
    calculation_method = parameters.get('calculation_method', 'DFT')

    setup_params = {
        'basis': parameters['basis_function'],
        'charge': parameters['charges'],
        'spin': parameters['spin'],
        'max_cycle': 150,
        'solvent_method': parameters['solvent_method'],
        'solvent': parameters['solvent'],
        'memory_mb': memory_mb,
    }

    # Add exchange-correlation functional for DFT and TDDFT
    if calculation_method in ['DFT', 'TDDFT']:
        setup_params['xc'] = parameters['exchange_correlation']

    # Add CCSD-specific parameters
    if calculation_method in ['CCSD', 'CCSD_T']:
        setup_params['frozen_core'] = parameters.get('frozen_core', True)
        setup_params['ccsd_t'] = (calculation_method == 'CCSD_T')

    # Add TDDFT-specific parameters
    if calculation_method == 'TDDFT':
        setup_params['nstates'] = parameters.get('tddft_nstates', 10)
        setup_params['tddft_method'] = parameters.get('tddft_method', 'TDDFT')
        setup_params['analyze_nto'] = parameters.get('tddft_analyze_nto', False)

    # Add CASCI/CASSCF-specific parameters
    if calculation_method in ['CASCI', 'CASSCF']:
        setup_params['ncas'] = parameters.get('ncas', 4)
        setup_params['nelecas'] = parameters.get('nelecas', 4)
        setup_params['natorb'] = parameters.get('natorb', True)
        setup_params['max_cycle_micro'] = parameters.get('max_cycle_micro', 4)

        if calculation_method == 'CASSCF':
            setup_params['max_cycle_macro'] = parameters.get('max_cycle_macro', 50)
            setup_params['conv_tol'] = parameters.get('conv_tol', 1e-6)
            setup_params['conv_tol_grad'] = parameters.get('conv_tol_grad', 1e-4)

    return setup_params


def _handle_calculation_error(error: Exception, calc_dir: str, repository,
                              calculation_method: str, memory_mb: int,
                              cpu_cores: int, process_logger) -> tuple:
    """
    Handle calculation errors, save error information, and return error details.
    Returns (success=False, error_message).
    """
    from quantum_calc import CalculationError, ConvergenceError, InputError

    error_message = str(error)
    error_type = type(error).__name__

    # Build error diagnosis
    error_info = {
        'error_type': error_type,
        'error_message': error_message,
        'calculation_method': calculation_method,
        'memory_mb': memory_mb,
        'cpu_cores': cpu_cores
    }

    # Add specific diagnoses based on error type and content
    if isinstance(error, (InputError, ConvergenceError, CalculationError)):
        if calculation_method in ['CASCI', 'CASSCF']:
            if 'import' in error_message.lower() or 'mcscf' in error_message.lower():
                error_info['diagnosis'] = 'PySCF mcscf module import failure - check PySCF installation'
                error_info['suggestion'] = 'Install PySCF with: conda install pyscf -c pyscf'
            elif 'memory' in error_message.lower():
                error_info['diagnosis'] = 'Insufficient memory for CASCI/CASSCF calculation'
                error_info['suggestion'] = f'Increase memory allocation (current: {memory_mb} MB, try: {memory_mb * 2} MB)'
            elif 'active' in error_message.lower() and 'space' in error_message.lower():
                error_info['diagnosis'] = 'Invalid active space configuration'
                error_info['suggestion'] = 'Check ncas and nelecas parameters'
    elif isinstance(error, ImportError):
        error_info['diagnosis'] = 'Python module import failure'
        error_info['suggestion'] = 'Check PySCF installation and dependencies'
        if 'mcscf' in error_message.lower():
            error_info['diagnosis'] = 'PySCF mcscf module not found'
            error_info['suggestion'] = 'Install complete PySCF package: conda install pyscf -c pyscf'
    else:
        # General error diagnosis
        error_str = error_message.lower()
        if 'pyscf' in error_str:
            error_info['diagnosis'] = 'PySCF library error'
            error_info['suggestion'] = 'Check PySCF installation and system compatibility'
        elif 'memory' in error_str or 'malloc' in error_str:
            error_info['diagnosis'] = 'Memory allocation error'
            error_info['suggestion'] = f'Increase available system memory or reduce memory_mb (current: {memory_mb} MB)'
        elif 'thread' in error_str or 'lock' in error_str:
            error_info['diagnosis'] = 'Threading/concurrency error'
            error_info['suggestion'] = 'Check system threading configuration'
        else:
            error_info['diagnosis'] = 'Unexpected error during calculation'
            error_info['suggestion'] = 'Check logs for more details'

    process_logger.error(f"Calculation error: {error_message}")
    process_logger.error(f"Error diagnosis: {error_info}")

    # Save error status and information
    repository.save_calculation_status(calc_dir, 'error')
    repository.save_calculation_results(calc_dir, {'error': error_message, 'diagnosis': error_info})

    return False, error_message


# ========== End of Private Helper Functions ==========


def calculation_worker(calculation_id: str, parameters: dict) -> tuple:
    """
    Worker function to run quantum chemistry calculations in a separate process.
    Returns (success: bool, error_message: str or None)
    """
    # Setup logging for this process
    process_logger = logging.getLogger(f'worker_{calculation_id}')
    process_logger.setLevel(logging.INFO)

    # Initialize variables for finally block
    original_threads = None

    # Setup environment and get configuration
    cpu_cores, memory_mb = _setup_worker_environment(parameters, process_logger)
    calculation_method = parameters.get('calculation_method', 'DFT')

    # Perform dependency checks for CASCI/CASSCF
    _check_casci_dependencies(calculation_method, process_logger)

    # Import calculator classes and exception types
    from quantum_calc._calculation_repository import CalculationRepository
    from quantum_calc import get_current_settings
    from quantum_calc.pause_manager import pause_manager
    from threadpoolctl import threadpool_info, threadpool_limits
    from pyscf import lib

    calculator_classes, exception_types = _import_calculator_classes(process_logger)

    # Load current settings to get calculations directory
    settings = get_current_settings()
    repository = CalculationRepository(base_dir=settings.calculations_directory)
    calc_dir = os.path.join(repository.get_base_directory(), calculation_id)

    try:
        # Update status to running
        repository.save_calculation_status(calc_dir, 'running')
        process_logger.info(f"Starting calculation {calculation_id} in process {os.getpid()}")
        process_logger.info(f"Using {cpu_cores} CPU cores and {memory_mb} MB memory")

        # Log detected threadpool libraries for debugging
        try:
            thread_info = threadpool_info()
            process_logger.info(f"Detected threadpool libraries: {len(thread_info)} found")
            for info in thread_info:
                process_logger.info(
                    f"  {info.get('user_api', 'unknown')}: {info.get('internal_api', 'unknown')} "
                    f"- threads: {info.get('num_threads', 'unknown')}"
                )
        except Exception as e:
            process_logger.warning(f"Could not get threadpool info: {e}")

        # Set PySCF thread count
        try:
            original_threads = lib.num_threads()
            lib.num_threads(int(cpu_cores))
            new_threads = lib.num_threads()
            process_logger.info(f"PySCF threads: {original_threads} -> {new_threads} (requested: {cpu_cores})")
        except Exception as e:
            process_logger.warning(f"Could not set PySCF threads: {e}")

        # Create calculator instance
        calculator = _create_calculator_instance(
            calculation_method, parameters, calc_dir, calculator_classes, process_logger
        )

        # Parse XYZ and setup calculation
        atoms = calculator.parse_xyz(parameters['xyz'])
        setup_params = _prepare_setup_parameters(parameters, memory_mb)
        calculator.setup_calculation(atoms, **setup_params)

        # Resume from checkpoint if this is a resumed calculation
        # IMPORTANT: Must be called AFTER setup_calculation() so that self.mf exists
        if parameters.get('resume_from_pause', False):
            pause_state = parameters.get('pause_state')
            process_logger.info(f"Resuming calculation {calculation_id} from checkpoint")
            if pause_state:
                process_logger.info(f"Pause state: {pause_state}")
            calculator.resume_from_checkpoint(pause_state)

        # Run calculation with controlled BLAS/LAPACK/OpenMP threading
        # Note: Not specifying user_api controls ALL threadpool libraries (blas, openmp, etc.)
        process_logger.info(f"Executing calculation with threadpool_limits(limits={cpu_cores})")
        with threadpool_limits(limits=int(cpu_cores)):
            results = calculator.run_calculation()

        # Save results and update status to completed
        repository.save_calculation_results(calc_dir, results)
        repository.save_calculation_status(calc_dir, 'completed')

        # Clean up pause state file if it exists (from previous pause/resume cycle)
        repository.delete_pause_state(calc_dir)
        process_logger.debug(f"Cleaned up pause state file for calculation {calculation_id}")

        process_logger.info(f"Calculation {calculation_id} completed successfully in process {os.getpid()}")
        return True, None

    except Exception as e:
        # Check if this is a pause request - handle it specially
        from quantum_calc.exceptions import PauseRequestedException
        if isinstance(e, PauseRequestedException):
            process_logger.info(f"Calculation {calculation_id} was paused by user request")

            # Save pause state information
            pause_state = {
                'calculation_phase': 'scf_calculation',  # Default to SCF phase
                'checkpoint_exists': os.path.exists(os.path.join(calc_dir, 'calculation.chk')),
            }

            # Check if there's a geometry trajectory file to determine optimization step
            trajectory_file = os.path.join(calc_dir, 'geom_opt_trajectory.xyz')
            if os.path.exists(trajectory_file):
                try:
                    with open(trajectory_file, 'r') as f:
                        content = f.read()
                        # Count the number of geometry steps
                        step_count = content.count('Optimization step')
                        if step_count > 0:
                            pause_state['optimization_step'] = step_count
                            pause_state['calculation_phase'] = 'geometry_optimization'
                except Exception as read_error:
                    process_logger.warning(f"Failed to read geometry trajectory: {read_error}")

            # Save pause state to file
            repository.save_pause_state(calc_dir, pause_state)

            # Update status to 'paused'
            repository.save_calculation_status(calc_dir, 'paused')

            # Remove pause flag file
            pause_manager.remove_pause_flag_file(calc_dir)

            process_logger.info(f"Calculation {calculation_id} paused successfully, state saved")

            # Re-raise exception so the parent process knows it was paused
            raise e

        # For all other exceptions, handle as errors
        return _handle_calculation_error(e, calc_dir, repository, calculation_method, memory_mb, cpu_cores, process_logger)

    finally:
        # Clean up pause flag file if it exists
        # (This is idempotent - safe to call even if already removed)
        try:
            pause_manager.remove_pause_flag_file(calc_dir)
            process_logger.debug(f"Cleaned up pause flag file for calculation {calculation_id}")
        except Exception as cleanup_error:
            process_logger.warning(f"Failed to clean up pause flag file: {cleanup_error}")

        # Restore original PySCF thread count
        if original_threads is not None:
            try:
                lib.num_threads(original_threads)
                process_logger.info(f"Restored PySCF threads to {original_threads}")
            except Exception as e:
                process_logger.warning(f"Could not restore PySCF threads: {e}")
