"""Process pool manager for CPU-bound quantum chemistry calculations."""

import os
import sys
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, Future
from typing import Dict, Any, Optional, Callable, List
import time
import threading
from datetime import datetime
from threadpoolctl import threadpool_limits
from queue import Queue
from dataclasses import dataclass
from .config_manager import get_memory_for_method
from .resource_manager import AllocationStatus
from .pause_manager import pause_manager
from .exceptions import PauseRequestedException

logger = logging.getLogger(__name__)


class ProcessManagerError(Exception):
    """Exception raised when process manager initialization or operation fails."""
    pass


@dataclass
class QueuedCalculation:
    """Represents a calculation waiting in the queue."""
    calculation_id: str
    parameters: dict
    created_at: datetime
    waiting_reason: Optional[str] = None
    
    def __lt__(self, other):
        """For priority queue ordering by creation time."""
        return self.created_at &lt; other.created_at


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


def _handle_calculation_error(error: Exception, calc_dir: str, file_manager,
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
    file_manager.save_calculation_status(calc_dir, 'error')
    file_manager.save_calculation_results(calc_dir, {'error': error_message, 'diagnosis': error_info})
    
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
    
    # Setup environment and get configuration
    cpu_cores, memory_mb = _setup_worker_environment(parameters, process_logger)
    calculation_method = parameters.get('calculation_method', 'DFT')
    
    # Perform dependency checks for CASCI/CASSCF
    _check_casci_dependencies(calculation_method, process_logger)
    
    # Import calculator classes and exception types
    from quantum_calc.file_manager import CalculationFileManager
    from quantum_calc import get_current_settings
    from threadpoolctl import threadpool_info, threadpool_limits
    from pyscf import lib

    calculator_classes, exception_types = _import_calculator_classes(process_logger)

    # Load current settings to get calculations directory
    settings = get_current_settings()
    file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
    
    try:
        # Update status to running
        file_manager.save_calculation_status(calc_dir, 'running')
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
        original_threads = None
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

        # Resume from checkpoint if this is a resumed calculation
        if parameters.get('resume_from_pause', False):
            pause_state = parameters.get('pause_state')
            process_logger.info(f"Resuming calculation {calculation_id} from checkpoint")
            if pause_state:
                process_logger.info(f"Pause state: {pause_state}")
            calculator.resume_from_checkpoint(pause_state)

        # Parse XYZ and setup calculation
        atoms = calculator.parse_xyz(parameters['xyz'])
        setup_params = _prepare_setup_parameters(parameters, memory_mb)
        calculator.setup_calculation(atoms, **setup_params)
        
        # Run calculation with controlled BLAS/LAPACK threading
        process_logger.info(f"Executing calculation with threadpool_limits(limits={cpu_cores}, user_api='blas')")
        with threadpool_limits(limits=int(cpu_cores), user_api='blas'):
            results = calculator.run_calculation()
        
        # Save results and update status to completed
        file_manager.save_calculation_results(calc_dir, results)
        file_manager.save_calculation_status(calc_dir, 'completed')

        # Clean up pause state file if it exists (from previous pause/resume cycle)
        file_manager.delete_pause_state(calc_dir)
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
            file_manager.save_pause_state(calc_dir, pause_state)

            # Update status to 'paused'
            file_manager.save_calculation_status(calc_dir, 'paused')

            # Remove pause flag file
            from quantum_calc.pause_manager import pause_manager
            pause_manager.remove_pause_flag_file(calc_dir)

            process_logger.info(f"Calculation {calculation_id} paused successfully, state saved")

            # Return success - the 'paused' status is saved in the file
            return True, None

        # For all other exceptions, handle as errors
        return _handle_calculation_error(e, calc_dir, file_manager, calculation_method, memory_mb, cpu_cores, process_logger)
    
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


class CalculationProcessManager:
    """Manages a process pool for quantum chemistry calculations with queueing support and resource management."""
    
    def __init__(self, max_workers: Optional[int] = None, max_parallel_instances: Optional[int] = None,
                 notification_callback: Optional[Callable] = None):
        """
        Initialize the process manager with a fixed process pool and queueing system.
        
        Args:
            max_workers: Maximum number of worker processes. If None, defaults to CPU count.
            max_parallel_instances: Maximum number of parallel calculation instances. If None, defaults to max_workers.
            notification_callback: Optional callback function for sending WebSocket notifications.
                                   Signature: callback(calculation_id, status, error_message)
        """
        self.notification_callback = notification_callback
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.max_parallel_instances = max_parallel_instances or self.max_workers
        self.active_futures: Dict[str, Future] = {}
        self.completion_callbacks: Dict[str, List[Callable]] = {}  # calculation_id -> list of callbacks
        self.calculation_queue: List[QueuedCalculation] = []  # Waiting calculations queue
        self._shutdown = False
        
        # Thread safety for queue processing
        self._queue_lock = threading.Lock()
        self._queue_processing = False
        
        # Initialize resource manager with error handling
        try:
            from .resource_manager import get_resource_manager
            self.resource_manager = get_resource_manager()
            logger.info("Resource manager initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize resource manager: {e}")
            # Create a minimal fallback resource manager
            self.resource_manager = None
            logger.warning("Continuing without resource manager - calculations may not be resource-constrained")
        
        # Initialize process pool executor with error handling
        self.executor = None
        try:
            self.executor = ProcessPoolExecutor(
                max_workers=self.max_workers,
                mp_context=multiprocessing.get_context('spawn')  # Better cross-platform compatibility
            )
            logger.info(f"ProcessPoolExecutor created successfully with {self.max_workers} workers")
        except Exception as e:
            logger.error(f"Failed to create ProcessPoolExecutor: {e}")
            raise ProcessManagerError(f"Cannot create process pool: {str(e)}")
        
        # Background resource monitoring (only if resource manager is available)
        self._resource_monitor_thread = None
        self._resource_monitor_stop_event = threading.Event()
        self._resource_monitor_interval = 10.0  # Check every 10 seconds
        
        if self.resource_manager is not None:
            try:
                self._start_resource_monitoring()
                logger.info("Background resource monitoring started")
            except Exception as e:
                logger.warning(f"Failed to start resource monitoring: {e}")
        else:
            logger.warning("Skipping resource monitoring due to resource manager initialization failure")
        
        logger.info(f"Initialized CalculationProcessManager with {self.max_workers} worker processes and {self.max_parallel_instances} max parallel instances")
    
    def _start_resource_monitoring(self):
        """Start the background resource monitoring thread."""
        if self._resource_monitor_thread is None or not self._resource_monitor_thread.is_alive():
            self._resource_monitor_stop_event.clear()
            self._resource_monitor_thread = threading.Thread(
                target=self._resource_monitoring_loop,
                name="ResourceMonitor",
                daemon=True
            )
            self._resource_monitor_thread.start()
            logger.info("Started background resource monitoring")
    
    def _stop_resource_monitoring(self):
        """Stop the background resource monitoring thread."""
        if self._resource_monitor_thread and self._resource_monitor_thread.is_alive():
            logger.info("Stopping background resource monitoring")
            self._resource_monitor_stop_event.set()
            self._resource_monitor_thread.join(timeout=5.0)
            if self._resource_monitor_thread.is_alive():
                logger.warning("Resource monitoring thread did not stop cleanly")
    
    def _resource_monitoring_loop(self):
        """Background loop that monitors system resources and processes queue when resources improve."""
        logger.info("Resource monitoring loop started")
        
        while not self._resource_monitor_stop_event.is_set():
            try:
                # Only process queue if there are queued calculations
                if self.calculation_queue and not self._shutdown:
                    # Check if resources have improved (if resource manager is available)
                    should_process = False
                    
                    if self.resource_manager is not None:
                        try:
                            has_improved, reason = self.resource_manager.has_resources_improved()
                            if has_improved:
                                logger.info(f"Resources improved: {reason}. Processing queue...")
                                should_process = True
                        except Exception as e:
                            logger.warning(f"Failed to check resource improvements: {e}. Processing queue anyway.")
                            should_process = True
                    else:
                        # No resource manager - periodically try to process queue
                        should_process = True
                    
                    if should_process:
                        self._process_queue()
                
                # Wait for the next check or until stop event is set
                self._resource_monitor_stop_event.wait(timeout=self._resource_monitor_interval)
                
            except Exception as e:
                logger.error(f"Error in resource monitoring loop: {e}")
                # Continue monitoring despite errors
                self._resource_monitor_stop_event.wait(timeout=self._resource_monitor_interval)
        
        logger.info("Resource monitoring loop stopped")
    

    def set_max_parallel_instances(self, max_instances: int) -> None:
        """
        Update the maximum number of parallel instances.
        
        Args:
            max_instances: New maximum number of parallel instances.
        """
        old_max = self.max_parallel_instances
        self.max_parallel_instances = max(1, min(max_instances, self.max_workers))
        logger.info(f"Updated max parallel instances from {old_max} to {self.max_parallel_instances}")
        
        # If we increased the limit, try to start queued calculations
        if self.max_parallel_instances > old_max:
            logger.info(f"Max parallel instances increased, processing queue (current queue size: {len(self.calculation_queue)})")
            self._process_queue()
    
    def submit_calculation(self, calculation_id: str, parameters: dict) -> tuple[bool, str, Optional[str]]:
        """
        Submit a calculation to the process pool or queue with resource checking.

        Args:
            calculation_id: Unique identifier for the calculation
            parameters: Calculation parameters

        Returns:
            Tuple of (success: bool, status: str, waiting_reason: Optional[str])
            where status is 'running', 'waiting', or 'error'
        """
        if self._shutdown:
            logger.error("Cannot submit calculation: process manager is shut down")
            return False, 'error', None

        if self.executor is None:
            logger.error("Process pool executor is not available")
            return False, 'error', None
        
        # Extract calculation parameters
        user_cpu_cores = parameters.get('cpu_cores') or 1
        user_memory_mb = parameters.get('memory_mb') or get_memory_for_method(parameters.get('calculation_method', 'DFT'))
        calculation_method = parameters.get('calculation_method', 'DFT')
        
        # Unified resource allocation check (skip if resource manager is not available)
        allocation_status = AllocationStatus.CAN_START
        reason = "Resource checking disabled"
        
        if self.resource_manager is not None:
            try:
                allocation_status, reason = self.resource_manager.check_allocation_status(
                    cpu_cores=user_cpu_cores,
                    memory_mb=user_memory_mb,
                    calculation_method=calculation_method,
                    active_count=len(self.active_futures),
                    max_slots=self.max_parallel_instances
                )
            except Exception as e:
                logger.warning(f"Resource allocation check failed: {e}. Proceeding with CAN_START.")
                allocation_status = AllocationStatus.CAN_START
                reason = "Resource checking failed - proceeding without constraints"
        
        # Handle allocation status
        if allocation_status == AllocationStatus.INSUFFICIENT_RESOURCES:
            # System resources are fundamentally insufficient - return error
            logger.error(f"System resources insufficient for calculation {calculation_id}: {reason}")
            return False, 'error', reason
        
        elif allocation_status == AllocationStatus.SHOULD_QUEUE:
            # Add to queue - resources temporarily unavailable or slots full
            queued_calc = QueuedCalculation(
                calculation_id=calculation_id,
                parameters=parameters,
                created_at=datetime.now(),
                waiting_reason=reason
            )
            self.calculation_queue.append(queued_calc)
            self.calculation_queue.sort(key=lambda x: x.created_at)

            logger.info(f"Added calculation {calculation_id} to queue (position {len(self.calculation_queue)}): {reason}")

            # Send WebSocket notification for waiting status
            self._send_websocket_notification(calculation_id, 'waiting', None)

            return True, 'waiting', reason
        
        # allocation_status == AllocationStatus.CAN_START - start immediately
        try:
            logger.info(f"About to submit calculation {calculation_id} to executor")
            
            # Register resources with the resource manager (if available)
            if self.resource_manager is not None:
                try:
                    self.resource_manager.register_calculation(
                        calculation_id=calculation_id,
                        cpu_cores=user_cpu_cores,
                        memory_mb=user_memory_mb,
                        calculation_method=calculation_method
                    )
                except Exception as e:
                    logger.warning(f"Failed to register calculation resources: {e}")
                    # Continue without resource registration
            
            # Submit to executor
            future = self.executor.submit(calculation_worker, calculation_id, parameters)
            self.active_futures[calculation_id] = future
            
            # Add callback to clean up completed futures and process queue
            future.add_done_callback(lambda f: self._cleanup_future(calculation_id, f))

            logger.info(f"Started calculation {calculation_id} immediately with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(self.active_futures)}/{self.max_parallel_instances} slots used)")

            # Send WebSocket notification for running status
            self._send_websocket_notification(calculation_id, 'running', None)

            return True, 'running', None
            
        except Exception as e:
            logger.error(f"Failed to submit calculation {calculation_id}: {e}")
            import traceback
            logger.error(f"Full traceback:\n{traceback.format_exc()}")
            
            # Unregister resources on failure (if resource manager is available)
            if self.resource_manager is not None:
                try:
                    self.resource_manager.unregister_calculation(calculation_id)
                except Exception as cleanup_error:
                    logger.warning(f"Failed to unregister calculation resources on failure: {cleanup_error}")
            return False, 'error', None
    
    def _cleanup_future(self, calculation_id: str, future: Future):
        """Clean up a completed future, call completion callbacks, and process queue."""
        success = False
        error_message = None
        
        try:
            if calculation_id in self.active_futures:
                del self.active_futures[calculation_id]
            
            # Unregister resources from resource manager (if available)
            if self.resource_manager is not None:
                try:
                    self.resource_manager.unregister_calculation(calculation_id)
                except Exception as e:
                    logger.warning(f"Failed to unregister calculation resources: {e}")
            
            # Determine result and log
            if future.exception():
                exception = future.exception()
                # Check if calculation was paused
                if isinstance(exception, PauseRequestedException):
                    logger.info(f"Calculation {calculation_id} was paused by user request")
                    # Get file manager to update status
                    from quantum_calc.file_manager import CalculationFileManager
                    from quantum_calc import get_current_settings
                    settings = get_current_settings()
                    file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
                    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)

                    # Update status to 'paused'
                    file_manager.save_calculation_status(calc_dir, 'paused')

                    # Clear pause request and remove flag file
                    pause_manager.clear_pause_request(calculation_id)
                    pause_manager.remove_pause_flag_file(calc_dir)

                    # Send WebSocket notification
                    self._send_websocket_notification(calculation_id, 'paused', None)
                else:
                    error_message = str(exception)
                    logger.error(f"Calculation {calculation_id} failed with exception: {error_message}")
                    # Send WebSocket notification for exception error
                    self._send_websocket_notification(calculation_id, 'error', error_message)
            else:
                success, calc_error = future.result()
                if success:
                    logger.info(f"Calculation {calculation_id} completed successfully")
                    # Send WebSocket notification for successful completion
                    self._send_websocket_notification(calculation_id, 'completed', None)
                else:
                    error_message = calc_error
                    logger.warning(f"Calculation {calculation_id} failed: {error_message}")
                    # Send WebSocket notification for calculation error
                    self._send_websocket_notification(calculation_id, 'error', error_message)
            
            # Call completion callbacks
            if calculation_id in self.completion_callbacks:
                callbacks = self.completion_callbacks.pop(calculation_id)
                for callback in callbacks:
                    try:
                        callback(calculation_id, success, error_message)
                    except Exception as callback_error:
                        logger.error(f"Error in completion callback for {calculation_id}: {callback_error}")
            
            # Process queue to start next calculation if available
            logger.info(f"Calculation {calculation_id} cleanup completed. Processing queue for waiting calculations...")
            self._process_queue()
                        
        except Exception as e:
            logger.error(f"Error cleaning up future for {calculation_id}: {e}")
    
    def _process_queue(self):
        """Process the calculation queue and start next calculations if slots and resources are available."""
        # Thread-safe queue processing with detailed logging
        with self._queue_lock:
            if self._queue_processing:
                logger.debug("Queue processing already in progress, skipping")
                return
            self._queue_processing = True
            
        try:
            logger.info(f"Processing calculation queue - Active: {len(self.active_futures)}/{self.max_parallel_instances}, Queued: {len(self.calculation_queue)}")
            started_count = 0
            
            # Process queue until no more calculations can start
            while self.calculation_queue and not self._shutdown:
                # Try to start the next calculation from queue
                started = False
                
                # Iterate through queue to find a calculation that can start
                # (some may need to wait for specific resources while others can start)
                for i, queued_calc in enumerate(self.calculation_queue):
                    user_cpu_cores = queued_calc.parameters.get('cpu_cores') or 1
                    user_memory_mb = queued_calc.parameters.get('memory_mb') or get_memory_for_method(queued_calc.parameters.get('calculation_method', 'DFT'))
                    calculation_method = queued_calc.parameters.get('calculation_method', 'DFT')
                    
                    # Unified resource allocation check (skip if resource manager is not available)
                    allocation_status = AllocationStatus.CAN_START
                    reason = "Resource checking disabled"
                    
                    if self.resource_manager is not None:
                        try:
                            allocation_status, reason = self.resource_manager.check_allocation_status(
                                cpu_cores=user_cpu_cores,
                                memory_mb=user_memory_mb,
                                calculation_method=calculation_method,
                                active_count=len(self.active_futures),
                                max_slots=self.max_parallel_instances
                            )
                        except Exception as e:
                            logger.warning(f"Resource allocation check failed for queued calculation: {e}")
                            allocation_status = AllocationStatus.CAN_START
                            reason = "Resource checking failed"
                    
                    # Handle allocation status
                    if allocation_status == AllocationStatus.INSUFFICIENT_RESOURCES:
                        # System resources fundamentally insufficient - remove from queue and mark as error
                        next_calc = self.calculation_queue.pop(i)
                        logger.warning(f"Removing calculation {next_calc.calculation_id} from queue due to insufficient resources: {reason}")
                        self._update_calculation_status_from_queue(next_calc.calculation_id, 'error', reason)
                        started = True  # Continue processing queue
                        break
                    
                    elif allocation_status == AllocationStatus.CAN_START:
                        # Can start this calculation - remove from queue and start it
                        next_calc = self.calculation_queue.pop(i)
                        
                        try:
                            # Register resources with the resource manager (if available)
                            if self.resource_manager is not None:
                                try:
                                    self.resource_manager.register_calculation(
                                        calculation_id=next_calc.calculation_id,
                                        cpu_cores=user_cpu_cores,
                                        memory_mb=user_memory_mb,
                                        calculation_method=calculation_method
                                    )
                                except Exception as e:
                                    logger.warning(f"Failed to register queued calculation resources: {e}")
                                    # Continue without resource registration
                            
                            # Start the calculation
                            # Check if executor has been shut down before submitting
                            if self._shutdown or self.executor._shutdown:
                                logger.warning(f"Executor has been shut down. Cannot schedule calculation {next_calc.calculation_id}")
                                self._update_calculation_status_from_queue(next_calc.calculation_id, 'error', 'Process manager has been shut down')
                                started = True
                                break

                            future = self.executor.submit(calculation_worker, next_calc.calculation_id, next_calc.parameters)
                            self.active_futures[next_calc.calculation_id] = future

                            # Add callback to clean up completed futures
                            future.add_done_callback(lambda f, calc_id=next_calc.calculation_id: self._cleanup_future(calc_id, f))

                            # Update calculation status from 'waiting' to 'running'
                            self._update_calculation_status_from_queue(next_calc.calculation_id, 'running')

                            started_count += 1
                            logger.info(f"Started queued calculation {next_calc.calculation_id} with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(self.active_futures)}/{self.max_parallel_instances} slots used, {len(self.calculation_queue)} remaining in queue)")
                            started = True
                            break

                        except RuntimeError as e:
                            if 'cannot schedule new futures after shutdown' in str(e):
                                logger.error(f"Cannot schedule calculation {next_calc.calculation_id}: executor has been shut down")
                                self._update_calculation_status_from_queue(next_calc.calculation_id, 'error', 'Process manager has been shut down')
                                # Unregister resources on failure (if resource manager is available)
                                if self.resource_manager is not None:
                                    try:
                                        self.resource_manager.unregister_calculation(next_calc.calculation_id)
                                    except Exception as cleanup_error:
                                        logger.warning(f"Failed to unregister queued calculation resources on failure: {cleanup_error}")
                                started = True
                                break
                            else:
                                raise
                        except Exception as e:
                            logger.error(f"Failed to start queued calculation {next_calc.calculation_id}: {e}")
                            # Unregister resources on failure (if resource manager is available)
                            if self.resource_manager is not None:
                                try:
                                    self.resource_manager.unregister_calculation(next_calc.calculation_id)
                                except Exception as cleanup_error:
                                    logger.warning(f"Failed to unregister queued calculation resources on failure: {cleanup_error}")
                            # Update status to error
                            self._update_calculation_status_from_queue(next_calc.calculation_id, 'error', str(e))
                            started = True  # Continue processing queue
                            break
                    
                    # If SHOULD_QUEUE, keep this calculation in queue and check next one
                    else:
                        logger.debug(f"Calculation {queued_calc.calculation_id} should remain in queue: {reason}")
                
                # If no calculation could be started in this iteration, stop processing
                if not started:
                    if self.calculation_queue:
                        logger.info(f"Queue processing paused - {len(self.calculation_queue)} calculations waiting for resources")
                        # Log first few waiting calculations for debugging
                        for i, queued_calc in enumerate(self.calculation_queue[:3]):
                            logger.debug(f"  Waiting calculation {i+1}: {queued_calc.calculation_id} (reason: {queued_calc.waiting_reason})")
                    break
                
        finally:
            # Always release the queue processing flag
            with self._queue_lock:
                self._queue_processing = False
            
            if started_count > 0:
                logger.info(f"Queue processing completed - started {started_count} calculations from queue")
    
    def _update_calculation_status_from_queue(self, calculation_id: str, status: str, error_message: Optional[str] = None):
        """Update calculation status when transitioning from queue."""
        try:
            from quantum_calc.file_manager import CalculationFileManager
            from quantum_calc import get_current_settings
            settings = get_current_settings()
            file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
            calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
            file_manager.save_calculation_status(calc_dir, status)
            
            # Save error information if this is an error status
            if status == 'error' and error_message:
                file_manager.save_calculation_results(calc_dir, {'error': error_message})
                logger.info(f"Updated status for calculation {calculation_id} to 'error': {error_message}")
            else:
                logger.info(f"Updated status for calculation {calculation_id} from 'waiting' to '{status}'")
            
            # Send immediate WebSocket notification for status changes from queue
            self._send_websocket_notification(calculation_id, status, error_message)
                
        except Exception as e:
            logger.error(f"Failed to update status for calculation {calculation_id}: {e}")
    
    def is_running(self, calculation_id: str) -> bool:
        """Check if a calculation is currently running."""
        future = self.active_futures.get(calculation_id)
        return future is not None and not future.done()
    
    def register_completion_callback(self, calculation_id: str, callback: Callable[[str, bool, Optional[str]], None]):
        """
        Register a callback to be called when a calculation completes.
        
        Args:
            calculation_id: ID of the calculation to monitor
            callback: Function to call when calculation completes. 
                     Signature: callback(calculation_id, success, error_message)
        """
        if calculation_id not in self.completion_callbacks:
            self.completion_callbacks[calculation_id] = []
        self.completion_callbacks[calculation_id].append(callback)
    
    def get_active_calculations(self) -> list:
        """Get list of active calculation IDs."""
        return [calc_id for calc_id, future in self.active_futures.items() if not future.done()]
    
    def get_queued_calculations(self) -> list:
        """Get list of queued calculation IDs."""
        return [calc.calculation_id for calc in self.calculation_queue]
    
    def get_queue_status(self) -> dict:
        """Get current queue status information."""
        return {
            'active_calculations': len(self.active_futures),
            'queued_calculations': len(self.calculation_queue),
            'max_parallel_instances': self.max_parallel_instances,
            'max_workers': self.max_workers
        }
    
    def pause_calculation(self, calculation_id: str) -> bool:
        """
        Request pause for a running calculation.

        Args:
            calculation_id: ID of the calculation to pause

        Returns:
            True if pause request was accepted, False otherwise

        Raises:
            ValueError: If calculation is not found or not in a pausable state
        """
        logger.info(f"Pause requested for calculation: {calculation_id}")

        # Get file manager instance
        from quantum_calc.file_manager import CalculationFileManager
        from quantum_calc import get_current_settings
        settings = get_current_settings()
        file_manager = CalculationFileManager(base_dir=settings.calculations_directory)

        # Check if calculation directory exists
        calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
        if not os.path.exists(calc_dir):
            raise ValueError(f"Calculation not found: {calculation_id}")

        # Check if calculation is running
        status, _ = file_manager.read_calculation_status_details(calc_dir)
        if status != 'running':
            raise ValueError(f"Calculation is not running (status: {status})")

        # Create pause flag file for worker process
        pause_manager.create_pause_flag_file(calc_dir)
        pause_manager.request_pause(calculation_id)

        # Update status to 'pausing' (intermediate state)
        file_manager.save_calculation_status(calc_dir, 'pausing')

        # Send WebSocket notification
        self._send_websocket_notification(calculation_id, 'pausing', None)

        logger.info(f"Pause request accepted for calculation: {calculation_id}")
        return True

    def resume_calculation(self, calculation_id: str) -> Dict[str, Any]:
        """
        Resume a paused calculation from checkpoint.

        Args:
            calculation_id: ID of the calculation to resume

        Returns:
            Calculation instance dictionary

        Raises:
            ValueError: If calculation cannot be resumed
        """
        logger.info(f"Resuming calculation: {calculation_id}")

        # Get file manager instance
        from quantum_calc.file_manager import CalculationFileManager
        from quantum_calc import get_current_settings
        settings = get_current_settings()
        file_manager = CalculationFileManager(base_dir=settings.calculations_directory)

        # Check if calculation directory exists
        calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
        if not os.path.exists(calc_dir):
            raise ValueError(f"Calculation not found: {calculation_id}")

        # Check if calculation is paused
        status, _ = file_manager.read_calculation_status_details(calc_dir)
        if status != 'paused':
            raise ValueError(f"Calculation is not paused (status: {status})")

        # Load pause state (optional - for future checkpoint resume)
        pause_state = file_manager.load_pause_state(calc_dir)

        # Load original parameters
        params = file_manager.read_calculation_parameters(calc_dir)
        if not params:
            raise ValueError("No calculation parameters found")

        # Add resume flag for future use
        params['resume_from_pause'] = True
        if pause_state:
            params['pause_state'] = pause_state

        # Resubmit calculation
        success, message, waiting_reason = self.submit_calculation(calculation_id, params)

        if not success:
            raise ValueError(f"Failed to resume calculation: {message}")

        logger.info(f"Calculation resumed: {calculation_id}")

        # Return simple confirmation - actual calculation instance
        # will be retrieved by service layer after this method returns
        return {'calculation_id': calculation_id, 'status': message}

    
    def shutdown(self, wait: bool = True, timeout: Optional[float] = None):
        """
        Shutdown the process pool.
        
        Args:
            wait: Whether to wait for running calculations to complete
            timeout: Maximum time to wait for shutdown
        """
        if self._shutdown:
            return
            
        self._shutdown = True
        
        # Stop resource monitoring first
        self._stop_resource_monitoring()
        
        if self.executor:
            logger.info(f"Shutting down process pool with {len(self.active_futures)} active calculations")
            
            try:
                self.executor.shutdown(wait=wait, cancel_futures=not wait)
                logger.info("Process pool shut down successfully")
            except Exception as e:
                logger.error(f"Error during process pool shutdown: {e}")
            finally:
                self.executor = None
                self.active_futures.clear()
    
    def _send_websocket_notification(self, calculation_id: str, status: str, error_message: Optional[str] = None):
        """Send immediate WebSocket notification for calculation status changes."""
        if self.notification_callback is None:
            logger.debug(f"WebSocket notification not available for calculation {calculation_id}")
            return
        
        try:
            self.notification_callback(calculation_id, status, error_message)
            logger.debug(f"Sent WebSocket notification for calculation {calculation_id} with status {status}")
        except Exception as e:
            logger.warning(f"Failed to send WebSocket notification for calculation {calculation_id}: {e}")
    
    def __del__(self):
        """Ensure proper cleanup when the manager is destroyed."""
        if not self._shutdown:
            self.shutdown(wait=False)


# Global process manager instance
_process_manager: Optional[CalculationProcessManager] = None


def initialize_process_manager_with_callback(notification_callback: Optional[Callable] = None,
                                             max_parallel_instances: Optional[int] = None) -> CalculationProcessManager:
    """
    Initialize the global process manager with a notification callback.
    This should be called once during application startup.
    
    Args:
        notification_callback: Function to call for WebSocket notifications.
                               Signature: callback(calculation_id, status, error_message)
        max_parallel_instances: Maximum number of parallel calculation instances.
    
    Returns:
        The initialized CalculationProcessManager instance.
    """
    global _process_manager
    
    if _process_manager is not None:
        logger.warning("Process manager already initialized. Updating notification callback.")
        _process_manager.notification_callback = notification_callback
        if max_parallel_instances is not None:
            _process_manager.set_max_parallel_instances(max_parallel_instances)
        return _process_manager
    
    # Determine max parallel instances
    if max_parallel_instances is None:
        max_parallel_instances = min(4, multiprocessing.cpu_count())
    
    logger.info(f"Initializing process manager with notification callback, max_parallel={max_parallel_instances}")
    
    try:
        _process_manager = CalculationProcessManager(
            max_parallel_instances=max_parallel_instances,
            notification_callback=notification_callback
        )
        logger.info("Process manager initialized successfully with notification callback")
    except Exception as e:
        logger.error(f"Failed to initialize process manager: {e}")
        _process_manager = None
        raise ProcessManagerError(f"Failed to initialize process manager: {str(e)}")
    
    return _process_manager


def get_process_manager() -> CalculationProcessManager:
    """
    Get the global process manager instance.
    If not already initialized, creates a basic instance without notification callback.
    For proper initialization with callback, use initialize_process_manager_with_callback().
    """
    global _process_manager
    if _process_manager is None:
        # Initialize with default values without callback
        # This is a fallback - proper initialization should use initialize_process_manager_with_callback()
        default_max_parallel = min(4, multiprocessing.cpu_count())
        logger.warning(f"Process manager not initialized with callback. Using default settings: max_parallel={default_max_parallel}")
        
        try:
            _process_manager = CalculationProcessManager(max_parallel_instances=default_max_parallel)
            logger.info("Process manager initialized with defaults (no callback)")
        except Exception as e:
            logger.error(f"Failed to initialize process manager: {e}")
            _process_manager = None
            raise ProcessManagerError(f"Failed to initialize process manager: {str(e)}")
    
    return _process_manager


def update_process_manager_settings():
    """Update process manager settings from settings manager (called after initialization)."""
    global _process_manager
    if _process_manager is not None:
        try:
            from .settings_manager import get_current_settings
            settings = get_current_settings()
            _process_manager.set_max_parallel_instances(settings.max_parallel_instances)
            logger.info(f"Updated process manager settings: max_parallel_instances={settings.max_parallel_instances}")
        except Exception as e:
            logger.warning(f"Failed to update process manager settings: {e}. Keeping current settings.")


def shutdown_process_manager():
    """Shutdown the global process manager."""
    global _process_manager
    if _process_manager is not None:
        _process_manager.shutdown()
        _process_manager = None