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
        return self.created_at < other.created_at


def calculation_worker(calculation_id: str, parameters: dict) -> tuple:
    """
    Worker function to run quantum chemistry calculations in a separate process.
    Returns (success: bool, error_message: str or None)
    """
    # Get user-specified CPU cores or default to 1
    cpu_cores = parameters.get('cpu_cores') or 1
    cpu_cores_str = str(int(cpu_cores))
    
    # Set all parallel processing environment variables to control CPU usage
    # This ensures that user-specified CPU count is respected by all underlying libraries
    os.environ['OMP_NUM_THREADS'] = cpu_cores_str
    os.environ['MKL_NUM_THREADS'] = cpu_cores_str       # Intel MKL
    os.environ['OPENBLAS_NUM_THREADS'] = cpu_cores_str  # OpenBLAS
    os.environ['BLIS_NUM_THREADS'] = cpu_cores_str      # BLIS
    os.environ['VECLIB_MAXIMUM_THREADS'] = cpu_cores_str # macOS Accelerate Framework
    os.environ['NUMEXPR_NUM_THREADS'] = cpu_cores_str   # NumExpr
    
    memory_mb = parameters.get('memory_mb') or 2000  # Noneや空の値をデフォルト値に置き換え

    # Import here to avoid issues with multiprocessing and module loading
    from quantum_calc import DFTCalculator, HFCalculator, MP2Calculator, CCSDCalculator, TDDFTCalculator
    from quantum_calc import CalculationError, ConvergenceError, InputError
    from quantum_calc.file_manager import CalculationFileManager
    from threadpoolctl import threadpool_info
    from pyscf import lib
    
    # Setup logging for this process
    process_logger = logging.getLogger(f'worker_{calculation_id}')
    process_logger.setLevel(logging.INFO)
    
    file_manager = CalculationFileManager()
    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)

    try:
        # Update status to running on the file system
        file_manager.save_calculation_status(calc_dir, 'running')
        process_logger.info(f"Starting calculation {calculation_id} in process {os.getpid()}")
        process_logger.info(f"Using {cpu_cores} CPU cores per process for calculation.")
        process_logger.info(f"Set environment variables: OMP_NUM_THREADS={cpu_cores_str}, MKL_NUM_THREADS={cpu_cores_str}, OPENBLAS_NUM_THREADS={cpu_cores_str}")
        
        # Log detected threadpool libraries for debugging
        try:
            thread_info = threadpool_info()
            process_logger.info(f"Detected threadpool libraries: {len(thread_info)} found")
            for info in thread_info:
                process_logger.info(f"  {info.get('user_api', 'unknown')}: {info.get('internal_api', 'unknown')} - threads: {info.get('num_threads', 'unknown')}")
        except Exception as e:
            process_logger.warning(f"Could not get threadpool info: {e}")
        
        # Set PySCF thread count using official API
        original_threads = None
        try:
            original_threads = lib.num_threads()
            lib.num_threads(int(cpu_cores))
            new_threads = lib.num_threads()
            process_logger.info(f"PySCF threads: {original_threads} -> {new_threads} (requested: {cpu_cores})")
        except Exception as e:
            process_logger.warning(f"Could not set PySCF threads: {e}")
        
        # Initialize calculator based on calculation method
        calculation_method = parameters.get('calculation_method', 'DFT')
        if calculation_method == 'HF':
            calculator = HFCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        elif calculation_method == 'MP2':
            calculator = MP2Calculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        elif calculation_method == 'CCSD':
            calculator = CCSDCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        elif calculation_method == 'CCSD_T':
            calculator = CCSDCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        elif calculation_method == 'TDDFT':
            calculator = TDDFTCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        else:  # Default to DFT
            calculator = DFTCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        
        # Parse XYZ and setup calculation
        atoms = calculator.parse_xyz(parameters['xyz'])
        
        # Prepare setup parameters
        setup_params = {
            'basis': parameters['basis_function'],
            'charge': parameters['charges'],
            'spin': parameters['spin'],
            'max_cycle': 150,
            'solvent_method': parameters['solvent_method'],
            'solvent': parameters['solvent'],
            'memory_mb': memory_mb,  # メモリ設定を渡す (必ず有効な値)
        }
        
        # Add exchange-correlation functional for DFT and TDDFT calculations
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
        
        calculator.setup_calculation(atoms, **setup_params)
        
        # Run calculation with controlled BLAS/LAPACK threading
        process_logger.info(f"Executing calculation with threadpool_limits(limits={cpu_cores}, user_api='blas')")
        with threadpool_limits(limits=int(cpu_cores), user_api='blas'):
            results = calculator.run_calculation()
        
        # Save results and update status to completed
        file_manager.save_calculation_results(calc_dir, results)
        file_manager.save_calculation_status(calc_dir, 'completed')
        
        process_logger.info(f"Calculation {calculation_id} completed successfully in process {os.getpid()}")
        return True, None

    except (InputError, ConvergenceError, CalculationError) as e:
        process_logger.error(f"Calculation {calculation_id} failed: {e}")
        file_manager.save_calculation_status(calc_dir, 'error')
        # Save error information to results.json
        file_manager.save_calculation_results(calc_dir, {'error': str(e)})
        return False, str(e)

    except Exception as e:
        process_logger.error(f"Unexpected error in calculation {calculation_id}: {e}", exc_info=True)
        file_manager.save_calculation_status(calc_dir, 'error')
        file_manager.save_calculation_results(calc_dir, {'error': 'An unexpected internal server error occurred.'})
        return False, 'An unexpected internal server error occurred.'
    
    finally:
        # Restore original PySCF thread count
        if original_threads is not None:
            try:
                lib.num_threads(original_threads)
                process_logger.info(f"Restored PySCF threads to {original_threads}")
            except Exception as e:
                process_logger.warning(f"Could not restore PySCF threads: {e}")


class CalculationProcessManager:
    """Manages a process pool for quantum chemistry calculations with queueing support and resource management."""
    
    def __init__(self, max_workers: Optional[int] = None, max_parallel_instances: Optional[int] = None):
        """
        Initialize the process manager with a fixed process pool and queueing system.
        
        Args:
            max_workers: Maximum number of worker processes. If None, defaults to CPU count.
            max_parallel_instances: Maximum number of parallel calculation instances. If None, defaults to max_workers.
        """
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
                # Only check for improvements if there are queued calculations and resource manager is available
                if self.calculation_queue and not self._shutdown and self.resource_manager is not None:
                    try:
                        # Check if resources have improved
                        has_improved, reason = self.resource_manager.has_resources_improved()
                        
                        if has_improved:
                            logger.info(f"Resources improved: {reason}. Processing queue...")
                            self._process_queue()
                            
                            # Also check for calculations that should be marked as errors
                            self._check_queue_for_resource_errors()
                    except Exception as e:
                        logger.warning(f"Failed to check resource improvements: {e}")
                        # Try to process queue anyway
                        if self.calculation_queue:
                            self._process_queue()
                    
                # Wait for the next check or until stop event is set
                self._resource_monitor_stop_event.wait(timeout=self._resource_monitor_interval)
                
            except Exception as e:
                logger.error(f"Error in resource monitoring loop: {e}")
                # Continue monitoring despite errors
                self._resource_monitor_stop_event.wait(timeout=self._resource_monitor_interval)
        
        logger.info("Resource monitoring loop stopped")
    
    def _check_queue_for_resource_errors(self):
        """Check queued calculations for those that should be marked as errors due to insufficient system resources."""
        if not self.calculation_queue:
            return
            
        calculations_to_error = []
        
        for queued_calc in self.calculation_queue:
            # Check if system resources are fundamentally insufficient (if resource manager is available)
            insufficient = False
            reason = None
            
            if self.resource_manager is not None:
                try:
                    insufficient, reason = self.resource_manager.check_if_system_resources_insufficient()
                except Exception as e:
                    logger.warning(f"Failed to check system resources insufficiency for queued calculation: {e}")
                    insufficient = False
            
            if insufficient:
                logger.warning(f"Marking calculation {queued_calc.calculation_id} as error: {reason}")
                calculations_to_error.append((queued_calc, reason))
        
        # Remove from queue and update status to error
        for calc, error_reason in calculations_to_error:
            self.calculation_queue.remove(calc)
            self._update_calculation_status_from_queue(calc.calculation_id, 'error', error_reason)
    
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
        
        user_cpu_cores = parameters.get('cpu_cores') or 1
        user_memory_mb = parameters.get('memory_mb') or 2000
        calculation_method = parameters.get('calculation_method', 'DFT')
        
        # Check resource availability (skip if resource manager is not available)
        can_allocate = True
        reason = "Resource checking disabled"
        
        if self.resource_manager is not None:
            try:
                can_allocate, reason = self.resource_manager.can_allocate_resources(
                    cpu_cores=user_cpu_cores,
                    memory_mb=user_memory_mb,
                    calculation_method=calculation_method
                )
            except Exception as e:
                logger.warning(f"Resource checking failed: {e}. Proceeding without resource constraints.")
                can_allocate = True
                reason = "Resource checking failed - proceeding without constraints"
        
        if not can_allocate:
            logger.warning(f"Cannot start calculation {calculation_id}: {reason}")
            
            # Check if this is a fundamental system resource insufficiency (no active calculations)
            # Only check if resource manager is available
            insufficient = False
            insufficient_reason = None
            
            if self.resource_manager is not None:
                try:
                    insufficient, insufficient_reason = self.resource_manager.check_if_system_resources_insufficient()
                except Exception as e:
                    logger.warning(f"Failed to check system resource insufficiency: {e}")
                    insufficient = False
            
            if insufficient:
                # System resources are fundamentally insufficient - return error instead of waiting
                logger.error(f"System resources insufficient for calculation {calculation_id}: {insufficient_reason}")
                return False, 'error', insufficient_reason
            else:
                # Add to queue if resource constraints prevent immediate execution but system is generally sufficient
                queued_calc = QueuedCalculation(
                    calculation_id=calculation_id,
                    parameters=parameters,
                    created_at=datetime.now(),
                    waiting_reason=reason
                )
                self.calculation_queue.append(queued_calc)
                self.calculation_queue.sort(key=lambda x: x.created_at)
                
                logger.info(f"Added calculation {calculation_id} to queue due to resource constraints: {reason}")
                return True, 'waiting', reason
        
        # Check if we can start the calculation immediately (slot availability)
        if len(self.active_futures) < self.max_parallel_instances:
            # Start calculation immediately
            try:
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
                
                future = self.executor.submit(calculation_worker, calculation_id, parameters)
                self.active_futures[calculation_id] = future
                
                # Add callback to clean up completed futures and process queue
                future.add_done_callback(lambda f: self._cleanup_future(calculation_id, f))
                
                logger.info(f"Started calculation {calculation_id} immediately with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(self.active_futures)}/{self.max_parallel_instances} slots used)")
                return True, 'running', None
                
            except Exception as e:
                logger.error(f"Failed to submit calculation {calculation_id}: {e}")
                # Unregister resources on failure (if resource manager is available)
                if self.resource_manager is not None:
                    try:
                        self.resource_manager.unregister_calculation(calculation_id)
                    except Exception as cleanup_error:
                        logger.warning(f"Failed to unregister calculation resources on failure: {cleanup_error}")
                return False, 'error', None
        else:
            # Add to queue when all slots are occupied
            slot_reason = f"All calculation slots occupied ({len(self.active_futures)}/{self.max_parallel_instances})"
            queued_calc = QueuedCalculation(
                calculation_id=calculation_id,
                parameters=parameters,
                created_at=datetime.now(),
                waiting_reason=slot_reason
            )
            self.calculation_queue.append(queued_calc)
            
            # Sort queue by creation time (FIFO)
            self.calculation_queue.sort(key=lambda x: x.created_at)
            
            logger.info(f"Added calculation {calculation_id} to queue (position {len(self.calculation_queue)}) - all {self.max_parallel_instances} slots occupied")
            return True, 'waiting', slot_reason
    
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
                error_message = str(future.exception())
                logger.error(f"Calculation {calculation_id} failed with exception: {error_message}")
            else:
                success, calc_error = future.result()
                if success:
                    logger.info(f"Calculation {calculation_id} completed successfully")
                else:
                    error_message = calc_error
                    logger.warning(f"Calculation {calculation_id} failed: {error_message}")
            
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
            processed_calculations = []
            
            while (len(self.active_futures) < self.max_parallel_instances and 
                   self.calculation_queue and 
                   not self._shutdown):
                
                # Find the first calculation in queue that can be started with current resources
                calc_found = False
                for i, queued_calc in enumerate(self.calculation_queue):
                    user_cpu_cores = queued_calc.parameters.get('cpu_cores') or 1
                    user_memory_mb = queued_calc.parameters.get('memory_mb') or 2000
                    calculation_method = queued_calc.parameters.get('calculation_method', 'DFT')
                    
                    # Check if resources are available for this calculation (if resource manager is available)
                    can_allocate = True
                    reason = "Resource checking disabled"
                    
                    if self.resource_manager is not None:
                        try:
                            can_allocate, reason = self.resource_manager.can_allocate_resources(
                                cpu_cores=user_cpu_cores,
                                memory_mb=user_memory_mb,
                                calculation_method=calculation_method
                            )
                        except Exception as e:
                            logger.warning(f"Resource checking failed in queue processing: {e}")
                            can_allocate = True
                            reason = "Resource checking failed"
                    
                    if can_allocate:
                        # Remove this calculation from queue and start it
                        next_calc = self.calculation_queue.pop(i)
                        calc_found = True
                        
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
                            future = self.executor.submit(calculation_worker, next_calc.calculation_id, next_calc.parameters)
                            self.active_futures[next_calc.calculation_id] = future
                            
                            # Add callback to clean up completed futures
                            future.add_done_callback(lambda f, calc_id=next_calc.calculation_id: self._cleanup_future(calc_id, f))
                            
                            # Update calculation status from 'waiting' to 'running'
                            self._update_calculation_status_from_queue(next_calc.calculation_id, 'running')
                            
                            logger.info(f"Started queued calculation {next_calc.calculation_id} with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(self.active_futures)}/{self.max_parallel_instances} slots used, {len(self.calculation_queue)} remaining in queue)")
                            break
                            
                        except Exception as e:
                            logger.error(f"Failed to start queued calculation {next_calc.calculation_id}: {e}")
                            # Unregister resources on failure (if resource manager is available)
                            if self.resource_manager is not None:
                                try:
                                    self.resource_manager.unregister_calculation(next_calc.calculation_id)
                                except Exception as cleanup_error:
                                    logger.warning(f"Failed to unregister queued calculation resources on failure: {cleanup_error}")
                            # Update status to error
                            self._update_calculation_status_from_queue(next_calc.calculation_id, 'error')
                            break
                    else:
                        logger.debug(f"Cannot start queued calculation {queued_calc.calculation_id}: {reason}")
                
                # If no calculation in the queue can be started due to resource constraints, break
                if not calc_found:
                    if self.calculation_queue:
                        logger.warning(f"No queued calculations can be started due to resource constraints. {len(self.calculation_queue)} calculations still waiting.")
                        # Log details of waiting calculations
                        for i, queued_calc in enumerate(self.calculation_queue[:3]):  # Log first 3 for debugging
                            logger.info(f"Waiting calculation {i+1}: {queued_calc.calculation_id} (reason: {queued_calc.waiting_reason})")
                    else:
                        logger.debug("Queue is empty, no calculations to process")
                    break
                
        finally:
            # Always release the queue processing flag
            with self._queue_lock:
                self._queue_processing = False
            
            if processed_calculations:
                logger.info(f"Queue processing completed - started {len(processed_calculations)} calculations from queue")
    
    def _update_calculation_status_from_queue(self, calculation_id: str, status: str, error_message: Optional[str] = None):
        """Update calculation status when transitioning from queue."""
        try:
            from quantum_calc.file_manager import CalculationFileManager
            file_manager = CalculationFileManager()
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
    
    def cancel_calculation(self, calculation_id: str) -> bool:
        """
        Attempt to cancel a calculation (either running or queued).
        
        Args:
            calculation_id: ID of the calculation to cancel
            
        Returns:
            True if successfully cancelled, False otherwise
        """
        # Check if calculation is running
        future = self.active_futures.get(calculation_id)
        if future and not future.done():
            cancelled = future.cancel()
            if cancelled:
                logger.info(f"Successfully cancelled running calculation {calculation_id}")
                # Update status to cancelled on the file system
                try:
                    from quantum_calc.file_manager import CalculationFileManager
                    file_manager = CalculationFileManager()
                    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
                    file_manager.save_calculation_status(calc_dir, 'cancelled')
                except Exception as e:
                    logger.error(f"Failed to update status for cancelled calculation {calculation_id}: {e}")
            return cancelled
        
        # Check if calculation is queued
        for i, calc in enumerate(self.calculation_queue):
            if calc.calculation_id == calculation_id:
                # Remove from queue
                removed_calc = self.calculation_queue.pop(i)
                logger.info(f"Successfully cancelled queued calculation {calculation_id} (was at position {i+1} in queue)")
                
                # Update status to cancelled on the file system
                try:
                    from quantum_calc.file_manager import CalculationFileManager
                    file_manager = CalculationFileManager()
                    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
                    file_manager.save_calculation_status(calc_dir, 'cancelled')
                except Exception as e:
                    logger.error(f"Failed to update status for cancelled calculation {calculation_id}: {e}")
                
                return True
        
        return False
    
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
        try:
            # Import here to avoid circular imports
            import sys
            app_module = sys.modules.get('app')
            if app_module and hasattr(app_module, 'send_immediate_websocket_notification'):
                app_module.send_immediate_websocket_notification(calculation_id, status, error_message)
                logger.debug(f"Sent WebSocket notification for calculation {calculation_id} with status {status}")
            else:
                logger.debug(f"WebSocket notification not available for calculation {calculation_id}")
        except Exception as e:
            logger.warning(f"Failed to send WebSocket notification for calculation {calculation_id}: {e}")
    
    def __del__(self):
        """Ensure proper cleanup when the manager is destroyed."""
        if not self._shutdown:
            self.shutdown(wait=False)


# Global process manager instance
_process_manager: Optional[CalculationProcessManager] = None


def get_process_manager() -> CalculationProcessManager:
    """Get the global process manager instance."""
    global _process_manager
    if _process_manager is None:
        # Initialize with default values to avoid circular imports
        # Settings will be applied later via update_process_manager_settings()
        default_max_parallel = min(4, multiprocessing.cpu_count())
        logger.info(f"Initializing process manager with default settings: max_parallel={default_max_parallel}")
        
        try:
            _process_manager = CalculationProcessManager(max_parallel_instances=default_max_parallel)
            logger.info("Process manager initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize process manager: {e}")
            # Create a minimal process manager that can handle errors gracefully
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