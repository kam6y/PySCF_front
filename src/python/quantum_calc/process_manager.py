"""Process pool manager for CPU-bound quantum chemistry calculations."""

import os
import sys
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, Future
from typing import Dict, Any, Optional, Callable, List
import time
from datetime import datetime
from threadpoolctl import threadpool_limits
from queue import Queue
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class QueuedCalculation:
    """Represents a calculation waiting in the queue."""
    calculation_id: str
    parameters: dict
    created_at: datetime
    
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
            'spin': parameters['spin_multiplicity'],
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
        
        # Initialize resource manager
        from .resource_manager import get_resource_manager
        self.resource_manager = get_resource_manager()
        
        # Create executor immediately with fixed worker count
        self.executor = ProcessPoolExecutor(
            max_workers=self.max_workers,
            mp_context=multiprocessing.get_context('spawn')  # Better cross-platform compatibility
        )
        
        logger.info(f"Initialized CalculationProcessManager with {self.max_workers} worker processes and {self.max_parallel_instances} max parallel instances")
    
    
    
    
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
            self._process_queue()
    
    def submit_calculation(self, calculation_id: str, parameters: dict) -> tuple[bool, str]:
        """
        Submit a calculation to the process pool or queue with resource checking.
        
        Args:
            calculation_id: Unique identifier for the calculation
            parameters: Calculation parameters
            
        Returns:
            Tuple of (success: bool, status: str) where status is 'running', 'waiting', or 'error'
        """
        if self._shutdown:
            logger.error("Cannot submit calculation: process manager is shut down")
            return False, 'error'
        
        if self.executor is None:
            logger.error("Process pool executor is not available")
            return False, 'error'
        
        user_cpu_cores = parameters.get('cpu_cores') or 1
        user_memory_mb = parameters.get('memory_mb') or 2000
        calculation_method = parameters.get('calculation_method', 'DFT')
        
        # Check resource availability
        can_allocate, reason = self.resource_manager.can_allocate_resources(
            cpu_cores=user_cpu_cores,
            memory_mb=user_memory_mb,
            calculation_method=calculation_method
        )
        
        if not can_allocate:
            logger.warning(f"Cannot start calculation {calculation_id}: {reason}")
            # Add to queue if resource constraints prevent immediate execution
            queued_calc = QueuedCalculation(
                calculation_id=calculation_id,
                parameters=parameters,
                created_at=datetime.now()
            )
            self.calculation_queue.append(queued_calc)
            self.calculation_queue.sort(key=lambda x: x.created_at)
            
            logger.info(f"Added calculation {calculation_id} to queue due to resource constraints: {reason}")
            return True, 'waiting'
        
        # Check if we can start the calculation immediately (slot availability)
        if len(self.active_futures) < self.max_parallel_instances:
            # Start calculation immediately
            try:
                # Register resources with the resource manager
                self.resource_manager.register_calculation(
                    calculation_id=calculation_id,
                    cpu_cores=user_cpu_cores,
                    memory_mb=user_memory_mb,
                    calculation_method=calculation_method
                )
                
                future = self.executor.submit(calculation_worker, calculation_id, parameters)
                self.active_futures[calculation_id] = future
                
                # Add callback to clean up completed futures and process queue
                future.add_done_callback(lambda f: self._cleanup_future(calculation_id, f))
                
                logger.info(f"Started calculation {calculation_id} immediately with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(self.active_futures)}/{self.max_parallel_instances} slots used)")
                return True, 'running'
                
            except Exception as e:
                logger.error(f"Failed to submit calculation {calculation_id}: {e}")
                # Unregister resources on failure
                self.resource_manager.unregister_calculation(calculation_id)
                return False, 'error'
        else:
            # Add to queue
            queued_calc = QueuedCalculation(
                calculation_id=calculation_id,
                parameters=parameters,
                created_at=datetime.now()
            )
            self.calculation_queue.append(queued_calc)
            
            # Sort queue by creation time (FIFO)
            self.calculation_queue.sort(key=lambda x: x.created_at)
            
            logger.info(f"Added calculation {calculation_id} to queue (position {len(self.calculation_queue)}) - all {self.max_parallel_instances} slots occupied")
            return True, 'waiting'
    
    def _cleanup_future(self, calculation_id: str, future: Future):
        """Clean up a completed future, call completion callbacks, and process queue."""
        success = False
        error_message = None
        
        try:
            if calculation_id in self.active_futures:
                del self.active_futures[calculation_id]
            
            # Unregister resources from resource manager
            self.resource_manager.unregister_calculation(calculation_id)
            
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
            self._process_queue()
                        
        except Exception as e:
            logger.error(f"Error cleaning up future for {calculation_id}: {e}")
    
    def _process_queue(self):
        """Process the calculation queue and start next calculations if slots and resources are available."""
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
                
                # Check if resources are available for this calculation
                can_allocate, reason = self.resource_manager.can_allocate_resources(
                    cpu_cores=user_cpu_cores,
                    memory_mb=user_memory_mb,
                    calculation_method=calculation_method
                )
                
                if can_allocate:
                    # Remove this calculation from queue and start it
                    next_calc = self.calculation_queue.pop(i)
                    calc_found = True
                    
                    try:
                        # Register resources with the resource manager
                        self.resource_manager.register_calculation(
                            calculation_id=next_calc.calculation_id,
                            cpu_cores=user_cpu_cores,
                            memory_mb=user_memory_mb,
                            calculation_method=calculation_method
                        )
                        
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
                        # Unregister resources on failure
                        self.resource_manager.unregister_calculation(next_calc.calculation_id)
                        # Update status to error
                        self._update_calculation_status_from_queue(next_calc.calculation_id, 'error')
                        break
                else:
                    logger.debug(f"Cannot start queued calculation {queued_calc.calculation_id}: {reason}")
            
            # If no calculation in the queue can be started due to resource constraints, break
            if not calc_found:
                logger.debug("No queued calculations can be started due to resource constraints")
                break
    
    def _update_calculation_status_from_queue(self, calculation_id: str, status: str):
        """Update calculation status when transitioning from queue to running."""
        try:
            from quantum_calc.file_manager import CalculationFileManager
            file_manager = CalculationFileManager()
            calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
            file_manager.save_calculation_status(calc_dir, status)
            logger.info(f"Updated status for calculation {calculation_id} from 'waiting' to '{status}'")
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
        # Get current settings to determine max parallel instances
        try:
            from .settings_manager import get_current_settings
            settings = get_current_settings()
            max_parallel = settings.max_parallel_instances
        except Exception as e:
            logger.warning(f"Failed to load settings for process manager: {e}. Using default.")
            max_parallel = min(4, multiprocessing.cpu_count())
        
        _process_manager = CalculationProcessManager(max_parallel_instances=max_parallel)
    return _process_manager


def shutdown_process_manager():
    """Shutdown the global process manager."""
    global _process_manager
    if _process_manager is not None:
        _process_manager.shutdown()
        _process_manager = None