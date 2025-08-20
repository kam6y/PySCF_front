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

logger = logging.getLogger(__name__)


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
            'spin': (parameters['spin_multiplicity'] - 1) // 2,
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
    """Manages a process pool for quantum chemistry calculations."""
    
    def __init__(self, max_workers: Optional[int] = None):
        """
        Initialize the process manager with a fixed process pool.
        
        Args:
            max_workers: Maximum number of worker processes. If None, defaults to CPU count.
        """
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.active_futures: Dict[str, Future] = {}
        self.completion_callbacks: Dict[str, List[Callable]] = {}  # calculation_id -> list of callbacks
        self._shutdown = False
        
        # Create executor immediately with fixed worker count
        self.executor = ProcessPoolExecutor(
            max_workers=self.max_workers,
            mp_context=multiprocessing.get_context('spawn')  # Better cross-platform compatibility
        )
        
        logger.info(f"Initialized CalculationProcessManager with {self.max_workers} worker processes")
    
    
    
    
    def submit_calculation(self, calculation_id: str, parameters: dict) -> bool:
        """
        Submit a calculation to the process pool.
        
        Args:
            calculation_id: Unique identifier for the calculation
            parameters: Calculation parameters
            
        Returns:
            True if successfully submitted, False otherwise
        """
        if self._shutdown:
            logger.error("Cannot submit calculation: process manager is shut down")
            return False
        
        if self.executor is None:
            logger.error("Process pool executor is not available")
            return False
        
        user_cpu_cores = parameters.get('cpu_cores') or 1
        logger.info(f"Submitting calculation {calculation_id} with {user_cpu_cores} CPU cores to process pool ({self.max_workers} workers available)")
        
        try:
            future = self.executor.submit(calculation_worker, calculation_id, parameters)
            self.active_futures[calculation_id] = future
            
            # Add callback to clean up completed futures
            future.add_done_callback(lambda f: self._cleanup_future(calculation_id, f))
            
            logger.info(f"Successfully submitted calculation {calculation_id} to process pool")
            return True
            
        except Exception as e:
            logger.error(f"Failed to submit calculation {calculation_id}: {e}")
            return False
    
    def _cleanup_future(self, calculation_id: str, future: Future):
        """Clean up a completed future and call completion callbacks."""
        success = False
        error_message = None
        
        try:
            if calculation_id in self.active_futures:
                del self.active_futures[calculation_id]
            
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
                        
        except Exception as e:
            logger.error(f"Error cleaning up future for {calculation_id}: {e}")
    
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
    
    def cancel_calculation(self, calculation_id: str) -> bool:
        """
        Attempt to cancel a calculation.
        
        Args:
            calculation_id: ID of the calculation to cancel
            
        Returns:
            True if successfully cancelled, False otherwise
        """
        future = self.active_futures.get(calculation_id)
        if future and not future.done():
            cancelled = future.cancel()
            if cancelled:
                logger.info(f"Successfully cancelled calculation {calculation_id}")
                # Update status to cancelled on the file system
                try:
                    from quantum_calc.file_manager import CalculationFileManager
                    file_manager = CalculationFileManager()
                    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
                    file_manager.save_calculation_status(calc_dir, 'cancelled')
                except Exception as e:
                    logger.error(f"Failed to update status for cancelled calculation {calculation_id}: {e}")
            return cancelled
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
        _process_manager = CalculationProcessManager()
    return _process_manager


def shutdown_process_manager():
    """Shutdown the global process manager."""
    global _process_manager
    if _process_manager is not None:
        _process_manager.shutdown()
        _process_manager = None