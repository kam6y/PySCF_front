"""Process pool manager for CPU-bound quantum chemistry calculations."""

import os
import sys
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, Future
from typing import Dict, Any, Optional
import time
from datetime import datetime

logger = logging.getLogger(__name__)


def calculation_worker(calculation_id: str, parameters: dict) -> tuple:
    """
    Worker function to run quantum chemistry calculations in a separate process.
    Returns (success: bool, error_message: str or None)
    """
    # Import here to avoid issues with multiprocessing and module loading
    from quantum_calc import DFTCalculator, HFCalculator, MP2Calculator, TDDFTCalculator
    from quantum_calc import CalculationError, ConvergenceError, InputError
    from quantum_calc.file_manager import CalculationFileManager
    
    # Setup logging for this process
    process_logger = logging.getLogger(f'worker_{calculation_id}')
    process_logger.setLevel(logging.INFO)
    
    file_manager = CalculationFileManager()
    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)

    try:
        # Update status to running on the file system
        file_manager.save_calculation_status(calc_dir, 'running')
        process_logger.info(f"Starting calculation {calculation_id} in process {os.getpid()}")
        
        # Configure resource settings
        cpu_cores = parameters.get('cpu_cores')
        memory_mb = parameters.get('memory_mb')
        
        # Set thread counts for parallel libraries if cpu_cores is specified
        if cpu_cores is not None and cpu_cores > 0:
            os.environ['OMP_NUM_THREADS'] = str(cpu_cores)
            os.environ['OPENBLAS_NUM_THREADS'] = str(cpu_cores)
            os.environ['MKL_NUM_THREADS'] = str(cpu_cores)
            process_logger.info(f"Set parallel library thread counts to {cpu_cores}")
        
        # Initialize calculator based on calculation method
        calculation_method = parameters.get('calculation_method', 'DFT')
        if calculation_method == 'HF':
            calculator = HFCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        elif calculation_method == 'MP2':
            calculator = MP2Calculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
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
            'cpu_cores': cpu_cores,
            'memory_mb': memory_mb
        }
        
        # Add exchange-correlation functional for DFT and TDDFT calculations
        if calculation_method in ['DFT', 'TDDFT']:
            setup_params['xc'] = parameters['exchange_correlation']
        
        # Add TDDFT-specific parameters
        if calculation_method == 'TDDFT':
            setup_params['nstates'] = parameters.get('tddft_nstates', 10)
            setup_params['tddft_method'] = parameters.get('tddft_method', 'TDDFT')
            setup_params['analyze_nto'] = parameters.get('tddft_analyze_nto', False)
        
        calculator.setup_calculation(atoms, **setup_params)
        
        # Run calculation
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


class CalculationProcessManager:
    """Manages a process pool for quantum chemistry calculations."""
    
    def __init__(self, max_workers: Optional[int] = None):
        """
        Initialize the process manager.
        
        Args:
            max_workers: Maximum number of worker processes. If None, defaults to CPU count.
        """
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.executor: Optional[ProcessPoolExecutor] = None
        self.active_futures: Dict[str, Future] = {}
        self._shutdown = False
        
        logger.info(f"Initializing CalculationProcessManager with {self.max_workers} max workers")
    
    def start(self):
        """Start the process pool."""
        if self.executor is None and not self._shutdown:
            self.executor = ProcessPoolExecutor(
                max_workers=self.max_workers,
                mp_context=multiprocessing.get_context('spawn')  # Better cross-platform compatibility
            )
            logger.info(f"Started ProcessPoolExecutor with {self.max_workers} workers")
    
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
            self.start()
        
        try:
            future = self.executor.submit(calculation_worker, calculation_id, parameters)
            self.active_futures[calculation_id] = future
            
            # Add callback to clean up completed futures
            future.add_done_callback(lambda f: self._cleanup_future(calculation_id, f))
            
            logger.info(f"Submitted calculation {calculation_id} to process pool")
            return True
            
        except Exception as e:
            logger.error(f"Failed to submit calculation {calculation_id}: {e}")
            return False
    
    def _cleanup_future(self, calculation_id: str, future: Future):
        """Clean up a completed future."""
        try:
            if calculation_id in self.active_futures:
                del self.active_futures[calculation_id]
            
            # Log the result
            if future.exception():
                logger.error(f"Calculation {calculation_id} failed with exception: {future.exception()}")
            else:
                success, error = future.result()
                if success:
                    logger.info(f"Calculation {calculation_id} completed successfully")
                else:
                    logger.warning(f"Calculation {calculation_id} failed: {error}")
                    
        except Exception as e:
            logger.error(f"Error cleaning up future for {calculation_id}: {e}")
    
    def is_running(self, calculation_id: str) -> bool:
        """Check if a calculation is currently running."""
        future = self.active_futures.get(calculation_id)
        return future is not None and not future.done()
    
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