"""
Quantum chemistry calculation service.

This service encapsulates all quantum chemistry calculation business logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
import os
import shutil
import time
from datetime import datetime
from typing import Optional, Dict, Any, Tuple, List

from quantum_calc import (
    get_process_manager, get_all_supported_parameters, get_current_settings,
    InputError, GeometryError, ProcessManagerError, CalculationError, FileManagerError
)
from quantum_calc.file_manager import CalculationFileManager
from quantum_calc.orbital_generator import MolecularOrbitalGenerator
from quantum_calc.ir_spectrum import create_ir_spectrum_from_calculation_results
from .exceptions import (
    ServiceError, NotFoundError, ValidationError,
    ResourceUnavailableError, InsufficientResourcesError, PermissionDeniedError
)

logger = logging.getLogger(__name__)


class QuantumService:
    """Service for quantum chemistry calculation operations."""

    def __init__(self):
        """Initialize QuantumService."""
        # Load calculations directory from settings
        settings = get_current_settings()
        calculations_dir = settings.calculations_directory
        self.file_manager = CalculationFileManager(base_dir=calculations_dir)

    def update_calculations_directory(self, new_directory: str) -> None:
        """
        Update the calculations directory.

        Args:
            new_directory: New directory path for calculations
        """
        logger.info(f"Updating QuantumService calculations directory to: {new_directory}")
        self.file_manager.set_base_directory(new_directory)

    def get_supported_parameters(self) -> Dict[str, Any]:
        """
        Get supported quantum chemistry parameters.
        
        Returns:
            Dict containing supported methods, basis sets, functionals, etc.
            
        Raises:
            ServiceError: If parameter retrieval fails
        """
        try:
            logger.info("Getting supported quantum chemistry parameters")
            parameters = get_all_supported_parameters()
            logger.info("Successfully retrieved supported parameters")
            return parameters
        except Exception as e:
            logger.error(f"Error getting supported parameters: {e}", exc_info=True)
            raise ServiceError(f'Failed to retrieve supported parameters: {str(e)}')
    
    def validate_calculation_parameters(self, params: Dict[str, Any]) -> Optional[str]:
        """
        Validate calculation parameters for compatibility and theoretical correctness.
        
        Args:
            params: Dictionary of calculation parameters
            
        Returns:
            None if validation passes, error message string if validation fails
        """
        calculation_method = params.get('calculation_method')
        
        # Check HF method specific constraints
        if calculation_method == 'HF':
            if params.get('exchange_correlation') and params.get('exchange_correlation') != 'B3LYP':
                logger.warning(f"HF calculation with non-default exchange_correlation: {params.get('exchange_correlation')}")
            
            if params.get('tddft_nstates') and params.get('tddft_nstates') > 10:
                logger.warning(f"TDDFT parameters specified for HF calculation - these will be ignored")
        
        # Check DFT method constraints
        if calculation_method == 'DFT':
            if not params.get('exchange_correlation'):
                return "DFT method requires an exchange-correlation functional to be specified"
        
        # Check TDDFT method constraints
        if calculation_method == 'TDDFT':
            if not params.get('exchange_correlation'):
                return "TDDFT method requires an exchange-correlation functional to be specified"
            if not params.get('tddft_nstates') or params.get('tddft_nstates') < 1:
                return "TDDFT method requires tddft_nstates to be specified and greater than 0"
        
        # Check CASCI/CASSCF method constraints
        if calculation_method in ['CASCI', 'CASSCF']:
            if not params.get('ncas') or params.get('ncas') < 1:
                return f"{calculation_method} method requires ncas (active space orbitals) to be specified and greater than 0"
            if not params.get('nelecas') or params.get('nelecas') < 1:
                return f"{calculation_method} method requires nelecas (active space electrons) to be specified and greater than 0"
            if params.get('nelecas') > 2 * params.get('ncas'):
                return f"{calculation_method} method: nelecas ({params.get('nelecas')}) cannot exceed 2 * ncas ({2 * params.get('ncas')})"
        
        # Check general parameter constraints
        if params.get('spin', 0) < 0:
            return "Spin multiplicity (2S) cannot be negative"
        
        if params.get('charges') is not None and abs(params.get('charges')) > 10:
            logger.warning(f"High molecular charge ({params.get('charges')}) - please verify this is correct")
        
        return None
    
    def start_calculation(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Start a new quantum chemistry calculation.
        
        Args:
            params: Dictionary of calculation parameters
            
        Returns:
            Dict containing calculation instance information
            
        Raises:
            ValidationError: If parameters are invalid
            ResourceUnavailableError: If process manager is unavailable
            InsufficientResourcesError: If insufficient system resources
            ServiceError: For other errors
        """
        try:
            # Validate calculation parameters
            validation_error = self.validate_calculation_parameters(params)
            if validation_error:
                logger.warning(f"Parameter validation failed: {validation_error}")
                raise ValidationError(f'Invalid parameters: {validation_error}')
            
            # Initialize calculation directory
            try:
                calc_dir = self.file_manager.create_calculation_dir(params['name'])
                calculation_id = os.path.basename(calc_dir)
                
                # Save initial parameters
                self.file_manager.save_calculation_parameters(calc_dir, params)
                logger.info(f"Created calculation directory and saved parameters for calculation {calculation_id}")
            except Exception as file_error:
                logger.error(f"Failed to set up calculation files: {file_error}")
                raise ServiceError(f'Failed to initialize calculation: {str(file_error)}')
            
            # Submit to process manager
            try:
                process_manager = get_process_manager()
                
                # Apply current settings
                try:
                    from quantum_calc import update_process_manager_settings
                    update_process_manager_settings()
                except Exception as settings_error:
                    logger.warning(f"Failed to update process manager settings: {settings_error}")
                
                # Submit calculation
                logger.info(f"About to submit calculation {calculation_id}")
                success, initial_status, waiting_reason = process_manager.submit_calculation(calculation_id, params)
                logger.info(f"Submit result: success={success}, status={initial_status}, reason={waiting_reason}")
                
                # Handle submission failure
                if not success:
                    error_message = waiting_reason if waiting_reason else 'Failed to submit calculation to process pool.'
                    self.file_manager.save_calculation_status(calc_dir, 'error')
                    self.file_manager.save_calculation_results(calc_dir, {'error': error_message})
                    logger.error(f"Failed to submit calculation {calculation_id} to process pool: {error_message}")
                    
                    # Create error instance but still return it (calculation was created)
                    error_instance = self._build_calculation_instance(calculation_id, params, 'error')
                    error_instance['error'] = error_message
                    return error_instance
                
                # Set initial status
                self.file_manager.save_calculation_status(calc_dir, initial_status, waiting_reason)
                logger.info(f"Queued calculation {calculation_id} for molecule '{params['name']}'")
                
                # Return initial calculation instance
                return self._build_calculation_instance(calculation_id, params, initial_status, waiting_reason)
                
            except ProcessManagerError as e:
                # Clean up created directory
                try:
                    shutil.rmtree(calc_dir, ignore_errors=True)
                except Exception:
                    pass
                logger.error(f"Process manager error: {e}")
                raise ResourceUnavailableError(f'System initialization error: Unable to initialize calculation system. Please check system resources and try again.')
            except Exception as submit_error:
                # Update status to error
                self.file_manager.save_calculation_status(calc_dir, 'error')
                error_message = f'Failed to submit calculation: {str(submit_error)}'
                self.file_manager.save_calculation_results(calc_dir, {'error': error_message})
                logger.error(f"Unexpected error during calculation submission: {submit_error}")
                raise ServiceError(error_message)
                
        except (InputError, GeometryError) as e:
            logger.warning(f"Invalid calculation parameters: {e}")
            raise ValidationError(f'Invalid input parameters: {str(e)}')
        except FileManagerError as e:
            logger.error(f"File management error during calculation setup: {e}")
            raise ServiceError('Failed to set up calculation files.')
        except OSError as e:
            logger.error(f"System error during calculation setup: {e}")
            raise InsufficientResourcesError('Insufficient system resources to start calculation.')
        except PermissionError as e:
            logger.error(f"Permission error during calculation setup: {e}")
            raise PermissionDeniedError('System permission error. Please contact administrator.')
    
    def list_calculations(
        self,
        name_query: Optional[str] = None,
        status: Optional[str] = None,
        calculation_method: Optional[str] = None,
        basis_function: Optional[str] = None,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        List available calculations with optional filtering.

        Args:
            name_query: Partial match search in calculation name (case-insensitive)
            status: Filter by status ("completed", "running", "error", etc.)
            calculation_method: Filter by calculation method ("DFT", "HF", "MP2", etc.)
            basis_function: Filter by basis set (case-insensitive)
            date_from: Start date for date range filtering (ISO format: YYYY-MM-DD)
            date_to: End date for date range filtering (ISO format: YYYY-MM-DD)

        Returns:
            Dict containing list of calculations and metadata

        Raises:
            ServiceError: If listing fails
        """
        try:
            calculations = self.file_manager.list_calculations(
                name_query=name_query,
                status=status,
                calculation_method=calculation_method,
                basis_function=basis_function,
                date_from=date_from,
                date_to=date_to
            )

            return {
                'base_directory': self.file_manager.get_base_directory(),
                'calculations': calculations,
                'count': len(calculations)
            }
        except FileManagerError as e:
            logger.error(f"File manager error while listing calculations: {e}")
            raise ServiceError('Unable to access calculation directory.')
        except PermissionError as e:
            logger.error(f"Permission error while listing calculations: {e}")
            raise PermissionDeniedError('Permission denied accessing calculation directory.')
        except OSError as e:
            logger.error(f"System error while listing calculations: {e}")
            raise ServiceError('System error accessing calculation files.')
    
    def get_calculation_details(self, calculation_id: str) -> Dict[str, Any]:
        """
        Get detailed information about a specific calculation.
        
        Args:
            calculation_id: Unique ID of the calculation
            
        Returns:
            Dict containing calculation details
            
        Raises:
            NotFoundError: If calculation not found
            ValidationError: If calculation_id is invalid
            ServiceError: For other errors
        """
        try:
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            # Read all calculation data from disk
            parameters = self.file_manager.read_calculation_parameters(calc_path) or {}
            results = self.file_manager.read_calculation_results(calc_path)
            status, waiting_reason = self.file_manager.read_calculation_status_details(calc_path)
            
            display_name = self.file_manager._get_display_name(calculation_id, parameters)
            creation_date = parameters.get('created_at', datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat())
            
            calculation_instance = {
                'id': calculation_id,
                'name': display_name,
                'status': status,
                'createdAt': creation_date,
                'updatedAt': datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat(),
                'workingDirectory': calc_path,
                'parameters': parameters,
                'results': results
            }
            
            # Include waiting reason if available
            if waiting_reason is not None:
                calculation_instance['waitingReason'] = waiting_reason
            
            return {
                'calculation': calculation_instance,
                'files': {
                    'checkpoint_exists': self.file_manager.file_exists(calc_path, 'calculation.chk'),
                    'parameters_file_exists': parameters is not None,
                    'results_file_exists': results is not None,
                }
            }
        except FileManagerError as e:
            logger.error(f"File manager error getting calculation {calculation_id}: {e}")
            raise ServiceError('Unable to access calculation data.')
        except OSError as e:
            logger.error(f"System error getting calculation {calculation_id}: {e}")
            raise ServiceError('System error accessing calculation files.')
        except ValueError as e:
            logger.warning(f"Invalid calculation ID {calculation_id}: {e}")
            raise ValidationError('Invalid calculation ID format.')
    
    def update_calculation(self, calculation_id: str, new_name: str) -> Dict[str, Any]:
        """
        Update calculation metadata (currently only name).
        
        Args:
            calculation_id: Unique ID of the calculation
            new_name: New display name for the calculation
            
        Returns:
            Dict containing success message and updated name
            
        Raises:
            NotFoundError: If calculation not found
            ValidationError: If inputs are invalid
            ServiceError: For other errors
        """
        try:
            result_id = self.file_manager.rename_calculation(calculation_id, new_name)
            if not result_id:
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            logger.info(f"Updated display name for calculation {calculation_id} to '{new_name}'")
            return {
                'message': 'Calculation renamed successfully.',
                'name': new_name
            }
        except FileManagerError as e:
            logger.error(f"File manager error updating calculation {calculation_id}: {e}")
            raise ServiceError('Unable to update calculation data.')
        except ValueError as e:
            logger.warning(f"Invalid input for calculation {calculation_id}: {e}")
            raise ValidationError('Invalid input data.')
        except OSError as e:
            logger.error(f"System error updating calculation {calculation_id}: {e}")
            raise ServiceError('System error updating calculation.')
    
    def cancel_calculation(self, calculation_id: str) -> Dict[str, Any]:
        """
        Cancel a running calculation.
        
        Args:
            calculation_id: Unique ID of the calculation
            
        Returns:
            Dict containing success message
            
        Raises:
            ValidationError: If calculation is not running
            ResourceUnavailableError: If process manager is unavailable
            ServiceError: For other errors
        """
        try:
            process_manager = get_process_manager()
            
            # Check if calculation is running
            if not process_manager.is_running(calculation_id):
                raise ValidationError(f'Calculation "{calculation_id}" is not currently running.')
            
            # Attempt to cancel
            cancelled = process_manager.cancel_calculation(calculation_id)
            
            if cancelled:
                logger.info(f"Successfully cancelled calculation {calculation_id}")
                return {
                    'message': f'Calculation "{calculation_id}" has been cancelled successfully',
                    'calculation_id': calculation_id
                }
            else:
                raise ValidationError(f'Failed to cancel calculation "{calculation_id}". It may have already completed.')
                
        except ProcessManagerError as e:
            logger.error(f"Process manager error cancelling calculation {calculation_id}: {e}")
            raise ResourceUnavailableError('Process manager unavailable.')
        except ValueError as e:
            logger.warning(f"Invalid calculation ID for cancellation {calculation_id}: {e}")
            raise ValidationError('Invalid calculation ID.')
    
    def delete_calculation(self, calculation_id: str) -> Dict[str, Any]:
        """
        Delete a calculation and its files.
        
        Args:
            calculation_id: Unique ID of the calculation
            
        Returns:
            Dict containing success message
            
        Raises:
            NotFoundError: If calculation not found
            ValidationError: If calculation_id is invalid
            ServiceError: For other errors
        """
        try:
            process_manager = get_process_manager()
            
            # Try to cancel if running
            if process_manager.is_running(calculation_id):
                logger.info(f"Cancelling running calculation {calculation_id} before deletion")
                process_manager.cancel_calculation(calculation_id)
                time.sleep(0.5)  # Wait for cancellation
            
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            # Delete the directory
            shutil.rmtree(calc_path)
            logger.info(f"Deleted calculation directory: {calc_path}")
            
            return {
                'message': f'Calculation "{calculation_id}" has been deleted successfully',
                'deleted_id': calculation_id
            }
            
        except ProcessManagerError as e:
            logger.error(f"Process manager error deleting calculation {calculation_id}: {e}")
            raise ResourceUnavailableError('Process manager unavailable.')
        except FileManagerError as e:
            logger.error(f"File manager error deleting calculation {calculation_id}: {e}")
            raise ServiceError('Unable to access calculation files.')
        except PermissionError as e:
            logger.error(f"Permission error deleting calculation {calculation_id}: {e}")
            raise PermissionDeniedError('Permission denied deleting calculation.')
        except OSError as e:
            logger.error(f"System error deleting calculation {calculation_id}: {e}")
            raise ServiceError('System error deleting calculation files.')
        except ValueError as e:
            logger.warning(f"Invalid calculation ID for deletion {calculation_id}: {e}")
            raise ValidationError('Invalid calculation ID.')
    
    def get_calculation_status(self) -> Dict[str, Any]:
        """
        Get status information about the calculation system.
        
        Returns:
            Dict containing process pool and system information
            
        Raises:
            ResourceUnavailableError: If process manager is unavailable
            ServiceError: For other errors
        """
        try:
            process_manager = get_process_manager()
            active_calculations = process_manager.get_active_calculations()
            
            return {
                'process_pool': {
                    'max_workers': process_manager.max_workers,
                    'active_calculations': active_calculations,
                    'active_count': len(active_calculations),
                    'is_shutdown': process_manager._shutdown
                },
                'system': {
                    'cpu_count': os.cpu_count()
                }
            }
        except ProcessManagerError as e:
            logger.error(f"Process manager error getting status: {e}")
            raise ResourceUnavailableError('Process manager is unavailable.')
        except AttributeError as e:
            logger.error(f"Process manager attribute error: {e}")
            raise ServiceError('Process manager not properly initialized.')
    
    def get_molecular_orbitals(self, calculation_id: str) -> Dict[str, Any]:
        """
        Get molecular orbital information for a calculation.
        
        Args:
            calculation_id: Unique ID of the calculation
            
        Returns:
            Dict containing orbital information
            
        Raises:
            NotFoundError: If calculation not found or orbital data unavailable
            ValidationError: If calculation is not completed
            ServiceError: For other errors
        """
        try:
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            # Check if completed
            status = self.file_manager.read_calculation_status(calc_path)
            if status != 'completed':
                raise ValidationError(f'Calculation "{calculation_id}" is not completed. Status: {status}')
            
            # Initialize orbital generator
            orbital_generator = MolecularOrbitalGenerator(calc_path)
            
            # Validate calculation data
            if not orbital_generator.validate_calculation():
                raise NotFoundError('Orbital data is not available or calculation data is invalid.')
            
            # Get orbital summary
            orbital_summary = orbital_generator.get_orbital_summary()
            
            logger.info(f"Retrieved orbital information for calculation {calculation_id}")
            logger.info(f"Total orbitals: {orbital_summary['total_orbitals']}, "
                       f"HOMO: {orbital_summary['homo_index']}, LUMO: {orbital_summary['lumo_index']}")
            
            return orbital_summary
            
        except CalculationError as e:
            logger.error(f"Calculation error getting orbitals for {calculation_id}: {e}")
            raise ServiceError(str(e))
        except FileManagerError as e:
            logger.error(f"File manager error getting orbitals for {calculation_id}: {e}")
            raise ServiceError('Unable to access calculation files.')
    
    def generate_orbital_cube(
        self,
        calculation_id: str,
        orbital_index: int,
        grid_size: int = 80,
        isovalue_pos: Optional[float] = None,
        isovalue_neg: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        Generate CUBE file for a specific molecular orbital.
        
        Args:
            calculation_id: Unique ID of the calculation
            orbital_index: Index of the orbital
            grid_size: Grid resolution
            isovalue_pos: Positive isovalue
            isovalue_neg: Negative isovalue
            
        Returns:
            Dict containing CUBE file data and metadata
            
        Raises:
            NotFoundError: If calculation not found or orbital data unavailable
            ValidationError: If parameters are invalid or calculation not completed
            ServiceError: For other errors
        """
        try:
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            # Check if completed
            status = self.file_manager.read_calculation_status(calc_path)
            if status != 'completed':
                raise ValidationError(f'Calculation "{calculation_id}" is not completed. Status: {status}')
            
            # Initialize orbital generator
            orbital_generator = MolecularOrbitalGenerator(calc_path)
            
            # Validate calculation data
            if not orbital_generator.validate_calculation():
                raise NotFoundError('Orbital data is not available or calculation data is invalid.')
            
            # Generate CUBE file
            logger.info(f"Generating CUBE file for calculation {calculation_id}, orbital {orbital_index}")
            logger.info(f"Parameters: grid_size={grid_size}, isovalue_pos={isovalue_pos}, isovalue_neg={isovalue_neg}")
            
            cube_data = orbital_generator.generate_cube_file(
                orbital_index=orbital_index,
                grid_size=grid_size,
                isovalue_pos=isovalue_pos,
                isovalue_neg=isovalue_neg,
                return_content=True,
                save_to_disk=True
            )
            
            if cube_data.get('cached', False):
                logger.info(f"Using cached CUBE file for calculation {calculation_id}, orbital {orbital_index}")
            else:
                logger.info(f"Successfully generated CUBE file for calculation {calculation_id}, orbital {orbital_index}")
            logger.info(f"File size: {cube_data['generation_params']['file_size_kb']:.1f} KB")
            
            return cube_data
            
        except CalculationError as e:
            logger.error(f"Calculation error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}")
            raise ServiceError(str(e))
        except FileManagerError as e:
            logger.error(f"File manager error generating CUBE for {calculation_id}, orbital {orbital_index}: {e}")
            raise ServiceError('Unable to access calculation files.')
        except ValueError as e:
            logger.warning(f"Invalid orbital index for calculation {calculation_id}: {e}")
            raise ValidationError('Invalid orbital index.')
    
    def list_cube_files(self, calculation_id: str) -> Dict[str, Any]:
        """
        List all CUBE files for a calculation.
        
        Args:
            calculation_id: Unique ID of the calculation
            
        Returns:
            Dict containing list of CUBE files and metadata
            
        Raises:
            NotFoundError: If calculation not found
            ServiceError: For other errors
        """
        try:
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            cube_files = self.file_manager.get_cube_files_info(calc_path)
            
            logger.info(f"Found {len(cube_files)} CUBE files for calculation {calculation_id}")
            
            return {
                'calculation_id': calculation_id,
                'cube_files': cube_files,
                'total_files': len(cube_files),
                'total_size_kb': sum(f['file_size_kb'] for f in cube_files)
            }
        except Exception as e:
            logger.error(f"Error listing CUBE files for {calculation_id}: {e}", exc_info=True)
            raise ServiceError('An internal error occurred.')
    
    def delete_cube_files(self, calculation_id: str, orbital_index: Optional[int] = None) -> Dict[str, Any]:
        """
        Delete CUBE files for a calculation.
        
        Args:
            calculation_id: Unique ID of the calculation
            orbital_index: Optional orbital index to delete specific file
            
        Returns:
            Dict containing deletion results
            
        Raises:
            NotFoundError: If calculation not found
            ServiceError: For other errors
        """
        try:
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            deleted_count = self.file_manager.delete_cube_files(calc_path, orbital_index)
            
            if deleted_count > 0:
                if orbital_index is not None:
                    logger.info(f"Deleted {deleted_count} CUBE files for orbital {orbital_index} in calculation {calculation_id}")
                    message = f"Deleted {deleted_count} CUBE files for orbital {orbital_index}."
                else:
                    logger.info(f"Deleted {deleted_count} CUBE files for calculation {calculation_id}")
                    message = f"Deleted {deleted_count} CUBE files."
            else:
                if orbital_index is not None:
                    message = f"No CUBE files found for orbital {orbital_index}."
                else:
                    message = "No CUBE files found."
            
            return {
                'calculation_id': calculation_id,
                'orbital_index': orbital_index,
                'deleted_files': deleted_count,
                'message': message
            }
        except Exception as e:
            logger.error(f"Error deleting CUBE files for {calculation_id}: {e}", exc_info=True)
            raise ServiceError('An internal error occurred.')
    
    def generate_ir_spectrum(
        self,
        calculation_id: str,
        broadening_fwhm: float = 100.0,
        x_min: float = 400.0,
        x_max: float = 4000.0,
        show_peaks: bool = True
    ) -> Dict[str, Any]:
        """
        Generate IR spectrum for a calculation.
        
        Args:
            calculation_id: Unique ID of the calculation
            broadening_fwhm: Full width at half maximum for broadening
            x_min: Minimum wavenumber
            x_max: Maximum wavenumber
            show_peaks: Whether to mark peaks
            
        Returns:
            Dict containing spectrum data and plot
            
        Raises:
            NotFoundError: If calculation not found or frequency data unavailable
            ValidationError: If parameters are invalid or calculation not completed
            ServiceError: For other errors
        """
        try:
            calc_path = os.path.join(self.file_manager.get_base_directory(), calculation_id)
            
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')
            
            # Check if completed
            status = self.file_manager.read_calculation_status(calc_path)
            if status != 'completed':
                raise ValidationError(f'Calculation is not completed (status: {status}). IR spectrum cannot be generated.')
            
            # Read results
            results = self.file_manager.read_calculation_results(calc_path)
            if not results:
                raise NotFoundError('Calculation results not found.')
            
            # Check for vibrational frequency data
            if not results.get('vibrational_frequencies'):
                raise ValidationError('No vibrational frequency data found. Frequency analysis may not have been performed or failed.')
            
            # Validate parameters
            if broadening_fwhm <= 0:
                raise ValidationError('Broadening FWHM must be positive.')
            
            if x_min >= x_max:
                raise ValidationError('x_min must be less than x_max.')
            
            logger.info(f"Generating IR spectrum for calculation {calculation_id}")
            logger.info(f"Parameters: FWHM={broadening_fwhm}, range=({x_min}, {x_max}), show_peaks={show_peaks}")
            
            # Generate spectrum
            ir_result = create_ir_spectrum_from_calculation_results(
                results,
                broadening_fwhm=broadening_fwhm,
                x_range=(x_min, x_max)
            )
            
            if not ir_result.get('success'):
                error_msg = ir_result.get('error', 'Unknown error occurred during IR spectrum generation')
                logger.error(f"IR spectrum generation failed for {calculation_id}: {error_msg}")
                raise ServiceError(f'IR spectrum generation failed: {error_msg}')
            
            spectrum_data = ir_result.get('spectrum_data', {})
            plot_image = ir_result.get('plot_image_base64')
            
            logger.info(f"IR spectrum generated successfully for {calculation_id}")
            logger.info(f"Spectrum contains {len(spectrum_data.get('peaks', []))} peaks")
            
            return {
                'calculation_id': calculation_id,
                'spectrum': {
                    'x_axis': spectrum_data.get('x_axis', []),
                    'y_axis': spectrum_data.get('spectrum', []),
                    'peaks': spectrum_data.get('peaks', []),
                    'metadata': spectrum_data.get('metadata', {})
                },
                'plot_image_base64': plot_image,
                'generation_info': {
                    'broadening_fwhm_cm': broadening_fwhm,
                    'frequency_range_cm': [x_min, x_max],
                    'peaks_marked': show_peaks,
                    'generated_at': datetime.now().isoformat()
                }
            }
        except ValueError as e:
            logger.error(f"Invalid parameters for IR spectrum generation ({calculation_id}): {e}")
            raise ValidationError(f'Invalid parameters: {str(e)}')
    
    def _build_calculation_instance(
        self,
        calculation_id: str,
        parameters: Dict[str, Any],
        status: str,
        waiting_reason: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Build calculation instance dict for responses.
        
        Args:
            calculation_id: Unique ID
            parameters: Calculation parameters
            status: Current status
            waiting_reason: Optional waiting reason
            
        Returns:
            Dict containing calculation instance
        """
        instance = {
            'id': calculation_id,
            'name': parameters['name'],
            'status': status,
            'createdAt': parameters['created_at'],
            'updatedAt': parameters['created_at'],
            'parameters': parameters
        }
        
        if waiting_reason is not None:
            instance['waitingReason'] = waiting_reason
        
        return instance
