"""
Analysis Tools Module

This module provides tools for reading, analyzing, and retrieving data
without modifying system state. These tools are safe for all agents to use,
especially those focused on analysis and reporting.

These tools include calculation data retrieval, molecular orbital analysis,
IR spectrum generation, and system status queries.
"""

import json
import logging
from typing import Optional
from datetime import datetime

from services import (
    get_quantum_service, get_pubchem_service, get_settings_service,
    get_system_service, ServiceError
)

# Logger configuration
logger = logging.getLogger(__name__)


def _validation_error(message: str) -> str:
    """
    Return a validation error message in JSON format.

    This helper function ensures consistent error response format for input validation
    errors by returning a JSON string with success=False and the error message.

    Args:
        message: The validation error message

    Returns:
        str: JSON-formatted error message with structure: {"success": false, "error": "message"}
    """
    return json.dumps({
        "success": False,
        "error": message
    }, ensure_ascii=False, indent=2)


def _handle_service_error(error: ServiceError, context: str = "") -> str:
    """
    Handle service errors and return appropriate error messages in JSON format.

    This function ensures consistent error response format for AI agents by always
    returning a JSON string with a structured error object containing success=False
    and an error message.

    Args:
        error: The ServiceError that occurred
        context: Additional context about where the error occurred

    Returns:
        str: JSON-formatted error message with structure: {"success": false, "error": "message"}
    """
    context_prefix = f"{context}: " if context else ""
    error_msg = f"{context_prefix}{error.message}"

    logger.error(f"Service error in {context}: {error_msg}")

    # Return error in consistent JSON format
    return json.dumps({
        "success": False,
        "error": error_msg
    }, ensure_ascii=False, indent=2)


def list_all_calculations() -> str:
    """
    Retrieves a list of all quantum chemistry calculations that have been executed and saved.

    Returns a list of calculation summaries including ID, name, status, and creation date.
    Use this when the user asks "show me calculation history" or "what calculations are available".

    Returns:
        str: JSON string of calculation list. On success, contains an array of calculation summaries.
             On error, contains structured text with error message.
    """
    try:
        logger.debug("Fetching all calculations from service")
        quantum_service = get_quantum_service()
        result = quantum_service.list_calculations()

        # Return in the expected format with success flag
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "list_all_calculations")
    except Exception as e:
        logger.error(f"Unexpected error in list_all_calculations: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def find_calculations(
    name_query: Optional[str] = None,
    status: Optional[str] = None,
    calculation_method: Optional[str] = None,
    basis_function: Optional[str] = None,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None
) -> str:
    """
    Search for calculations using flexible filtering criteria.

    This function searches through all saved calculations and filters them based on
    the provided criteria. All parameters are optional - you can combine them as needed.
    Performs case-insensitive partial matching for name_query.

    Args:
        name_query (str, optional): Partial match search in calculation name.
                                   Case-insensitive. Example: "水" matches "水分子のDFT計算"
        status (str, optional): Filter by status. Options: "completed", "running", "error"
        calculation_method (str, optional): Filter by calculation method.
                                           Options: "DFT", "HF", "MP2", "CCSD", "TDDFT", "CASCI", "CASSCF"
        basis_function (str, optional): Filter by basis set. Case-insensitive.
                                       Examples: "6-31G(d)", "cc-pVDZ", "def2-SVP"
        date_from (str, optional): Start date for date range filtering (ISO format: YYYY-MM-DD)
        date_to (str, optional): End date for date range filtering (ISO format: YYYY-MM-DD)

    Returns:
        str: JSON string containing search results.

             Single match:
             {
               "success": true,
               "count": 1,
               "calculation": {...}  // Detailed calculation info
             }

             Multiple matches:
             {
               "success": true,
               "multiple_matches": true,
               "count": N,
               "message": "N件の計算が見つかりました。いずれかを選択してください:",
               "calculations": [
                 {
                   "id": "calc_...",
                   "name": "...",
                   "status": "...",
                   "calculation_method": "...",
                   "basis_function": "...",
                   "created_at": "..."
                 },
                 ...
               ]
             }

             No matches:
             {
               "success": false,
               "error": "検索条件に一致する計算が見つかりませんでした。"
             }
    """
    try:
        # Get all calculations
        logger.debug("Fetching all calculations for filtering")
        all_calculations_json = list_all_calculations()
        all_calculations_result = json.loads(all_calculations_json)

        if not all_calculations_result.get('success'):
            return all_calculations_json  # Return error from base function

        # Extract calculations list from the nested data structure
        # Structure: {"success": true, "data": {"calculations": [...], "count": N, "base_directory": "..."}}
        data = all_calculations_result.get('data', {})
        calculations = data.get('calculations', [])
        logger.info(f"Retrieved {len(calculations)} total calculations for filtering")

        # Apply filters
        filtered = calculations

        # Filter by name (case-insensitive partial match)
        if name_query:
            name_lower = name_query.lower()
            filtered = [
                calc for calc in filtered
                if name_lower in calc.get('name', '').lower()
            ]
            logger.debug(f"After name filter '{name_query}': {len(filtered)} calculations")

        # Filter by status (exact match)
        if status:
            filtered = [
                calc for calc in filtered
                if calc.get('status', '').lower() == status.lower()
            ]
            logger.debug(f"After status filter '{status}': {len(filtered)} calculations")

        # Filter by calculation_method (case-insensitive exact match)
        if calculation_method:
            filtered = [
                calc for calc in filtered
                if calc.get('calculation_method', '').lower() == calculation_method.lower()
            ]
            logger.debug(f"After method filter '{calculation_method}': {len(filtered)} calculations")

        # Filter by basis_function (case-insensitive exact match)
        if basis_function:
            filtered = [
                calc for calc in filtered
                if calc.get('basis_function', '').lower() == basis_function.lower()
            ]
            logger.debug(f"After basis filter '{basis_function}': {len(filtered)} calculations")

        # Filter by date range
        if date_from or date_to:
            def parse_date(date_str: str) -> Optional[datetime]:
                """Parse ISO date string to datetime object."""
                try:
                    return datetime.fromisoformat(date_str.replace('Z', '+00:00'))
                except (ValueError, AttributeError):
                    return None

            date_from_dt = parse_date(date_from) if date_from else None
            date_to_dt = parse_date(date_to) if date_to else None

            filtered_by_date = []
            for calc in filtered:
                created_at = calc.get('created_at')
                if not created_at:
                    continue

                calc_date = parse_date(created_at)
                if not calc_date:
                    continue

                # Check date range
                if date_from_dt and calc_date < date_from_dt:
                    continue
                if date_to_dt and calc_date > date_to_dt:
                    continue

                filtered_by_date.append(calc)

            filtered = filtered_by_date
            logger.debug(f"After date filter (from={date_from}, to={date_to}): {len(filtered)} calculations")

        # Process results
        count = len(filtered)

        if count == 0:
            # No matches
            logger.info("No calculations matched the search criteria")
            return json.dumps({
                "success": False,
                "error": "検索条件に一致する計算が見つかりませんでした。条件を変更して再度お試しください。"
            }, ensure_ascii=False, indent=2)

        elif count == 1:
            # Single match - return detailed info
            logger.info(f"Found single matching calculation: {filtered[0].get('id')}")
            return json.dumps({
                "success": True,
                "count": 1,
                "calculation": filtered[0]
            }, ensure_ascii=False, indent=2)

        else:
            # Multiple matches - return list for user selection
            logger.info(f"Found {count} matching calculations")

            # Sort by created_at (newest first)
            sorted_calcs = sorted(
                filtered,
                key=lambda x: x.get('created_at', ''),
                reverse=True
            )

            # Prepare simplified list for presentation
            calc_list = []
            for calc in sorted_calcs:
                calc_list.append({
                    "id": calc.get('id'),
                    "name": calc.get('name'),
                    "status": calc.get('status'),
                    "calculation_method": calc.get('calculation_method'),
                    "basis_function": calc.get('basis_function'),
                    "created_at": calc.get('created_at')
                })

            return json.dumps({
                "success": True,
                "multiple_matches": True,
                "count": count,
                "message": f"{count}件の計算が見つかりました。いずれかを選択してください:",
                "calculations": calc_list
            }, ensure_ascii=False, indent=2)

    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error in find_calculations: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"計算データの解析中にエラーが発生しました: {str(e)}"
        }, ensure_ascii=False, indent=2)

    except Exception as e:
        logger.error(f"Unexpected error in find_calculations: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"計算検索中に予期しないエラーが発生しました: {str(e)}"
        }, ensure_ascii=False, indent=2)


def get_calculation_details(calculation_id: str) -> str:
    """
    Retrieves detailed results for a completed calculation with the specified ID.

    Returns all physical properties and analysis results obtained from the calculation,
    including energy, dipole moment, orbital energies, optimized structures, etc.

    Args:
        calculation_id (str): Unique ID of the calculation to retrieve. Required parameter.
                             Example: "calc_20240101_120000_abcd1234"

    Returns:
        str: JSON string of calculation details. On success, contains complete calculation result data.
             On error, contains structured text with error message.
    """
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    try:
        logger.debug(f"Fetching calculation details for ID: {calculation_id}")
        quantum_service = get_quantum_service()
        result = quantum_service.get_calculation_details(calculation_id)

        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"get_calculation_details(id={calculation_id})")
    except Exception as e:
        logger.error(f"Unexpected error in get_calculation_details: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def get_supported_parameters() -> str:
    """
    Retrieves lists of calculation methods, basis sets, functionals, and solvents available in the application.

    Provides information useful for users when selecting calculation settings.
    Includes calculation methods such as DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF,
    and lists of available basis sets and exchange-correlation functionals for each method.

    Returns:
        str: JSON string of supported parameters. On success, contains parameter lists.
             On error, contains structured text with error message.
    """
    try:
        logger.debug("Fetching supported parameters from service")
        quantum_service = get_quantum_service()
        parameters = quantum_service.get_supported_parameters()

        return json.dumps({
            'success': True,
            'data': parameters
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "get_supported_parameters")
    except Exception as e:
        logger.error(f"Unexpected error in get_supported_parameters: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def get_app_settings() -> str:
    """
    Retrieve current application settings.

    This function returns the current configuration settings for the application,
    including parallel processing limits, resource utilization constraints, and API keys.

    Returns:
        str: JSON string containing application settings.
             On success, includes:
             - max_parallel_instances: Maximum number of parallel calculations
             - max_cpu_utilization_percent: Maximum CPU usage percentage
             - max_memory_utilization_percent: Maximum memory usage percentage
             - system_total_cores: Total CPU cores in the system
             - system_total_memory_mb: Total system memory in MB
             - gemini_api_key: Google Gemini API key (if configured)
             On error, contains structured error message.
    """
    try:
        logger.debug("Fetching application settings")
        settings_service = get_settings_service()
        settings = settings_service.get_settings()

        logger.info("Successfully retrieved application settings")
        return json.dumps({
            'success': True,
            'data': {'settings': settings}
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "get_app_settings")
    except Exception as e:
        logger.error(f"Unexpected error in get_app_settings: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def get_system_resources() -> str:
    """
    Retrieve current system resource status and allocation information.

    This function provides detailed information about system resources including
    CPU and memory usage, resource constraints from settings, and current allocation
    to running calculations. Useful for monitoring system health and understanding
    why new calculations might be queued.

    Returns:
        str: JSON string containing system resource information.
             On success, includes:
             - system_info: Total and available CPU/memory, current usage percentages
             - resource_constraints: Maximum allowed resource utilization from settings
             - allocated_resources: Resources currently allocated to active calculations
             On error, contains structured error message.
    """
    try:
        logger.debug("Fetching system resource status")
        system_service = get_system_service()
        result = system_service.get_resource_status()

        logger.info("Successfully retrieved system resource status")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "get_system_resources")
    except Exception as e:
        logger.error(f"Unexpected error in get_system_resources: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def get_molecular_orbitals(calculation_id: str) -> str:
    """
    Retrieve molecular orbital information for a completed calculation.

    This function returns detailed information about all molecular orbitals including
    their energies, occupancies, and types (HOMO, LUMO, core, virtual).
    This data is useful for understanding electronic structure and preparing orbital visualizations.

    Args:
        calculation_id (str): Unique ID of the calculation.
                             Example: "calc_20240101_120000_abcd1234"

    Returns:
        str: JSON string containing orbital information.
             On success, includes:
             - List of all orbitals with energies (Hartree and eV), occupancy, and labels
             - HOMO and LUMO indices
             - Total number of occupied and virtual orbitals
             On error, contains structured error message.
    """
    # Input validation
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    if len(calculation_id.strip()) == 0:
        return _validation_error("calculation_id cannot be empty or contain only whitespace.")

    try:
        logger.debug(f"Fetching molecular orbital information for calculation: {calculation_id}")
        quantum_service = get_quantum_service()
        orbital_summary = quantum_service.get_molecular_orbitals(calculation_id)

        num_orbitals = orbital_summary.get('total_orbitals', 0)
        logger.info(f"Successfully retrieved {num_orbitals} molecular orbitals for calculation {calculation_id}")
        return json.dumps({
            'success': True,
            'data': orbital_summary
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"get_molecular_orbitals(id={calculation_id})")
    except Exception as e:
        logger.error(f"Unexpected error in get_molecular_orbitals: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def generate_orbital_cube(
    calculation_id: str,
    orbital_index: int,
    grid_size: int = 80,
    isovalue_pos: float = 0.02,
    isovalue_neg: float = -0.02
) -> str:
    """
    Generate a CUBE file for visualizing a specific molecular orbital.

    CUBE files contain 3D volumetric data for orbital visualization in molecular viewers.
    This function generates the CUBE file data for a specified orbital with customizable
    grid resolution and isosurface values.

    Args:
        calculation_id (str): Unique ID of the calculation.
                             Example: "calc_20240101_120000_abcd1234"
        orbital_index (int): Index of the molecular orbital to visualize (0-based).
                            Use get_molecular_orbitals() to find available orbital indices.
        grid_size (int): Grid resolution for the CUBE file (40-120).
                        Higher values give better resolution but larger files.
                        Default: 80
        isovalue_pos (float): Positive isovalue for orbital surface (0.001-0.1).
                             Default: 0.02
        isovalue_neg (float): Negative isovalue for orbital surface (-0.1 to -0.001).
                             Default: -0.02

    Returns:
        str: JSON string containing CUBE file data and metadata.
             On success, includes:
             - cube_data: CUBE file content as string
             - orbital_info: Information about the orbital
             - generation_params: Parameters used for generation
             - file_path: Path to saved file (if cached)
             On error, contains structured error message.
    """
    # Input validation
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    if len(calculation_id.strip()) == 0:
        return _validation_error("calculation_id cannot be empty or contain only whitespace.")

    if not isinstance(orbital_index, int) or orbital_index < 0:
        return _validation_error("orbital_index must be a non-negative integer.")

    if not isinstance(grid_size, int) or grid_size < 40 or grid_size > 120:
        return _validation_error("grid_size must be an integer between 40 and 120.")

    if not isinstance(isovalue_pos, (int, float)) or isovalue_pos < 0.001 or isovalue_pos > 0.1:
        return _validation_error("isovalue_pos must be a number between 0.001 and 0.1.")

    if not isinstance(isovalue_neg, (int, float)) or isovalue_neg < -0.1 or isovalue_neg > -0.001:
        return _validation_error("isovalue_neg must be a number between -0.1 and -0.001.")

    try:
        logger.debug(f"Generating CUBE file for orbital {orbital_index} in calculation {calculation_id}")
        quantum_service = get_quantum_service()
        cube_data = quantum_service.generate_orbital_cube(
            calculation_id,
            orbital_index,
            grid_size=grid_size,
            isovalue_pos=isovalue_pos,
            isovalue_neg=isovalue_neg
        )

        logger.info(f"Successfully generated CUBE file for orbital {orbital_index}")
        return json.dumps({
            'success': True,
            'data': cube_data
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"generate_orbital_cube(id={calculation_id}, orbital={orbital_index})")
    except Exception as e:
        logger.error(f"Unexpected error in generate_orbital_cube: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def list_cube_files(calculation_id: str) -> str:
    """
    List all generated CUBE files for a calculation.

    This function retrieves a list of all CUBE files that have been previously generated
    for molecular orbitals in a specific calculation. Useful for managing cached files
    and understanding what visualizations are available.

    Args:
        calculation_id (str): Unique ID of the calculation.
                             Example: "calc_20240101_120000_abcd1234"

    Returns:
        str: JSON string containing list of CUBE files.
             On success, includes:
             - Array of CUBE file information (filename, path, orbital index, size, etc.)
             - Total number of files
             - Total size in KB
             On error, contains structured error message.
    """
    # Input validation
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    if len(calculation_id.strip()) == 0:
        return _validation_error("calculation_id cannot be empty or contain only whitespace.")

    try:
        logger.debug(f"Listing CUBE files for calculation: {calculation_id}")
        quantum_service = get_quantum_service()
        result = quantum_service.list_cube_files(calculation_id)

        num_files = result.get('total_files', 0)
        logger.info(f"Found {num_files} CUBE files for calculation {calculation_id}")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"list_cube_files(id={calculation_id})")
    except Exception as e:
        logger.error(f"Unexpected error in list_cube_files: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def generate_ir_spectrum(
    calculation_id: str,
    broadening_fwhm: float = 100.0,
    x_min: float = 400.0,
    x_max: float = 4000.0,
    show_peaks: bool = True
) -> str:
    """
    Generate theoretical IR (infrared) spectrum from vibrational frequency data.

    This function creates an IR spectrum plot with Lorentzian broadening and scale factor corrections
    based on the calculation method and basis set. The spectrum is useful for comparing with
    experimental IR spectroscopy data and identifying functional groups.

    Note: This function requires that the calculation included vibrational frequency analysis.

    Args:
        calculation_id (str): Unique ID of the calculation.
                             Example: "calc_20240101_120000_abcd1234"
        broadening_fwhm (float): Full width at half maximum for Lorentzian broadening in cm⁻¹.
                                Higher values create broader, smoother peaks.
                                Range: 0.1-1000.0, Default: 100.0
        x_min (float): Minimum wavenumber for spectrum range in cm⁻¹.
                      Range: 0.0+, Default: 400.0
        x_max (float): Maximum wavenumber for spectrum range in cm⁻¹.
                      Range: up to 10000.0, Default: 4000.0
        show_peaks (bool): Whether to mark individual peaks in the plot.
                          Default: True

    Returns:
        str: JSON string containing IR spectrum data.
             On success, includes:
             - spectrum: x_axis (wavenumbers), y_axis (intensities), peaks, metadata
             - generation_info: Parameters used for generation
             On error, contains structured error message.
    """
    # Input validation
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    if len(calculation_id.strip()) == 0:
        return _validation_error("calculation_id cannot be empty or contain only whitespace.")

    if not isinstance(broadening_fwhm, (int, float)) or broadening_fwhm < 0.1 or broadening_fwhm > 1000.0:
        return _validation_error("broadening_fwhm must be a number between 0.1 and 1000.0.")

    if not isinstance(x_min, (int, float)) or x_min < 0.0:
        return _validation_error("x_min must be a non-negative number.")

    if not isinstance(x_max, (int, float)) or x_max > 10000.0:
        return _validation_error("x_max must be a number up to 10000.0.")

    if x_min >= x_max:
        return _validation_error("x_min must be less than x_max.")

    if not isinstance(show_peaks, bool):
        return _validation_error("show_peaks must be a boolean value.")

    try:
        logger.debug(f"Generating IR spectrum for calculation: {calculation_id}")
        quantum_service = get_quantum_service()
        result = quantum_service.generate_ir_spectrum(
            calculation_id,
            broadening_fwhm=broadening_fwhm,
            x_min=x_min,
            x_max=x_max,
            show_peaks=show_peaks
        )

        num_peaks = len(result.get('spectrum', {}).get('peaks', []))
        logger.info(f"Successfully generated IR spectrum with {num_peaks} peaks for calculation {calculation_id}")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"generate_ir_spectrum(id={calculation_id})")
    except Exception as e:
        logger.error(f"Unexpected error in generate_ir_spectrum: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


__all__ = [
    'list_all_calculations',
    'find_calculations',
    'get_calculation_details',
    'get_supported_parameters',
    'get_app_settings',
    'get_system_resources',
    'get_molecular_orbitals',
    'generate_orbital_cube',
    'list_cube_files',
    'generate_ir_spectrum',
]
