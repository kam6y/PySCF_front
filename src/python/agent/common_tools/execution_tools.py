"""
Execution Tools Module

This module provides tools for executing actions that modify system state,
including calculation execution, deletion, molecular structure conversion,
and settings management.

These tools are designed for agents responsible for active operations and state changes.
"""

import json
import logging
from typing import Optional, Dict, Any

from services import (
    get_quantum_service, get_pubchem_service, get_smiles_service,
    get_settings_service, ServiceError
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


def start_quantum_calculation(
    xyz: str,
    calculation_method: str = "DFT",
    basis_function: str = "6-31G(d)",
    charges: int = 0,
    spin: int = 0,
    exchange_correlation: Optional[str] = "B3LYP",
    solvent_method: str = "none",
    solvent: str = "-",
    name: str = "AI Generated Calculation",
    cpu_cores: Optional[int] = None,
    memory_mb: Optional[int] = None,
    tddft_nstates: int = 10,
    tddft_method: str = "TDDFT",
    tddft_analyze_nto: bool = False,
    ncas: int = 4,
    nelecas: int = 4,
    max_cycle_macro: int = 50,
    max_cycle_micro: int = 3,
    natorb: bool = True,
    conv_tol: float = 0.000001,
    conv_tol_grad: float = 0.0001,
    optimize_geometry: bool = True
) -> str:
    """
    Start a new quantum chemistry calculation with the specified molecular structure and parameters.

    This function initiates an asynchronous quantum chemistry calculation using PySCF.
    The calculation will run in the background and can be monitored using get_calculation_details.

    Args:
        xyz (str): XYZ molecular structure data (atomic symbols and coordinates).
                  Example: "3\nWater molecule\nO 0.0 0.0 0.0\nH 0.7570 0.5860 0.0\nH -0.7570 0.5860 0.0"
        calculation_method (str): Quantum chemistry method. Options: 'DFT', 'HF', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF'
        basis_function (str): Basis set for calculation. Examples: 'STO-3G', '6-31G(d)', '6-31+G(d,p)', 'cc-pVDZ', 'aug-cc-pVTZ', 'def2-SVP'
        charges (int): Molecular charge (-10 to 10). Default: 0 (neutral molecule)
        spin (int): Number of unpaired electrons (2S). Default: 0 (closed shell singlet)
        exchange_correlation (str): Exchange-correlation functional for DFT/TDDFT. Examples: 'B3LYP', 'PBE0', 'M06-2X', 'CAM-B3LYP'. Ignored for HF method.
        solvent_method (str): Solvent effect method. Options: 'none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo'. Default: 'none'
        solvent (str): Solvent type or dielectric constant. Predefined: 'water', 'methanol', 'ethanol', 'acetone', etc., or numeric value >1.0. Default: '-'
        name (str): Display name for the calculation. Used for identification in calculation list
        cpu_cores (int): Number of CPU cores to use (1-32). If None, system will auto-allocate based on availability
        memory_mb (int): Memory in MB (512-32768). If None, system will auto-allocate based on availability
        tddft_nstates (int): Number of excited states to calculate (TDDFT only, 1-50). Default: 10
        tddft_method (str): TDDFT calculation method. Options: 'TDDFT', 'TDA' (Tamm-Dancoff approximation). Default: 'TDDFT'
        tddft_analyze_nto (bool): Perform Natural Transition Orbital analysis (TDDFT only). Default: False
        ncas (int): Number of active space orbitals (CASCI/CASSCF only, 1-50). Default: 4
        nelecas (int): Number of active space electrons (CASCI/CASSCF only, 1-100). Default: 4
        max_cycle_macro (int): Maximum CASSCF macro iterations (CASSCF only, 1-200). Default: 50
        max_cycle_micro (int): Maximum CI solver micro iterations (CASCI/CASSCF, 1-20). Default: 3
        natorb (bool): Transform to natural orbitals in active space (CASCI/CASSCF only). Default: True
        conv_tol (float): Energy convergence tolerance (CASSCF only, 1e-12 to 1e-3). Default: 0.000001
        conv_tol_grad (float): Gradient convergence tolerance (CASSCF only, 1e-8 to 1e-2). Default: 0.0001
        optimize_geometry (bool): Whether to perform geometry optimization before the main calculation. Default: True

    Returns:
        str: JSON string of the created calculation instance with ID and initial status.
             On error, contains structured error message.
    """
    # Input validation
    if not xyz or not isinstance(xyz, str):
        return _validation_error("xyz parameter is required and must be a non-empty string containing molecular structure data.")

    if not isinstance(calculation_method, str) or calculation_method not in ['DFT', 'HF', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF']:
        return _validation_error("calculation_method must be one of: 'DFT', 'HF', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF'")

    if not isinstance(charges, int) or charges < -10 or charges > 10:
        return _validation_error("charges must be an integer between -10 and 10.")

    if not isinstance(spin, int) or spin < 0 or spin > 10:
        return _validation_error("spin must be an integer between 0 and 10 (number of unpaired electrons).")

    if not isinstance(name, str) or len(name.strip()) == 0 or len(name) > 100:
        return _validation_error("name must be a non-empty string with maximum 100 characters.")

    # Validate solvent parameters
    if not isinstance(solvent_method, str) or solvent_method not in ['none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo']:
        return _validation_error("solvent_method must be one of: 'none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo'")

    # Validate resource parameters
    if cpu_cores is not None and (not isinstance(cpu_cores, int) or cpu_cores < 1 or cpu_cores > 32):
        return _validation_error("cpu_cores must be an integer between 1 and 32, or None for auto-allocation.")

    if memory_mb is not None and (not isinstance(memory_mb, int) or memory_mb < 512 or memory_mb > 32768):
        return _validation_error("memory_mb must be an integer between 512 and 32768, or None for auto-allocation.")

    # Validate TDDFT parameters
    if calculation_method == 'TDDFT':
        if not isinstance(tddft_nstates, int) or tddft_nstates < 1 or tddft_nstates > 50:
            return _validation_error("tddft_nstates must be an integer between 1 and 50 for TDDFT calculations.")

        if not isinstance(tddft_method, str) or tddft_method not in ['TDDFT', 'TDA']:
            return _validation_error("tddft_method must be either 'TDDFT' or 'TDA'.")

    # Validate CASCI/CASSCF parameters
    if calculation_method in ['CASCI', 'CASSCF']:
        if not isinstance(ncas, int) or ncas < 1 or ncas > 50:
            return _validation_error("ncas must be an integer between 1 and 50 for CASCI/CASSCF calculations.")

        if not isinstance(nelecas, int) or nelecas < 1 or nelecas > 100:
            return _validation_error("nelecas must be an integer between 1 and 100 for CASCI/CASSCF calculations.")

        if nelecas > 2 * ncas:
            return _validation_error(f"nelecas ({nelecas}) cannot exceed 2 * ncas ({2 * ncas}).")

        if calculation_method == 'CASSCF':
            if not isinstance(max_cycle_macro, int) or max_cycle_macro < 1 or max_cycle_macro > 200:
                return _validation_error("max_cycle_macro must be an integer between 1 and 200 for CASSCF calculations.")

            if not isinstance(conv_tol, (int, float)) or conv_tol < 1e-12 or conv_tol > 1e-3:
                return _validation_error("conv_tol must be a number between 1e-12 and 1e-3 for CASSCF calculations.")

            if not isinstance(conv_tol_grad, (int, float)) or conv_tol_grad < 1e-8 or conv_tol_grad > 1e-2:
                return _validation_error("conv_tol_grad must be a number between 1e-8 and 1e-2 for CASSCF calculations.")

        if not isinstance(max_cycle_micro, int) or max_cycle_micro < 1 or max_cycle_micro > 20:
            return _validation_error("max_cycle_micro must be an integer between 1 and 20 for CASCI/CASSCF calculations.")

    # Build parameters dictionary
    from datetime import datetime

    parameters = {
        'xyz': xyz.strip(),
        'calculation_method': calculation_method,
        'basis_function': basis_function,
        'charges': charges,
        'spin': spin,
        'solvent_method': solvent_method,
        'solvent': solvent,
        'name': name.strip(),
        'tddft_nstates': tddft_nstates,
        'tddft_method': tddft_method,
        'tddft_analyze_nto': tddft_analyze_nto,
        'ncas': ncas,
        'nelecas': nelecas,
        'max_cycle_macro': max_cycle_macro,
        'max_cycle_micro': max_cycle_micro,
        'natorb': natorb,
        'conv_tol': conv_tol,
        'conv_tol_grad': conv_tol_grad,
        'optimize_geometry': optimize_geometry,
        'created_at': datetime.now().isoformat()
    }

    # Add exchange_correlation only for methods that use it (not HF)
    if calculation_method != "HF" and exchange_correlation:
        parameters["exchange_correlation"] = exchange_correlation
    else:
        parameters["exchange_correlation"] = None

    # Add resource parameters if specified
    if cpu_cores is not None:
        parameters["cpu_cores"] = cpu_cores

    if memory_mb is not None:
        parameters["memory_mb"] = memory_mb

    try:
        logger.debug(f"Starting quantum calculation: {calculation_method}/{basis_function} for '{name}'")
        quantum_service = get_quantum_service()
        result = quantum_service.start_calculation(parameters)

        logger.info(f"Successfully started calculation with ID: {result.get('id', 'unknown')}")
        return json.dumps({
            'success': True,
            'data': {'calculation': result}
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "start_quantum_calculation")
    except Exception as e:
        logger.error(f"Unexpected error in start_quantum_calculation: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def delete_calculation(calculation_id: str) -> str:
    """
    Request deletion of a calculation and all its associated files.

    IMPORTANT: This is a DESTRUCTIVE operation that cannot be undone.
    This tool does NOT directly delete the calculation. Instead, it returns
    a confirmation request that must be approved by the user through the UI.

    The actual deletion will only occur after explicit user confirmation.

    Args:
        calculation_id (str): Unique ID of the calculation to delete.
                             Example: "calc_20240101_120000_abcd1234"

    Returns:
        str: JSON string containing confirmation request details.
             This should be presented to the user for approval before deletion.
    """
    # Input validation
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    if len(calculation_id.strip()) == 0:
        return _validation_error("calculation_id cannot be empty or contain only whitespace.")

    try:
        # First, get calculation details to provide context for the confirmation
        logger.debug(f"Fetching calculation details for deletion request: {calculation_id}")
        quantum_service = get_quantum_service()
        calc_data = quantum_service.get_calculation_details(calculation_id)

        calculation_name = calc_data.get('calculation', {}).get('name', 'Unknown Calculation')

        # Return a structured confirmation request instead of deleting directly
        confirmation_request = {
            "requires_confirmation": True,
            "action": "delete_calculation",
            "calculation_id": calculation_id.strip(),
            "calculation_name": calculation_name,
            "message": f"Are you sure you want to permanently delete the calculation '{calculation_name}' (ID: {calculation_id})? This action cannot be undone."
        }

        logger.info(f"Returning confirmation request for deletion of calculation: {calculation_id}")
        return json.dumps(confirmation_request, ensure_ascii=False, indent=2)

    except ServiceError as e:
        return _handle_service_error(e, f"delete_calculation(id={calculation_id})")
    except Exception as e:
        logger.error(f"Unexpected error in delete_calculation: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def _execute_confirmed_deletion(calculation_id: str) -> Dict[str, Any]:
    """
    Internal function to execute a confirmed calculation deletion.

    This function should ONLY be called after user confirmation has been obtained
    through the Human-in-the-Loop confirmation mechanism. It actually performs
    the deletion by calling the service layer directly.

    Args:
        calculation_id (str): Unique ID of the calculation to delete.

    Returns:
        Dict[str, Any]: Result of the deletion operation.
                       - success (bool): Whether deletion succeeded
                       - message (str): Human-readable result message
                       - calculation_id (str): ID of the deleted calculation

    Raises:
        ServiceError: If the deletion fails
    """
    if not calculation_id or not isinstance(calculation_id, str):
        return {
            "success": False,
            "message": "Invalid calculation_id parameter.",
            "calculation_id": None
        }

    try:
        logger.info(f"Executing confirmed deletion for calculation: {calculation_id}")
        quantum_service = get_quantum_service()
        result = quantum_service.delete_calculation(calculation_id)

        logger.info(f"Successfully deleted calculation: {calculation_id}")

        return {
            "success": True,
            "message": result.get('message', 'Calculation deleted successfully'),
            "calculation_id": calculation_id
        }

    except ServiceError as e:
        error_msg = f"Error deleting calculation: {e.message}"
        logger.error(error_msg)
        return {
            "success": False,
            "message": error_msg,
            "calculation_id": calculation_id
        }
    except Exception as e:
        error_msg = f"Unexpected error occurred while deleting calculation: {str(e)}"
        logger.error(error_msg, exc_info=True)
        return {
            "success": False,
            "message": error_msg,
            "calculation_id": calculation_id
        }


def delete_cube_files(calculation_id: str, orbital_index: Optional[int] = None) -> str:
    """
    Delete CUBE files for a calculation.

    This function removes generated CUBE files to free up disk space.
    You can delete a specific orbital's CUBE file or all CUBE files for the calculation.

    Args:
        calculation_id (str): Unique ID of the calculation.
                             Example: "calc_20240101_120000_abcd1234"
        orbital_index (int, optional): Specific orbital index to delete.
                                      If None, deletes all CUBE files for the calculation.
                                      Default: None (delete all)

    Returns:
        str: JSON string containing deletion result.
             On success, includes:
             - Number of files deleted
             - Confirmation message
             On error, contains structured error message.
    """
    # Input validation
    if not calculation_id or not isinstance(calculation_id, str):
        return _validation_error("calculation_id is a required string parameter.")

    if len(calculation_id.strip()) == 0:
        return _validation_error("calculation_id cannot be empty or contain only whitespace.")

    if orbital_index is not None and (not isinstance(orbital_index, int) or orbital_index < 0):
        return _validation_error("orbital_index must be a non-negative integer or None.")

    try:
        if orbital_index is not None:
            logger.debug(f"Deleting CUBE file for orbital {orbital_index} in calculation {calculation_id}")
        else:
            logger.debug(f"Deleting all CUBE files for calculation {calculation_id}")

        quantum_service = get_quantum_service()
        result = quantum_service.delete_cube_files(calculation_id, orbital_index)

        num_deleted = result.get('deleted_files', 0)
        logger.info(f"Successfully deleted {num_deleted} CUBE file(s)")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        context = f"orbital {orbital_index}" if orbital_index is not None else "all orbitals"
        return _handle_service_error(e, f"delete_cube_files(id={calculation_id}, {context})")
    except Exception as e:
        logger.error(f"Unexpected error in delete_cube_files: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def convert_smiles_to_xyz(smiles: str) -> str:
    """
    Convert a SMILES string to 3D XYZ molecular structure format.

    This function takes a SMILES (Simplified Molecular Input Line Entry System) string
    and converts it to a 3D molecular structure in XYZ format suitable for quantum chemistry calculations.
    The conversion includes 3D coordinate generation and geometry optimization.

    Args:
        smiles (str): SMILES string representing the molecular structure.
                     Examples: 'CCO' (ethanol), 'c1ccccc1' (benzene), 'CC(=O)O' (acetic acid)

    Returns:
        str: JSON string containing the converted XYZ structure.
             On success, includes xyz field with the molecular structure.
             On error, contains structured error message.
    """
    # Input validation
    if not smiles or not isinstance(smiles, str):
        return _validation_error("smiles parameter is required and must be a non-empty string.")

    if len(smiles.strip()) == 0:
        return _validation_error("smiles string cannot be empty or contain only whitespace.")

    # SMILES strings typically don't exceed 500 characters for reasonable molecules
    if len(smiles) > 500:
        return _validation_error("smiles string is too long. Maximum length is 500 characters.")

    try:
        logger.debug(f"Converting SMILES to XYZ: '{smiles}'")
        smiles_service = get_smiles_service()
        result = smiles_service.convert_smiles(smiles.strip())

        logger.info(f"Successfully converted SMILES to XYZ structure")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"convert_smiles_to_xyz(smiles={smiles})")
    except Exception as e:
        logger.error(f"Unexpected error in convert_smiles_to_xyz: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def validate_xyz_format(xyz: str) -> str:
    """
    Validate the format of an XYZ molecular structure string.

    This function checks whether a given string is properly formatted as XYZ molecular structure data.
    It validates the number of atoms, element symbols, and coordinate values.

    Args:
        xyz (str): XYZ format string to validate.
                  Expected format:
                  Line 1: Number of atoms
                  Line 2: Comment line (molecule name)
                  Lines 3+: Element X Y Z (coordinates in Angstroms)

    Returns:
        str: JSON string containing validation results.
             On success, includes valid=true and parsed structure information.
             On validation failure, includes valid=false and error message.
    """
    # Input validation
    if not xyz or not isinstance(xyz, str):
        return _validation_error("xyz parameter is required and must be a non-empty string.")

    if len(xyz.strip()) == 0:
        return _validation_error("xyz string cannot be empty or contain only whitespace.")

    try:
        logger.debug("Validating XYZ format")
        pubchem_service = get_pubchem_service()
        result = pubchem_service.validate_xyz(xyz.strip())

        is_valid = result.get('valid', False)
        logger.info(f"XYZ validation result: {'valid' if is_valid else 'invalid'}")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "validate_xyz_format")
    except Exception as e:
        logger.error(f"Unexpected error in validate_xyz_format: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def search_pubchem_by_name(
    compound_name: str,
    search_type: str = "name"
) -> str:
    """
    Search PubChem database for molecular compounds by name, CID, or formula.

    This function queries the PubChem database to find molecular structures
    and returns 3D coordinates in XYZ format that can be used for quantum chemistry calculations.

    Args:
        compound_name (str): Name, CID, or formula of the compound to search.
                           Examples: 'water', 'benzene', 'aspirin', 'C6H6', '241' (CID)
        search_type (str): Type of search to perform. Options: 'name', 'cid', 'formula'.
                          Default: 'name' (search by compound name)

    Returns:
        str: JSON string containing compound information and XYZ coordinates.
             On success, includes molecular structure data suitable for calculations.
             On error, contains structured error message.
    """
    # Input validation
    if not compound_name or not isinstance(compound_name, str):
        return _validation_error("compound_name parameter is required and must be a non-empty string.")

    if not isinstance(search_type, str) or search_type not in ['name', 'cid', 'formula']:
        return _validation_error("search_type must be one of: 'name', 'cid', 'formula'")

    if len(compound_name.strip()) == 0:
        return _validation_error("compound_name cannot be empty or contain only whitespace.")

    try:
        logger.debug(f"Searching PubChem for compound: '{compound_name}' (type: {search_type})")
        pubchem_service = get_pubchem_service()
        result = pubchem_service.search_compound(compound_name.strip(), search_type)

        logger.info(f"Successfully found compound")
        return json.dumps({
            'success': True,
            'data': result
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, f"search_pubchem_by_name(compound={compound_name})")
    except Exception as e:
        logger.error(f"Unexpected error in search_pubchem_by_name: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


def update_app_settings(
    max_parallel_instances: Optional[int] = None,
    max_cpu_utilization_percent: Optional[float] = None,
    max_memory_utilization_percent: Optional[float] = None,
    gemini_api_key: Optional[str] = None
) -> str:
    """
    Update application settings.

    This function modifies the application configuration settings. All parameters are optional;
    only provided parameters will be updated. This allows for partial updates without
    affecting other settings.

    Args:
        max_parallel_instances (int, optional): Maximum number of parallel calculation instances.
                                               Range: 1-32
        max_cpu_utilization_percent (float, optional): Maximum CPU utilization percentage.
                                                       Range: 10.0-100.0
        max_memory_utilization_percent (float, optional): Maximum memory utilization percentage.
                                                          Range: 10.0-100.0
        gemini_api_key (str, optional): Google Gemini API key for AI agent functionality.
                                       Pass empty string to remove the key.

    Returns:
        str: JSON string containing updated settings.
             On success, includes all application settings with updated values.
             On error, contains structured error message.
    """
    # Input validation
    if max_parallel_instances is not None:
        if not isinstance(max_parallel_instances, int) or max_parallel_instances < 1 or max_parallel_instances > 32:
            return _validation_error("max_parallel_instances must be an integer between 1 and 32.")

    if max_cpu_utilization_percent is not None:
        if not isinstance(max_cpu_utilization_percent, (int, float)) or max_cpu_utilization_percent < 10.0 or max_cpu_utilization_percent > 100.0:
            return _validation_error("max_cpu_utilization_percent must be a number between 10.0 and 100.0.")

    if max_memory_utilization_percent is not None:
        if not isinstance(max_memory_utilization_percent, (int, float)) or max_memory_utilization_percent < 10.0 or max_memory_utilization_percent > 100.0:
            return _validation_error("max_memory_utilization_percent must be a number between 10.0 and 100.0.")

    if gemini_api_key is not None and not isinstance(gemini_api_key, str):
        return _validation_error("gemini_api_key must be a string.")

    # First, get current settings
    try:
        logger.debug("Fetching current settings for update")
        settings_service = get_settings_service()
        current_settings = settings_service.get_settings()
    except ServiceError as e:
        return _handle_service_error(e, "update_app_settings (fetching current settings)")
    except Exception as e:
        logger.error(f"Unexpected error fetching current settings: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)

    # Prepare update data with current values as defaults
    update_data = {
        "max_parallel_instances": current_settings.get('max_parallel_instances', 4),
        "max_cpu_utilization_percent": current_settings.get('max_cpu_utilization_percent', 95.0),
        "max_memory_utilization_percent": current_settings.get('max_memory_utilization_percent', 95.0),
        "system_total_cores": current_settings.get('system_total_cores', 8),
        "system_total_memory_mb": current_settings.get('system_total_memory_mb', 16384),
        "gemini_api_key": current_settings.get('gemini_api_key')
    }

    # Update with provided values
    if max_parallel_instances is not None:
        update_data["max_parallel_instances"] = max_parallel_instances

    if max_cpu_utilization_percent is not None:
        update_data["max_cpu_utilization_percent"] = max_cpu_utilization_percent

    if max_memory_utilization_percent is not None:
        update_data["max_memory_utilization_percent"] = max_memory_utilization_percent

    if gemini_api_key is not None:
        update_data["gemini_api_key"] = gemini_api_key if gemini_api_key else None

    try:
        logger.debug("Updating application settings")
        settings_service = get_settings_service()
        updated_settings = settings_service.update_settings(update_data)

        logger.info("Successfully updated application settings")
        return json.dumps({
            'success': True,
            'data': {'settings': updated_settings}
        }, ensure_ascii=False, indent=2)
    except ServiceError as e:
        return _handle_service_error(e, "update_app_settings")
    except Exception as e:
        logger.error(f"Unexpected error in update_app_settings: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"An unexpected error occurred: {str(e)}"
        }, ensure_ascii=False, indent=2)


__all__ = [
    'start_quantum_calculation',
    'delete_calculation',
    '_execute_confirmed_deletion',
    'delete_cube_files',
    'convert_smiles_to_xyz',
    'validate_xyz_format',
    'search_pubchem_by_name',
    'update_app_settings',
]
