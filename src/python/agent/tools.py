"""
AI Agent Tool Wrapper Module

This module provides Python wrapper functions that enable AI agents to 
utilize PySCF_front application API endpoints.

Each function is optimized for Gemini SDK's Function Calling feature and 
serves as an API wrapper that strictly adheres to OpenAPI specifications.
"""

import json
import logging
import os
from typing import Optional, Dict, Any

import requests
from flask import current_app

# Logger configuration
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PORT = 5000
API_TIMEOUT = 30  # seconds
API_HOST = "127.0.0.1"


def get_api_base_url() -> str:
    """
    Get the API base URL with dynamic port detection.

    The port number is retrieved in the following priority order:
    1. Flask app.config['SERVER_PORT'] (if within Flask application context)
    2. Environment variable PYSCF_SERVER_PORT (fallback for backward compatibility)
    3. DEFAULT_PORT constant (5000)

    This ensures a single source of truth for the server port across the application.

    Returns:
        str: The base URL for API requests (e.g., "http://127.0.0.1:5000")
    """
    server_port = DEFAULT_PORT  # Default fallback

    try:
        # Priority 1: Try to get port from Flask app.config (single source of truth)
        # Accessing current_app.config will raise RuntimeError if outside of Flask context
        configured_port = current_app.config.get('SERVER_PORT')
        if configured_port is not None:
            server_port = configured_port
            logger.debug(f"Using port from app.config: {server_port}")
        else:
            # Fallback to environment variable if not set in app.config
            env_port = os.getenv('PYSCF_SERVER_PORT')
            if env_port:
                server_port = int(env_port)
                logger.debug(f"Using port from environment variable: {server_port}")
            else:
                logger.debug(f"Using default port: {server_port}")
    except RuntimeError:
        # Outside of Flask application context - fallback to environment variable
        env_port = os.getenv('PYSCF_SERVER_PORT')
        if env_port:
            server_port = int(env_port)
            logger.debug(f"Outside Flask context, using port from environment variable: {server_port}")
        else:
            logger.debug(f"Outside Flask context, using default port: {server_port}")

    return f"http://{API_HOST}:{server_port}"


def _handle_request_error(error: Exception, context: str = "") -> str:
    """
    Handle request errors and return appropriate error messages.
    
    Args:
        error: The exception that occurred
        context: Additional context about where the error occurred
        
    Returns:
        str: Formatted error message
    """
    context_prefix = f"{context}: " if context else ""
    
    if isinstance(error, requests.exceptions.HTTPError):
        status_code = error.response.status_code
        if status_code == 404:
            return f"Error: {context_prefix}Resource not found. Please verify the request parameters."
        elif status_code == 400:
            try:
                error_detail = error.response.json()
                return f"Error: {context_prefix}Invalid request - {error_detail.get('message', 'Bad request')}"
            except:
                return f"Error: {context_prefix}Invalid request parameters. Please check your input values."
        elif status_code >= 500:
            return f"Error: {context_prefix}Internal server error occurred. Please try again after a moment."
        else:
            return f"Error: {context_prefix}HTTP error occurred (status code {status_code})."
    
    elif isinstance(error, requests.exceptions.ConnectionError):
        logger.error(f"Connection error in {context}: {error}")
        return "Error: Could not connect to API server. Please check if the server is running."
    
    elif isinstance(error, requests.exceptions.Timeout):
        logger.error(f"Timeout error in {context}: {error}")
        return "Error: API request timed out. Please try again after a moment."
    
    elif isinstance(error, requests.exceptions.RequestException):
        logger.error(f"Request error in {context}: {error}")
        return f"Error: Network error occurred: {str(error)}"
    
    elif isinstance(error, json.JSONDecodeError):
        logger.error(f"JSON decode error in {context}: {error}")
        return "Error: Failed to parse API response."
    
    else:
        logger.error(f"Unexpected error in {context}: {error}")
        return f"Error: An unexpected error occurred: {str(error)}"


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
        logger.debug("Fetching all calculations from API")
        response = requests.get(f"{get_api_base_url()}/api/quantum/calculations", timeout=API_TIMEOUT)
        response.raise_for_status()
        return json.dumps(response.json(), ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "list_all_calculations")


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
        return "Error: calculation_id is a required string parameter."
    
    try:
        logger.debug(f"Fetching calculation details for ID: {calculation_id}")
        response = requests.get(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}", 
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        return json.dumps(response.json(), ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"get_calculation_details(id={calculation_id})")


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
        logger.debug("Fetching supported parameters from API")
        response = requests.get(f"{get_api_base_url()}/api/quantum/supported-parameters", timeout=API_TIMEOUT)
        response.raise_for_status()
        return json.dumps(response.json(), ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "get_supported_parameters")


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
                  Example: "3
Water molecule
O 0.0 0.0 0.0
H 0.7570 0.5860 0.0
H -0.7570 0.5860 0.0"
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
        return "Error: xyz parameter is required and must be a non-empty string containing molecular structure data."
    
    if not isinstance(calculation_method, str) or calculation_method not in ['DFT', 'HF', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF']:
        return "Error: calculation_method must be one of: 'DFT', 'HF', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF'"
    
    if not isinstance(charges, int) or charges < -10 or charges > 10:
        return "Error: charges must be an integer between -10 and 10."
    
    if not isinstance(spin, int) or spin < 0 or spin > 10:
        return "Error: spin must be an integer between 0 and 10 (number of unpaired electrons)."
    
    if not isinstance(name, str) or len(name.strip()) == 0 or len(name) > 100:
        return "Error: name must be a non-empty string with maximum 100 characters."
    
    # Validate solvent parameters
    if not isinstance(solvent_method, str) or solvent_method not in ['none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo']:
        return "Error: solvent_method must be one of: 'none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo'"
    
    # Validate resource parameters
    if cpu_cores is not None and (not isinstance(cpu_cores, int) or cpu_cores < 1 or cpu_cores > 32):
        return "Error: cpu_cores must be an integer between 1 and 32, or None for auto-allocation."
    
    if memory_mb is not None and (not isinstance(memory_mb, int) or memory_mb < 512 or memory_mb > 32768):
        return "Error: memory_mb must be an integer between 512 and 32768, or None for auto-allocation."
    
    # Validate TDDFT parameters
    if calculation_method == 'TDDFT':
        if not isinstance(tddft_nstates, int) or tddft_nstates < 1 or tddft_nstates > 50:
            return "Error: tddft_nstates must be an integer between 1 and 50 for TDDFT calculations."
        
        if not isinstance(tddft_method, str) or tddft_method not in ['TDDFT', 'TDA']:
            return "Error: tddft_method must be either 'TDDFT' or 'TDA'."
    
    # Validate CASCI/CASSCF parameters
    if calculation_method in ['CASCI', 'CASSCF']:
        if not isinstance(ncas, int) or ncas < 1 or ncas > 50:
            return "Error: ncas must be an integer between 1 and 50 for CASCI/CASSCF calculations."
        
        if not isinstance(nelecas, int) or nelecas < 1 or nelecas > 100:
            return "Error: nelecas must be an integer between 1 and 100 for CASCI/CASSCF calculations."
        
        if nelecas > 2 * ncas:
            return f"Error: nelecas ({nelecas}) cannot exceed 2 * ncas ({2 * ncas})."
        
        if calculation_method == 'CASSCF':
            if not isinstance(max_cycle_macro, int) or max_cycle_macro < 1 or max_cycle_macro > 200:
                return "Error: max_cycle_macro must be an integer between 1 and 200 for CASSCF calculations."
            
            if not isinstance(conv_tol, (int, float)) or conv_tol < 1e-12 or conv_tol > 1e-3:
                return "Error: conv_tol must be a number between 1e-12 and 1e-3 for CASSCF calculations."
            
            if not isinstance(conv_tol_grad, (int, float)) or conv_tol_grad < 1e-8 or conv_tol_grad > 1e-2:
                return "Error: conv_tol_grad must be a number between 1e-8 and 1e-2 for CASSCF calculations."
        
        if not isinstance(max_cycle_micro, int) or max_cycle_micro < 1 or max_cycle_micro > 20:
            return "Error: max_cycle_micro must be an integer between 1 and 20 for CASCI/CASSCF calculations."
    
    # Prepare request data
    request_data = {
        "xyz": xyz.strip(),
        "calculation_method": calculation_method,
        "basis_function": basis_function,
        "charges": charges,
        "spin": spin,
        "solvent_method": solvent_method,
        "solvent": solvent,
        "name": name.strip(),
        "tddft_nstates": tddft_nstates,
        "tddft_method": tddft_method,
        "tddft_analyze_nto": tddft_analyze_nto,
        "ncas": ncas,
        "nelecas": nelecas,
        "max_cycle_macro": max_cycle_macro,
        "max_cycle_micro": max_cycle_micro,
        "natorb": natorb,
        "conv_tol": conv_tol,
        "conv_tol_grad": conv_tol_grad,
        "optimize_geometry": optimize_geometry
    }
    
    # Add exchange_correlation only for methods that use it (not HF)
    if calculation_method != "HF" and exchange_correlation:
        request_data["exchange_correlation"] = exchange_correlation
    
    # Add resource parameters if specified
    if cpu_cores is not None:
        request_data["cpu_cores"] = cpu_cores
    
    if memory_mb is not None:
        request_data["memory_mb"] = memory_mb
    
    try:
        logger.debug(f"Starting quantum calculation: {calculation_method}/{basis_function} for '{name}'")
        response = requests.post(
            f"{get_api_base_url()}/api/quantum/calculate",
            json=request_data,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info(f"Successfully started calculation with ID: {result.get('id', 'unknown')}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "start_quantum_calculation")


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
        return "Error: compound_name parameter is required and must be a non-empty string."
    
    if not isinstance(search_type, str) or search_type not in ['name', 'cid', 'formula']:
        return "Error: search_type must be one of: 'name', 'cid', 'formula'"
    
    if len(compound_name.strip()) == 0:
        return "Error: compound_name cannot be empty or contain only whitespace."
    
    # Prepare request data
    request_data = {
        "query": compound_name.strip(),
        "searchType": search_type
    }
    
    try:
        logger.debug(f"Searching PubChem for compound: '{compound_name}' (type: {search_type})")
        response = requests.post(
            f"{get_api_base_url()}/api/pubchem/search",
            json=request_data,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info(f"Successfully found compound: {result.get('name', compound_name)}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"search_pubchem_by_name(compound={compound_name})")


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
        return "Error: calculation_id is a required string parameter."

    if len(calculation_id.strip()) == 0:
        return "Error: calculation_id cannot be empty or contain only whitespace."

    try:
        # First, get calculation details to provide context for the confirmation
        logger.debug(f"Fetching calculation details for deletion request: {calculation_id}")
        response = requests.get(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}",
            timeout=API_TIMEOUT
        )
        response.raise_for_status()

        calc_data = response.json()
        calculation_name = calc_data.get('data', {}).get('calculation', {}).get('name', 'Unknown Calculation')

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

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            return f"Error: Calculation with ID '{calculation_id}' not found. Please verify the calculation ID."
        else:
            return _handle_request_error(e, f"delete_calculation(id={calculation_id})")
    except Exception as e:
        return _handle_request_error(e, f"delete_calculation(id={calculation_id})")


def _execute_confirmed_deletion(calculation_id: str) -> Dict[str, Any]:
    """
    Internal function to execute a confirmed calculation deletion.

    This function should ONLY be called after user confirmation has been obtained
    through the Human-in-the-Loop confirmation mechanism. It actually performs
    the deletion by calling the DELETE API endpoint.

    Args:
        calculation_id (str): Unique ID of the calculation to delete.

    Returns:
        Dict[str, Any]: Result of the deletion operation.
                       - success (bool): Whether deletion succeeded
                       - message (str): Human-readable result message
                       - calculation_id (str): ID of the deleted calculation

    Raises:
        requests.exceptions.HTTPError: If the API request fails
        requests.exceptions.RequestException: If a network error occurs
    """
    if not calculation_id or not isinstance(calculation_id, str):
        return {
            "success": False,
            "message": "Invalid calculation_id parameter.",
            "calculation_id": None
        }

    try:
        logger.info(f"Executing confirmed deletion for calculation: {calculation_id}")
        response = requests.delete(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}",
            timeout=API_TIMEOUT
        )
        response.raise_for_status()

        result = response.json()
        logger.info(f"Successfully deleted calculation: {calculation_id}")

        return {
            "success": True,
            "message": result.get('data', {}).get('message', 'Calculation deleted successfully'),
            "calculation_id": calculation_id
        }

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            error_msg = f"Calculation with ID '{calculation_id}' not found."
            logger.warning(error_msg)
            return {
                "success": False,
                "message": error_msg,
                "calculation_id": calculation_id
            }
        else:
            error_msg = f"HTTP error occurred while deleting calculation: {str(e)}"
            logger.error(error_msg)
            return {
                "success": False,
                "message": error_msg,
                "calculation_id": calculation_id
            }
    except requests.exceptions.RequestException as e:
        error_msg = f"Network error occurred while deleting calculation: {str(e)}"
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
        return "Error: smiles parameter is required and must be a non-empty string."
    
    if len(smiles.strip()) == 0:
        return "Error: smiles string cannot be empty or contain only whitespace."
    
    # SMILES strings typically don't exceed 500 characters for reasonable molecules
    if len(smiles) > 500:
        return "Error: smiles string is too long. Maximum length is 500 characters."
    
    # Prepare request data
    request_data = {
        "smiles": smiles.strip()
    }
    
    try:
        logger.debug(f"Converting SMILES to XYZ: '{smiles}'")
        response = requests.post(
            f"{get_api_base_url()}/api/smiles/convert",
            json=request_data,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info(f"Successfully converted SMILES to XYZ structure")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"convert_smiles_to_xyz(smiles={smiles})")


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
        return "Error: xyz parameter is required and must be a non-empty string."
    
    if len(xyz.strip()) == 0:
        return "Error: xyz string cannot be empty or contain only whitespace."
    
    # Prepare request data
    request_data = {
        "xyz": xyz.strip()
    }
    
    try:
        logger.debug("Validating XYZ format")
        response = requests.post(
            f"{get_api_base_url()}/api/pubchem/validate",
            json=request_data,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        is_valid = result.get('data', {}).get('valid', False)
        logger.info(f"XYZ validation result: {'valid' if is_valid else 'invalid'}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "validate_xyz_format")



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
        return "Error: calculation_id is a required string parameter."
    
    if len(calculation_id.strip()) == 0:
        return "Error: calculation_id cannot be empty or contain only whitespace."
    
    try:
        logger.debug(f"Fetching molecular orbital information for calculation: {calculation_id}")
        response = requests.get(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}/orbitals",
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        num_orbitals = result.get('data', {}).get('total_orbitals', 0)
        logger.info(f"Successfully retrieved {num_orbitals} molecular orbitals for calculation {calculation_id}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"get_molecular_orbitals(id={calculation_id})")


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
        return "Error: calculation_id is a required string parameter."
    
    if len(calculation_id.strip()) == 0:
        return "Error: calculation_id cannot be empty or contain only whitespace."
    
    if not isinstance(orbital_index, int) or orbital_index < 0:
        return "Error: orbital_index must be a non-negative integer."
    
    if not isinstance(grid_size, int) or grid_size < 40 or grid_size > 120:
        return "Error: grid_size must be an integer between 40 and 120."
    
    if not isinstance(isovalue_pos, (int, float)) or isovalue_pos < 0.001 or isovalue_pos > 0.1:
        return "Error: isovalue_pos must be a number between 0.001 and 0.1."
    
    if not isinstance(isovalue_neg, (int, float)) or isovalue_neg < -0.1 or isovalue_neg > -0.001:
        return "Error: isovalue_neg must be a number between -0.1 and -0.001."
    
    # Prepare query parameters
    params = {
        "gridSize": grid_size,
        "isovaluePos": isovalue_pos,
        "isovalueNeg": isovalue_neg
    }
    
    try:
        logger.debug(f"Generating CUBE file for orbital {orbital_index} in calculation {calculation_id}")
        response = requests.get(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}/orbitals/{orbital_index}/cube",
            params=params,
            timeout=API_TIMEOUT * 2  # CUBE generation can take longer
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info(f"Successfully generated CUBE file for orbital {orbital_index}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"generate_orbital_cube(id={calculation_id}, orbital={orbital_index})")


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
        return "Error: calculation_id is a required string parameter."
    
    if len(calculation_id.strip()) == 0:
        return "Error: calculation_id cannot be empty or contain only whitespace."
    
    try:
        logger.debug(f"Listing CUBE files for calculation: {calculation_id}")
        response = requests.get(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}/orbitals/cube-files",
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        num_files = result.get('data', {}).get('total_files', 0)
        logger.info(f"Found {num_files} CUBE files for calculation {calculation_id}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"list_cube_files(id={calculation_id})")


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
        return "Error: calculation_id is a required string parameter."
    
    if len(calculation_id.strip()) == 0:
        return "Error: calculation_id cannot be empty or contain only whitespace."
    
    if orbital_index is not None and (not isinstance(orbital_index, int) or orbital_index < 0):
        return "Error: orbital_index must be a non-negative integer or None."
    
    # Prepare query parameters
    params = {}
    if orbital_index is not None:
        params["orbitalIndex"] = orbital_index
    
    try:
        if orbital_index is not None:
            logger.debug(f"Deleting CUBE file for orbital {orbital_index} in calculation {calculation_id}")
        else:
            logger.debug(f"Deleting all CUBE files for calculation {calculation_id}")
        
        response = requests.delete(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}/orbitals/cube-files",
            params=params,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        num_deleted = result.get('data', {}).get('deleted_files', 0)
        logger.info(f"Successfully deleted {num_deleted} CUBE file(s)")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        context = f"orbital {orbital_index}" if orbital_index is not None else "all orbitals"
        return _handle_request_error(e, f"delete_cube_files(id={calculation_id}, {context})")


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
        return "Error: calculation_id is a required string parameter."
    
    if len(calculation_id.strip()) == 0:
        return "Error: calculation_id cannot be empty or contain only whitespace."
    
    if not isinstance(broadening_fwhm, (int, float)) or broadening_fwhm < 0.1 or broadening_fwhm > 1000.0:
        return "Error: broadening_fwhm must be a number between 0.1 and 1000.0."
    
    if not isinstance(x_min, (int, float)) or x_min < 0.0:
        return "Error: x_min must be a non-negative number."
    
    if not isinstance(x_max, (int, float)) or x_max > 10000.0:
        return "Error: x_max must be a number up to 10000.0."
    
    if x_min >= x_max:
        return "Error: x_min must be less than x_max."
    
    if not isinstance(show_peaks, bool):
        return "Error: show_peaks must be a boolean value."
    
    # Prepare query parameters
    params = {
        "broadeningFwhm": broadening_fwhm,
        "xMin": x_min,
        "xMax": x_max,
        "showPeaks": show_peaks
    }
    
    try:
        logger.debug(f"Generating IR spectrum for calculation: {calculation_id}")
        response = requests.get(
            f"{get_api_base_url()}/api/quantum/calculations/{calculation_id}/ir-spectrum",
            params=params,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        num_peaks = len(result.get('data', {}).get('spectrum', {}).get('peaks', []))
        logger.info(f"Successfully generated IR spectrum with {num_peaks} peaks for calculation {calculation_id}")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, f"generate_ir_spectrum(id={calculation_id})")


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
        response = requests.get(
            f"{get_api_base_url()}/api/settings",
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info("Successfully retrieved application settings")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "get_app_settings")


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
            return "Error: max_parallel_instances must be an integer between 1 and 32."
    
    if max_cpu_utilization_percent is not None:
        if not isinstance(max_cpu_utilization_percent, (int, float)) or max_cpu_utilization_percent < 10.0 or max_cpu_utilization_percent > 100.0:
            return "Error: max_cpu_utilization_percent must be a number between 10.0 and 100.0."
    
    if max_memory_utilization_percent is not None:
        if not isinstance(max_memory_utilization_percent, (int, float)) or max_memory_utilization_percent < 10.0 or max_memory_utilization_percent > 100.0:
            return "Error: max_memory_utilization_percent must be a number between 10.0 and 100.0."
    
    if gemini_api_key is not None and not isinstance(gemini_api_key, str):
        return "Error: gemini_api_key must be a string."
    
    # First, get current settings
    try:
        logger.debug("Fetching current settings for update")
        get_response = requests.get(
            f"{get_api_base_url()}/api/settings",
            timeout=API_TIMEOUT
        )
        get_response.raise_for_status()
        current_settings = get_response.json().get('data', {}).get('settings', {})
    except Exception as e:
        return _handle_request_error(e, "update_app_settings (fetching current settings)")
    
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
        response = requests.put(
            f"{get_api_base_url()}/api/settings",
            json=update_data,
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info("Successfully updated application settings")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "update_app_settings")


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
        response = requests.get(
            f"{get_api_base_url()}/api/system/resource-status",
            timeout=API_TIMEOUT
        )
        response.raise_for_status()
        
        result = response.json()
        logger.info("Successfully retrieved system resource status")
        return json.dumps(result, ensure_ascii=False, indent=2)
    except Exception as e:
        return _handle_request_error(e, "get_system_resources")



