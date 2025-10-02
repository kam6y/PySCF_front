"""
AI Agent Tool Wrapper Module

This module provides Python wrapper functions that enable AI agents to 
utilize PySCF_front application API endpoints.

Each function is optimized for Gemini SDK's Function Calling feature and 
serves as an API wrapper that strictly adheres to OpenAPI specifications.
"""

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

# Logger configuration
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PORT = 5000
API_TIMEOUT = 30  # seconds
API_HOST = "127.0.0.1"


def get_api_base_url() -> str:
    """Get the API base URL with dynamic port detection."""
    server_port = os.getenv('PYSCF_SERVER_PORT', str(DEFAULT_PORT))
    return f"http://{API_HOST}:{server_port}"


API_BASE_URL = get_api_base_url()


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
        response = requests.get(f"{API_BASE_URL}/api/quantum/calculations", timeout=API_TIMEOUT)
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
            f"{API_BASE_URL}/api/quantum/calculations/{calculation_id}", 
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
        response = requests.get(f"{API_BASE_URL}/api/quantum/supported-parameters", timeout=API_TIMEOUT)
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
    name: str = "AI Generated Calculation"
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
        exchange_correlation (str): Exchange-correlation functional for DFT. Examples: 'B3LYP', 'PBE0', 'M06-2X', 'CAM-B3LYP'
        name (str): Display name for the calculation. Used for identification in calculation list
    
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
    
    # Prepare request data
    request_data = {
        "xyz": xyz.strip(),
        "calculation_method": calculation_method,
        "basis_function": basis_function,
        "charges": charges,
        "spin": spin,
        "name": name.strip()
    }
    
    # Add exchange_correlation only for methods that use it (not HF)
    if calculation_method != "HF" and exchange_correlation:
        request_data["exchange_correlation"] = exchange_correlation
    
    try:
        logger.debug(f"Starting quantum calculation: {calculation_method}/{basis_function} for '{name}'")
        response = requests.post(
            f"{API_BASE_URL}/api/quantum/calculate",
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
            f"{API_BASE_URL}/api/pubchem/search",
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
            f"{API_BASE_URL}/api/quantum/calculations/{calculation_id}",
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
            f"{API_BASE_URL}/api/quantum/calculations/{calculation_id}",
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