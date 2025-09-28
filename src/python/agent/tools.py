"""
AI Agent Tool Wrapper Module

This module provides Python wrapper functions that enable AI agents to 
utilize PySCF_front application API endpoints.

Each function is optimized for Gemini SDK's Function Calling feature and 
serves as an API wrapper that strictly adheres to OpenAPI specifications.
"""

import json
import logging
from typing import Optional

import requests

# ロガー設定
logger = logging.getLogger(__name__)

# Local API server base URL
API_BASE_URL = "http://127.0.0.1:5000"


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
        response = requests.get(f"{API_BASE_URL}/api/quantum/calculations", timeout=30)
        response.raise_for_status()
        
        # レスポンスをJSON文字列として返す
        return json.dumps(response.json(), ensure_ascii=False, indent=2)
        
    except requests.exceptions.HTTPError as e:
        error_msg = f"HTTP error occurred during API call: {e.response.status_code}"
        if e.response.status_code == 404:
            error_msg = "Calculation list not found."
        elif e.response.status_code >= 500:
            error_msg = "Internal server error occurred. Please try again after a moment."
        logger.error(f"HTTP error in list_all_calculations: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.ConnectionError as e:
        error_msg = "Could not connect to API server. Please check if the server is running."
        logger.error(f"Connection error in list_all_calculations: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.Timeout as e:
        error_msg = "API request timed out. Please try again after a moment."
        logger.error(f"Timeout error in list_all_calculations: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.RequestException as e:
        error_msg = f"Network error occurred: {str(e)}"
        logger.error(f"Request error in list_all_calculations: {e}")
        return f"Error: {error_msg}"
        
    except json.JSONDecodeError as e:
        error_msg = "Failed to parse API response."
        logger.error(f"JSON decode error in list_all_calculations: {e}")
        return f"Error: {error_msg}"
        
    except Exception as e:
        error_msg = f"An unexpected error occurred: {str(e)}"
        logger.error(f"Unexpected error in list_all_calculations: {e}")
        return f"Error: {error_msg}"


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
            timeout=30
        )
        response.raise_for_status()
        
        # レスポンスをJSON文字列として返す
        return json.dumps(response.json(), ensure_ascii=False, indent=2)
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            error_msg = f"Calculation with ID '{calculation_id}' not found. Please verify the ID is correct."
        elif e.response.status_code >= 500:
            error_msg = "Internal server error occurred. Please try again after a moment."
        else:
            error_msg = f"HTTP error occurred during API call: {e.response.status_code}"
        logger.error(f"HTTP error in get_calculation_details: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.ConnectionError as e:
        error_msg = "Could not connect to API server. Please check if the server is running."
        logger.error(f"Connection error in get_calculation_details: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.Timeout as e:
        error_msg = "API request timed out. Please try again after a moment."
        logger.error(f"Timeout error in get_calculation_details: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.RequestException as e:
        error_msg = f"Network error occurred: {str(e)}"
        logger.error(f"Request error in get_calculation_details: {e}")
        return f"Error: {error_msg}"
        
    except json.JSONDecodeError as e:
        error_msg = "Failed to parse API response."
        logger.error(f"JSON decode error in get_calculation_details: {e}")
        return f"Error: {error_msg}"
        
    except Exception as e:
        error_msg = f"An unexpected error occurred: {str(e)}"
        logger.error(f"Unexpected error in get_calculation_details: {e}")
        return f"Error: {error_msg}"


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
        response = requests.get(f"{API_BASE_URL}/api/quantum/supported-parameters", timeout=30)
        response.raise_for_status()
        
        # レスポンスをJSON文字列として返す
        return json.dumps(response.json(), ensure_ascii=False, indent=2)
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code >= 500:
            error_msg = "Internal server error occurred. Please try again after a moment."
        else:
            error_msg = f"HTTP error occurred during API call: {e.response.status_code}"
        logger.error(f"HTTP error in get_supported_parameters: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.ConnectionError as e:
        error_msg = "Could not connect to API server. Please check if the server is running."
        logger.error(f"Connection error in get_supported_parameters: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.Timeout as e:
        error_msg = "API request timed out. Please try again after a moment."
        logger.error(f"Timeout error in get_supported_parameters: {e}")
        return f"Error: {error_msg}"
        
    except requests.exceptions.RequestException as e:
        error_msg = f"Network error occurred: {str(e)}"
        logger.error(f"Request error in get_supported_parameters: {e}")
        return f"Error: {error_msg}"
        
    except json.JSONDecodeError as e:
        error_msg = "Failed to parse API response."
        logger.error(f"JSON decode error in get_supported_parameters: {e}")
        return f"Error: {error_msg}"
        
    except Exception as e:
        error_msg = f"An unexpected error occurred: {str(e)}"
        logger.error(f"Unexpected error in get_supported_parameters: {e}")
        return f"Error: {error_msg}"