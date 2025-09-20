"""
PubChem API endpoints.
Handles molecular data retrieval from PubChem database and XYZ validation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate
from pydantic import ValidationError

from pubchem.client import PubChemClient, PubChemError, PubChemNotFoundError
from pubchem import parser as xyz_parser
from generated_models import PubChemSearchRequest, XYZValidateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create pubchem blueprint
pubchem_bp = Blueprint('pubchem', __name__)

# Initialize PubChem client
pubchem_client = PubChemClient(timeout=30)  # Default timeout, can be made configurable later


@pubchem_bp.route('/api/pubchem/search', methods=['POST'])
@validate()
def search_pubchem(body: PubChemSearchRequest):
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    try:
        query = body.query
        search_type = body.searchType.value
        
        logger.info(f"Searching PubChem for '{query}' (type: {search_type})")
        
        compound_data = pubchem_client.search_compound(query, search_type)
        if not compound_data or not compound_data.atoms:
            return jsonify({'success': False, 'error': f'No compound with a 3D structure found for query: {query}'}), 404
        
        logger.info(f"Found CID {compound_data.cid} with {len(compound_data.atoms)} atoms.")
        
        title = xyz_parser.format_compound_title(compound_data, query)
        xyz_string = xyz_parser.atoms_to_xyz(compound_data.atoms, title)
        
        logger.info(f"Successfully generated XYZ for CID {compound_data.cid}")
        
        return jsonify({
            'success': True,
            'data': {
                'xyz': xyz_string,
                'compound_info': {
                    'cid': compound_data.cid,
                    'iupac_name': compound_data.iupac_name,
                    'molecular_formula': compound_data.molecular_formula,
                    'molecular_weight': compound_data.molecular_weight,
                    'synonyms': compound_data.synonyms,
                },
                'atom_count': len(compound_data.atoms)
            }
        })
            
    except PubChemNotFoundError as e:
        logger.warning(f"PubChem search failed (Not Found): {e}")
        return jsonify({'success': False, 'error': str(e)}), 404
    except PubChemError as e:
        logger.error(f"A PubChem API error occurred: {e}", exc_info=True)
        return jsonify({'success': False, 'error': str(e)}), e.status_code or 500
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@pubchem_bp.route('/api/pubchem/validate', methods=['POST'])
@validate()
def validate_xyz_endpoint(body: XYZValidateRequest):
    """Validate an XYZ format string."""
    try:
        xyz_string = body.xyz
        
        # Check for empty or None input
        if not xyz_string or not xyz_string.strip():
            return jsonify({'success': False, 'error': 'XYZ string cannot be empty.'}), 400
        
        validation_result = xyz_parser.validate_xyz(xyz_string)
        
        return jsonify({'success': True, 'data': validation_result})
        
    except (TypeError, AttributeError) as e:
        logger.warning(f"Invalid input data for XYZ validation: {e}")
        return jsonify({'success': False, 'error': 'Invalid input format for XYZ validation.'}), 400
    except MemoryError as e:
        logger.error(f"Memory error during XYZ validation: {e}")
        return jsonify({'success': False, 'error': 'XYZ string too large to process.'}), 413
    except Exception as e:
        logger.error(f"Unexpected error during XYZ validation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500