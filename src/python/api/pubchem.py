"""
PubChem API endpoints.
Handles molecular data retrieval from PubChem database and XYZ validation.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from services import get_pubchem_service, ServiceError
from generated_models import PubChemSearchRequest, XYZValidateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create pubchem blueprint
pubchem_bp = Blueprint('pubchem', __name__)


@pubchem_bp.route('/api/pubchem/search', methods=['POST'])
@validate()
def search_pubchem(body: PubChemSearchRequest):
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    try:
        pubchem_service = get_pubchem_service()
        
        query = body.query
        search_type = body.searchType.value
        
        # Call service layer
        result = pubchem_service.search_compound(query, search_type)
        
        return jsonify({
            'success': True,
            'data': result
        })
            
    except ServiceError as e:
        logger.error(f"Service error during PubChem search: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@pubchem_bp.route('/api/pubchem/validate', methods=['POST'])
@validate()
def validate_xyz_endpoint(body: XYZValidateRequest):
    """Validate an XYZ format string."""
    try:
        pubchem_service = get_pubchem_service()
        
        xyz_string = body.xyz
        
        # Call service layer
        validation_result = pubchem_service.validate_xyz(xyz_string)
        
        return jsonify({'success': True, 'data': validation_result})
        
    except ServiceError as e:
        logger.error(f"Service error during XYZ validation: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"Unexpected error during XYZ validation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500