"""
PubChem API endpoints.
Handles molecular data retrieval from PubChem database and XYZ validation.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from services import get_pubchem_service
from generated_models import PubChemSearchRequest, XYZValidateRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create pubchem blueprint
pubchem_bp = Blueprint('pubchem', __name__)


@pubchem_bp.route('/api/pubchem/search', methods=['POST'])
@validate()
def search_pubchem(body: PubChemSearchRequest):
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    pubchem_service = get_pubchem_service()

    query = body.query
    search_type = body.searchType.value

    # Call service layer
    result = pubchem_service.search_compound(query, search_type)

    return jsonify({
        'success': True,
        'data': result
    })


@pubchem_bp.route('/api/pubchem/validate', methods=['POST'])
@validate()
def validate_xyz_endpoint(body: XYZValidateRequest):
    """Validate an XYZ format string."""
    pubchem_service = get_pubchem_service()

    xyz_string = body.xyz

    # Call service layer
    validation_result = pubchem_service.validate_xyz(xyz_string)

    return jsonify({'success': True, 'data': validation_result})
