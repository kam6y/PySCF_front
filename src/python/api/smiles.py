"""
SMILES conversion API endpoints.
Handles conversion of SMILES strings to XYZ format.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from services import get_smiles_service, ServiceError
from generated_models import SMILESConvertRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create smiles blueprint
smiles_bp = Blueprint('smiles', __name__)


@smiles_bp.route('/api/smiles/convert', methods=['POST'])
@validate()
def convert_smiles(body: SMILESConvertRequest):
    """Converts a SMILES string to XYZ format."""
    try:
        smiles_service = get_smiles_service()
        
        smiles = body.smiles
        
        # Call service layer
        result = smiles_service.convert_smiles(smiles)
        
        return jsonify({'success': True, 'data': result})

    except ServiceError as e:
        logger.error(f"Service error during SMILES conversion: {e}")
        return jsonify({'success': False, 'error': e.message}), e.status_code
    except Exception as e:
        logger.error(f"An unexpected error occurred during SMILES conversion: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500