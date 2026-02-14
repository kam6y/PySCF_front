"""
SMILES conversion API endpoints.
Handles conversion of SMILES strings to XYZ format.
"""

import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from services import get_smiles_service
from generated_models import SMILESConvertRequest

# Set up logging
logger = logging.getLogger(__name__)

# Create smiles blueprint
smiles_bp = Blueprint('smiles', __name__)


@smiles_bp.route('/api/smiles/convert', methods=['POST'])
@validate()
def convert_smiles(body: SMILESConvertRequest):
    """Converts a SMILES string to XYZ format."""
    smiles_service = get_smiles_service()

    smiles = body.smiles

    # Call service layer
    result = smiles_service.convert_smiles(smiles)

    return jsonify({'success': True, 'data': result})
