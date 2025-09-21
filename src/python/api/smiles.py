"""
SMILES conversion API endpoints.
Handles conversion of SMILES strings to XYZ format.
"""


import logging
from flask import Blueprint, jsonify
from flask_pydantic import validate

from SMILES.smiles_converter import smiles_to_xyz, SMILESError
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
        smiles = body.smiles
        
        logger.info(f"Converting SMILES: {smiles}")

        xyz_string = smiles_to_xyz(smiles, title=f"Molecule from SMILES: {smiles}")

        return jsonify({'success': True, 'data': {'xyz': xyz_string}})

    except SMILESError as e:
        logger.error(f"SMILES conversion failed: {e}")
        return jsonify({'success': False, 'error': str(e)}), 400
    except Exception as e:
        logger.error(f"An unexpected error occurred during SMILES conversion: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500