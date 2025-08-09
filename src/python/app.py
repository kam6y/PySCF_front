"""Flask API server for PubChem and SMILES integration."""

import logging
import os
from flask import Flask, request, jsonify
from flask_cors import CORS

from pubchem.client import PubChemClient, PubChemError
from pubchem import parser as xyz_parser
from SMILES.smiles_converter import smiles_to_xyz, SMILESError

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize Flask app and extensions
app = Flask(__name__)
CORS(app)  # Enable CORS for cross-origin requests

# Initialize PubChem client
pubchem_client = PubChemClient(timeout=30)


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'ok',
        'service': 'pyscf-native-api',
        'version': '0.3.0'
    })


@app.route('/api/pubchem/search', methods=['POST'])
def search_pubchem():
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'success': False, 'error': 'Invalid request: No JSON data provided or incorrect Content-Type.'}), 400
        
        query = data.get('query', '').strip()
        if not query:
            return jsonify({'success': False, 'error': 'Query parameter is required.'}), 400
        
        search_type = data.get('search_type', 'name').lower()
        if search_type not in ['name', 'cid', 'formula']:
            return jsonify({'success': False, 'error': 'Invalid search_type. Must be "name", "cid", or "formula".'}), 400
        
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
            
    except PubChemError as e:
        logger.error(f"A PubChem API error occurred: {e}")
        return jsonify({'success': False, 'error': str(e)}), 500
        
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/smiles/convert', methods=['POST'])
def convert_smiles():
    """Converts a SMILES string to XYZ format."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'success': False, 'error': 'Invalid request: No JSON data provided.'}), 400

        smiles = data.get('smiles', '').strip()
        if not smiles:
            return jsonify({'success': False, 'error': 'SMILES string is required.'}), 400

        logger.info(f"Converting SMILES: {smiles}")

        xyz_string = smiles_to_xyz(smiles, title=f"Molecule from SMILES: {smiles}")

        return jsonify({'success': True, 'data': {'xyz': xyz_string}})

    except SMILESError as e:
        logger.error(f"SMILES conversion failed: {e}")
        return jsonify({'success': False, 'error': str(e)}), 400
    except Exception as e:
        logger.error(f"An unexpected error occurred during SMILES conversion: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/pubchem/validate', methods=['POST'])
def validate_xyz_endpoint():
    """Validate an XYZ format string."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'success': False, 'error': 'Invalid request: No JSON data provided or incorrect Content-Type.'}), 400
        
        xyz_string = data.get('xyz', '')
        if not xyz_string:
            return jsonify({'success': False, 'error': 'XYZ string is required.'}), 400
        
        validation_result = xyz_parser.validate_xyz(xyz_string)
        
        return jsonify({'success': True, 'data': validation_result})
        
    except Exception as e:
        logger.error(f"Error during XYZ validation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.errorhandler(404)
def not_found(error):
    return jsonify({'success': False, 'error': 'Endpoint not found.'}), 404

@app.errorhandler(405)
def method_not_allowed(error):
    return jsonify({'success': False, 'error': 'Method not allowed for this endpoint.'}), 405

if __name__ == '__main__':
    host = os.environ.get('FLASK_RUN_HOST', '127.0.0.1')
    port = int(os.environ.get('FLASK_RUN_PORT', 5000))
    debug = os.environ.get('FLASK_ENV') == 'development'
    
    logger.info(f"Starting API server on http://{host}:{port} (Debug: {debug})")
    app.run(host=host, port=port, debug=debug)