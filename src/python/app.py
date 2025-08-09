"""Flask API server for PubChem integration."""

import logging
import os
import sys
from flask import Flask, request, jsonify
from flask_cors import CORS

# The following sys.path modification is removed as we assume the package
# is installed in editable mode (`uv pip install -e .`)
# sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pubchem.client import PubChemClient, PubChemError, CompoundData
from pubchem import parser as xyz_parser # Import the module directly

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
        'service': 'pyscf-pubchem-api',
        'version': '0.2.0' # Version updated
    })


@app.route('/api/pubchem/search', methods=['POST'])
def search_pubchem():
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    try:
        # Using request.get_json() without force=True is safer
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
        
        # Generate XYZ string
        title = xyz_parser.format_compound_title(compound_data, query)
        xyz_string = xyz_parser.atoms_to_xyz(compound_data.atoms, title)
        
        # The parser now includes robust validation, but an extra check here is fine.
        validation = xyz_parser.validate_xyz(xyz_string)
        if not validation['valid']:
            logger.error(f"Internal error: Generated invalid XYZ for CID {compound_data.cid}. Reason: {validation['error']}")
            return jsonify({'success': False, 'error': 'Internal server error: Failed to generate valid XYZ format.'}), 500
        
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
    # For production, use a proper WSGI server like Gunicorn or uWSGI
    # Configuration should be loaded from environment variables
    host = os.environ.get('FLASK_RUN_HOST', '127.0.0.1')
    port = int(os.environ.get('FLASK_RUN_PORT', 5000))
    debug = os.environ.get('FLASK_ENV') == 'development'
    
    logger.info(f"Starting API server on http://{host}:{port} (Debug: {debug})")
    app.run(host=host, port=port, debug=debug)