"""Flask API server for PubChem integration."""

import logging
import sys
import os
from flask import Flask, request, jsonify
from flask_cors import CORS
from typing import Dict, Any

# Add the current directory to the path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pubchem.client import PubChemClient, PubChemError
from pubchem.parser import XYZParser

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)
CORS(app)  # Enable CORS for Electron integration

# Initialize PubChem client
pubchem_client = PubChemClient(timeout=30)
xyz_parser = XYZParser()


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'ok',
        'service': 'pyscf-pubchem-api',
        'version': '0.1.0'
    })


@app.route('/api/pubchem/search', methods=['POST'])
def search_pubchem():
    """Search PubChem for a compound and return XYZ coordinates.
    
    Expected JSON payload:
    {
        "query": "water" or "962",
        "search_type": "name" or "cid" (optional, defaults to "name")
    }
    
    Returns:
    {
        "success": true/false,
        "data": {
            "xyz": "xyz format string",
            "compound_info": {...}
        },
        "error": "error message" (if success is false)
    }
    """
    try:
        # Parse request data
        data = request.get_json(force=True, silent=True)
        if not data:
            return jsonify({
                'success': False,
                'error': 'No JSON data provided'
            }), 400
        
        query = data.get('query', '').strip()
        if not query:
            return jsonify({
                'success': False,
                'error': 'Query parameter is required'
            }), 400
        
        search_type = data.get('search_type', 'name').lower()
        if search_type not in ['name', 'cid', 'formula']:
            return jsonify({
                'success': False,
                'error': 'Invalid search_type. Must be "name", "cid", or "formula"'
            }), 400
        
        logger.info(f"Searching PubChem for: {query} (type: {search_type})")
        
        # Search for compound
        compound = pubchem_client.search_compound(query, search_type)
        if not compound:
            return jsonify({
                'success': False,
                'error': f'No compound found for query: {query}'
            }), 404
        
        logger.info(f"Found compound CID: {compound.cid}")
        
        # Get compound information
        compound_info = pubchem_client.get_compound_info(compound)
        
        # Try to get 3D structure
        atoms = pubchem_client.get_3d_structure(compound)
        if not atoms:
            return jsonify({
                'success': False,
                'error': f'No 3D structure available for {query} (CID: {compound.cid})'
            }), 404
        
        logger.info(f"Retrieved 3D structure with {len(atoms)} atoms")
        
        # Convert to XYZ format
        title = xyz_parser.format_compound_title(compound_info, query)
        xyz_string = xyz_parser.atoms_to_xyz(atoms, title)
        
        # Validate the generated XYZ
        validation = xyz_parser.validate_xyz(xyz_string)
        if not validation['valid']:
            logger.error(f"Generated invalid XYZ: {validation['error']}")
            return jsonify({
                'success': False,
                'error': 'Generated invalid XYZ format'
            }), 500
        
        logger.info(f"Successfully converted to XYZ format")
        
        return jsonify({
            'success': True,
            'data': {
                'xyz': xyz_string,
                'compound_info': compound_info,
                'atom_count': len(atoms)
            }
        })
        
    except PubChemError as e:
        logger.error(f"PubChem API error: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500
        
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        return jsonify({
            'success': False,
            'error': f'Data validation error: {e}'
        }), 400
        
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': 'Internal server error'
        }), 500


@app.route('/api/pubchem/validate', methods=['POST'])
def validate_xyz():
    """Validate XYZ format string.
    
    Expected JSON payload:
    {
        "xyz": "xyz format string"
    }
    
    Returns:
    {
        "success": true/false,
        "data": {
            "valid": true/false,
            "num_atoms": int (if valid),
            "error": "error message" (if invalid)
        }
    }
    """
    try:
        data = request.get_json(force=True, silent=True)
        if not data:
            return jsonify({
                'success': False,
                'error': 'No JSON data provided'
            }), 400
        
        xyz_string = data.get('xyz', '')
        if not xyz_string:
            return jsonify({
                'success': False,
                'error': 'XYZ string is required'
            }), 400
        
        # Validate XYZ format
        validation = xyz_parser.validate_xyz(xyz_string)
        
        return jsonify({
            'success': True,
            'data': validation
        })
        
    except Exception as e:
        logger.error(f"Error validating XYZ: {e}")
        return jsonify({
            'success': False,
            'error': 'Internal server error'
        }), 500


@app.errorhandler(404)
def not_found(error):
    """Handle 404 errors."""
    return jsonify({
        'success': False,
        'error': 'Endpoint not found'
    }), 404


@app.errorhandler(500)
def internal_error(error):
    """Handle 500 errors."""
    return jsonify({
        'success': False,
        'error': 'Internal server error'
    }), 500


if __name__ == '__main__':
    # Configuration
    host = '127.0.0.1'
    port = 5000
    debug = not os.environ.get('FLASK_ENV') == 'production'
    
    logger.info(f"Starting PySCF PubChem API server on {host}:{port}")
    logger.info(f"Debug mode: {debug}")
    
    try:
        app.run(host=host, port=port, debug=debug, threaded=True)
    except KeyboardInterrupt:
        logger.info("Server stopped by user")
    except Exception as e:
        logger.error(f"Server error: {e}")
        sys.exit(1)