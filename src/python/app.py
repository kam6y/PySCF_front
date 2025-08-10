"""Flask API server for PubChem and SMILES integration."""

import logging
import os
from flask import Flask, request, jsonify
from flask_cors import CORS

from pubchem.client import PubChemClient, PubChemError
from pubchem import parser as xyz_parser
from SMILES.smiles_converter import smiles_to_xyz, SMILESError
from quantum_calc import DFTCalculator, CalculationError, ConvergenceError, InputError
from quantum_calc.file_manager import CalculationFileManager

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


@app.route('/api/quantum/calculate', methods=['POST'])
def quantum_calculate():
    """Perform quantum chemistry calculation using PySCF."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'success': False, 'error': 'Invalid request: No JSON data provided.'}), 400
        
        # Extract calculation parameters
        xyz_data = data.get('xyz', '').strip()
        if not xyz_data:
            return jsonify({'success': False, 'error': 'XYZ molecular structure is required.'}), 400
        
        calc_method = data.get('calculation_method', 'DFT')
        basis = data.get('basis_function', '6-31G(d)')
        xc_functional = data.get('exchange_correlation', 'B3LYP')
        charge = int(data.get('charges', 0))
        spin_mult = int(data.get('spin_multiplicity', 1))
        
        # Convert spin multiplicity to spin (2S+1 -> S)
        spin = (spin_mult - 1) // 2
        
        logger.info(f"Starting {calc_method} calculation with basis={basis}, xc={xc_functional}")
        
        # Validate calculation method
        if calc_method != 'DFT':
            return jsonify({'success': False, 'error': f'Calculation method "{calc_method}" not yet supported. Only DFT is available.'}), 400
        
        # Extract molecule name for file organization
        molecule_name = data.get('molecule_name', '').strip() or None
        
        # Initialize calculator with persistent file storage
        calculator = DFTCalculator(keep_files=True, molecule_name=molecule_name)
        
        # Parse XYZ and setup calculation
        atoms = calculator.parse_xyz(xyz_data)
        calculator.setup_calculation(
            atoms,
            basis=basis,
            xc=xc_functional,
            charge=charge,
            spin=spin,
            max_cycle=150
        )
        
        # Run calculation
        results = calculator.run_calculation()
        
        # Keep files for analysis (don't clean up immediately)
        calculator.cleanup(keep_files=True)
        
        logger.info(f"Calculation completed successfully. HOMO: {results['homo_index']}, LUMO: {results['lumo_index']}")
        
        return jsonify({
            'success': True,
            'data': {
                'calculation_results': results,
                'calculation_parameters': {
                    'method': calc_method,
                    'basis': basis,
                    'xc_functional': xc_functional,
                    'charge': charge,
                    'spin_multiplicity': spin_mult
                }
            }
        })
        
    except InputError as e:
        logger.error(f"Input validation error: {e}")
        return jsonify({'success': False, 'error': str(e)}), 400
    except ConvergenceError as e:
        logger.error(f"Calculation convergence error: {e}")
        return jsonify({'success': False, 'error': f'Calculation failed to converge: {str(e)}'}), 422
    except CalculationError as e:
        logger.error(f"Calculation error: {e}")
        return jsonify({'success': False, 'error': f'Calculation error: {str(e)}'}), 500
    except Exception as e:
        logger.error(f"Unexpected error in quantum calculation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred during calculation.'}), 500


@app.route('/api/quantum/calculations', methods=['GET'])
def list_calculations():
    """List all available calculation directories."""
    try:
        file_manager = CalculationFileManager()
        calculations = file_manager.list_calculations()
        
        return jsonify({
            'success': True,
            'data': {
                'base_directory': file_manager.get_base_directory(),
                'calculations': calculations,
                'count': len(calculations)
            }
        })
        
    except Exception as e:
        logger.error(f"Error listing calculations: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to list calculations.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['GET'])
def get_calculation_details(calculation_id):
    """Get detailed information about a specific calculation."""
    try:
        file_manager = CalculationFileManager()
        
        # Find the calculation directory
        calculations = file_manager.list_calculations()
        target_calc = None
        for calc in calculations:
            if calc['name'] == calculation_id:
                target_calc = calc
                break
        
        if not target_calc:
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        # Read calculation info if available
        calc_info = file_manager.read_calculation_info(target_calc['path'])
        
        return jsonify({
            'success': True,
            'data': {
                'calculation': {
                    'id': calculation_id,
                    'name': calculation_id.replace('_\\d{8}_\\d{6}$', ''),  # Remove timestamp suffix
                    'path': target_calc['path'],
                    'date': target_calc['date'],
                    'has_checkpoint': target_calc['has_checkpoint'],
                    'info': calc_info
                },
                'files': {
                    'checkpoint_exists': target_calc['has_checkpoint'],
                    'info_file_exists': calc_info is not None,
                    'geometry_file_exists': file_manager.file_exists(target_calc['path'], 'optimized_geometry.xyz')
                }
            }
        })
        
    except Exception as e:
        logger.error(f"Error getting calculation details for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to get calculation details.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['PUT'])
def update_calculation(calculation_id):
    """Update calculation metadata (name, etc.)."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'success': False, 'error': 'Invalid request: No JSON data provided.'}), 400
        
        # For now, this is a placeholder since renaming directories is complex
        # We'll return success but note that actual renaming isn't implemented
        new_name = data.get('name', '')
        if not new_name:
            return jsonify({'success': False, 'error': 'New name is required.'}), 400
        
        logger.info(f"Update request for calculation {calculation_id} with new name: {new_name}")
        
        return jsonify({
            'success': True,
            'data': {
                'message': f'Calculation name update noted (actual filesystem rename not implemented)',
                'calculation_id': calculation_id,
                'new_name': new_name
            }
        })
        
    except Exception as e:
        logger.error(f"Error updating calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to update calculation.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['DELETE'])
def delete_calculation(calculation_id):
    """Delete a calculation and its files."""
    try:
        import shutil
        
        file_manager = CalculationFileManager()
        
        # Find the calculation directory
        calculations = file_manager.list_calculations()
        target_calc = None
        for calc in calculations:
            if calc['name'] == calculation_id:
                target_calc = calc
                break
        
        if not target_calc:
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        # Delete the calculation directory
        calc_path = target_calc['path']
        if os.path.exists(calc_path):
            shutil.rmtree(calc_path)
            logger.info(f"Deleted calculation directory: {calc_path}")
        
        return jsonify({
            'success': True,
            'data': {
                'message': f'Calculation "{calculation_id}" has been deleted successfully',
                'deleted_path': calc_path
            }
        })
        
    except Exception as e:
        logger.error(f"Error deleting calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to delete calculation.'}), 500


@app.route('/api/quantum/test-simple', methods=['POST'])
def test_simple_calculation():
    """Test calculation with a simple water molecule."""
    try:
        # Simple water molecule XYZ data
        water_xyz = """3
water molecule test
O  0.0000   0.0000   0.0000
H  0.7571   0.5861   0.0000  
H -0.7571   0.5861   0.0000"""
        
        logger.info("Starting simple water molecule test calculation")
        
        # Initialize calculator
        calculator = DFTCalculator(keep_files=True, molecule_name="water_test")
        
        # Parse XYZ and setup calculation
        atoms = calculator.parse_xyz(water_xyz)
        calculator.setup_calculation(
            atoms,
            basis='STO-3G',  # Smaller basis for faster calculation
            xc='B3LYP',
            charge=0,
            spin=0,
            max_cycle=50  # Fewer cycles for test
        )
        
        # Run calculation
        results = calculator.run_calculation()
        
        # Keep files for analysis
        calculator.cleanup(keep_files=True)
        
        logger.info(f"Test calculation completed. Files saved to: {results['working_directory']}")
        
        return jsonify({
            'success': True,
            'data': {
                'message': 'Simple test calculation completed successfully',
                'working_directory': results['working_directory'],
                'homo_index': results['homo_index'],
                'lumo_index': results['lumo_index'],
                'scf_energy': results['scf_energy'],
                'converged': results['converged']
            }
        })
        
    except Exception as e:
        logger.error(f"Test calculation failed: {e}", exc_info=True)
        return jsonify({'success': False, 'error': f'Test calculation failed: {str(e)}'}), 500


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