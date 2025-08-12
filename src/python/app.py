# src/python/app.py

import logging
import os
import sys
import json  # jsonをインポート
import time  # timeをインポート
from datetime import datetime
from flask import Flask, request, jsonify
from flask_cors import CORS
from flask_sock import Sock  # flask-sockをインポート
from flask_pydantic import validate
from pydantic import ValidationError
import threading
import socket
import shutil
from typing import Dict
from werkzeug.serving import make_server

from pubchem.client import PubChemClient, PubChemError, PubChemNotFoundError
from pubchem import parser as xyz_parser
from SMILES.smiles_converter import smiles_to_xyz, SMILESError
from quantum_calc import DFTCalculator, CalculationError, ConvergenceError, InputError
from quantum_calc.file_manager import CalculationFileManager
from generated_models import (
    PubChemSearchRequest, SMILESConvertRequest, XYZValidateRequest,
    QuantumCalculationRequest, CalculationUpdateRequest
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize Flask app and extensions
app = Flask(__name__)
CORS(app)  # Enable CORS for cross-origin requests
sock = Sock(app)  # Sockインスタンスを作成

# Initialize PubChem client
pubchem_client = PubChemClient(timeout=30)



def run_calculation_in_background(calculation_id: str, parameters: dict):
    """
    Worker function to run the DFT calculation in a separate thread.
    Updates the calculation status and results directly on the file system.
    """
    file_manager = CalculationFileManager()
    calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)

    try:
        # Update status to running on the file system
        file_manager.save_calculation_status(calc_dir, 'running')
        
        # Initialize calculator
        calculator = DFTCalculator(working_dir=calc_dir, keep_files=True, molecule_name=parameters['name'])
        
        # Parse XYZ and setup calculation
        atoms = calculator.parse_xyz(parameters['xyz'])
        calculator.setup_calculation(
            atoms,
            basis=parameters['basis_function'],
            xc=parameters['exchange_correlation'],
            charge=parameters['charges'],
            spin=(parameters['spin_multiplicity'] - 1) // 2,
            max_cycle=150,
            solvent_method=parameters['solvent_method'],
            solvent=parameters['solvent']
        )
        
        # Run calculation
        results = calculator.run_calculation()
        
        # Save results and update status to completed
        file_manager.save_calculation_results(calc_dir, results)
        file_manager.save_calculation_status(calc_dir, 'completed')
        
        logger.info(f"Calculation {calculation_id} completed successfully.")

    except (InputError, ConvergenceError, CalculationError) as e:
        logger.error(f"Calculation {calculation_id} failed: {e}")
        file_manager.save_calculation_status(calc_dir, 'error')
        # エラー情報をresults.jsonに保存することも検討できる
        file_manager.save_calculation_results(calc_dir, {'error': str(e)})

    except Exception as e:
        logger.error(f"Unexpected error in calculation {calculation_id}: {e}", exc_info=True)
        file_manager.save_calculation_status(calc_dir, 'error')
        file_manager.save_calculation_results(calc_dir, {'error': 'An unexpected internal server error occurred.'})


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'ok',
        'service': 'pyscf-front-api',
        'version': '0.3.0'
    })


@app.route('/api/pubchem/search', methods=['POST'])
@validate()
def search_pubchem(body: PubChemSearchRequest):
    """Search PubChem for a compound and return its 3D structure in XYZ format."""
    try:
        query = body.query
        search_type = body.search_type.value
        
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


@app.route('/api/smiles/convert', methods=['POST'])
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


@app.route('/api/pubchem/validate', methods=['POST'])
@validate()
def validate_xyz_endpoint(body: XYZValidateRequest):
    """Validate an XYZ format string."""
    try:
        xyz_string = body.xyz
        
        validation_result = xyz_parser.validate_xyz(xyz_string)
        
        return jsonify({'success': True, 'data': validation_result})
        
    except Exception as e:
        logger.error(f"Error during XYZ validation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'An internal server error occurred.'}), 500


@app.route('/api/quantum/calculate', methods=['POST'])
@validate()
def quantum_calculate(body: QuantumCalculationRequest):
    """
    Starts a quantum chemistry calculation in the background.
    Immediately returns a calculation ID to track the job.
    """
    try:
        # Prepare parameters using validated data from Pydantic model
        parameters = {
            'calculation_method': body.calculation_method.value,
            'basis_function': body.basis_function,
            'exchange_correlation': body.exchange_correlation,
            'charges': body.charges,
            'spin_multiplicity': body.spin_multiplicity,
            'solvent_method': body.solvent_method.value,
            'solvent': body.solvent,
            'xyz': body.xyz,
            'name': body.name,
            'cpu_cores': body.cpu_cores,
            'memory_mb': body.memory_mb,
            'created_at': datetime.now().isoformat()
        }
        
        # Initialize file manager and create directory
        file_manager = CalculationFileManager()
        calc_dir = file_manager.create_calculation_dir(parameters['name'])
        calculation_id = os.path.basename(calc_dir)
        
        # Save initial parameters and status
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, 'pending')

        # Start background thread for the calculation
        thread = threading.Thread(
            target=run_calculation_in_background,
            args=(calculation_id, parameters)
        )
        thread.daemon = True
        thread.start()
        
        logger.info(f"Queued calculation {calculation_id} for molecule '{parameters['name']}'")

        # Return immediately with the new calculation instance
        initial_instance = {
            'id': calculation_id,
            'name': parameters['name'],
            'status': 'running', # Frontend polls for actual status, so this is fine
            'createdAt': parameters['created_at'],
            'updatedAt': parameters['created_at'],
            'parameters': parameters
        }

        return jsonify({'success': True, 'data': {'calculation': initial_instance}}), 202

    except Exception as e:
        logger.error(f"Error queuing calculation: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to queue calculation.'}), 500


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
    """
    Get detailed information about a specific calculation.
    This now exclusively reads from the file system, making it the single source of truth.
    """
    try:
        
        file_manager = CalculationFileManager()
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)

        if not os.path.isdir(calc_path):
             return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404

        # Read all calculation data from disk
        parameters = file_manager.read_calculation_parameters(calc_path) or {}
        results = file_manager.read_calculation_results(calc_path)
        status = file_manager.read_calculation_status(calc_path)
        
        display_name = file_manager._get_display_name(calculation_id, parameters)
        creation_date = parameters.get('created_at', datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat())

        calculation_instance = {
            'id': calculation_id,
            'name': display_name,
            'status': status,
            'createdAt': creation_date,
            'updatedAt': datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat(),
            'workingDirectory': calc_path,
            'parameters': parameters,
            'results': results
        }
        
        return jsonify({
            'success': True,
            'data': {
                'calculation': calculation_instance,
                'files': {
                    'checkpoint_exists': file_manager.file_exists(calc_path, 'calculation.chk'),
                    'parameters_file_exists': parameters is not None,
                    'results_file_exists': results is not None,
                }
            }
        })
        
    except Exception as e:
        logger.error(f"Error getting calculation details for {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to get calculation details.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['PUT'])
@validate()
def update_calculation(calculation_id, body: CalculationUpdateRequest):
    """Update calculation metadata (currently only name)."""
    try:
        file_manager = CalculationFileManager()
        
        new_name = body.name
        
        # Update the display name (calculation_id remains the same)
        result_id = file_manager.rename_calculation(calculation_id, new_name)
        if not result_id:
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404
        
        logger.info(f"Updated display name for calculation {calculation_id} to '{new_name}'")
        return jsonify({
            'success': True,
            'data': {
                'message': 'Calculation renamed successfully.',
                'name': new_name
            }
        })

    except Exception as e:
        logger.error(f"Error updating calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to update calculation.'}), 500


@app.route('/api/quantum/calculations/<calculation_id>', methods=['DELETE'])
def delete_calculation(calculation_id):
    """Delete a calculation and its files."""
    try:
        file_manager = CalculationFileManager()
        
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)

        if not os.path.isdir(calc_path):
            return jsonify({'success': False, 'error': f'Calculation "{calculation_id}" not found.'}), 404

        # Delete the calculation directory
        shutil.rmtree(calc_path)
        logger.info(f"Deleted calculation directory: {calc_path}")
        

        return jsonify({
            'success': True,
            'data': {
                'message': f'Calculation "{calculation_id}" has been deleted successfully',
                'deleted_id': calculation_id
            }
        })
        
    except Exception as e:
        logger.error(f"Error deleting calculation {calculation_id}: {e}", exc_info=True)
        return jsonify({'success': False, 'error': 'Failed to delete calculation.'}), 500


#
# ===== ここから新しいWebSocketエンドポイントを追加 =====
#
@sock.route('/ws/calculations/<calculation_id>')
def calculation_status_socket(ws, calculation_id):
    """
    WebSocketエンドポイント。計算のステータスを監視し、クライアントに送信する。
    """
    logger.info(f"WebSocket connection established for calculation {calculation_id}")
    file_manager = CalculationFileManager()
    calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
    last_known_status = None

    # 計算ディレクトリの存在確認
    if not os.path.isdir(calc_path):
        logger.warning(f"Calculation directory not found: {calc_path}")
        try:
            ws.send(json.dumps({
                'error': f'Calculation "{calculation_id}" not found.',
                'id': calculation_id
            }))
        except Exception as send_error:
            logger.error(f"Failed to send error message: {send_error}")
        finally:
            try:
                ws.close()
            except:
                pass
        return

    try:
        while True:
            try:
                # ファイルシステムから現在の状態を取得
                current_status = file_manager.read_calculation_status(calc_path)
                
                if current_status is None:
                    logger.warning(f"Could not read status for calculation {calculation_id}")
                    current_status = 'unknown'

                # ステータスが変更された場合のみデータを送信
                if current_status != last_known_status:
                    logger.info(f"Status changed for {calculation_id} to '{current_status}'. Sending update.")
                    
                    try:
                        # get_calculation_detailsと同様のロジックで詳細データを取得
                        parameters = file_manager.read_calculation_parameters(calc_path) or {}
                        results = file_manager.read_calculation_results(calc_path)
                        display_name = file_manager._get_display_name(calculation_id, parameters)
                        
                        # 作成日時の安全な取得
                        try:
                            creation_date = parameters.get('created_at', datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat())
                            updated_date = datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat()
                        except OSError:
                            current_time = datetime.now().isoformat()
                            creation_date = current_time
                            updated_date = current_time

                        # CalculationInstance と同じ形式のデータを作成
                        calculation_instance = {
                            'id': calculation_id,
                            'name': display_name,
                            'status': current_status,
                            'createdAt': creation_date,
                            'updatedAt': updated_date,
                            'parameters': parameters,
                            'results': results,
                            'workingDirectory': calc_path,
                        }
                        
                        # JSON送信をtry-catchで包む
                        try:
                            ws.send(json.dumps(calculation_instance))
                            last_known_status = current_status
                        except Exception as send_error:
                            logger.error(f"Failed to send WebSocket message: {send_error}")
                            break  # 送信に失敗した場合は接続を終了
                            
                    except Exception as data_error:
                        logger.error(f"Error preparing calculation data for {calculation_id}: {data_error}")
                        try:
                            ws.send(json.dumps({
                                'error': 'Failed to read calculation data',
                                'id': calculation_id,
                                'status': current_status
                            }))
                        except:
                            pass  # エラーメッセージ送信に失敗してもログ出力は行わない

                # 計算が完了またはエラーになったらループを終了
                if current_status in ['completed', 'error']:
                    logger.info(f"Calculation {calculation_id} finished with status '{current_status}'. Closing WebSocket.")
                    break
                
                # 1秒待機して次のチェックへ
                time.sleep(1)
                
            except Exception as loop_error:
                logger.error(f"Error in WebSocket loop for {calculation_id}: {loop_error}")
                break
            
    except Exception as e:
        logger.error(f"WebSocket error for {calculation_id}: {e}", exc_info=True)
    finally:
        # 接続が閉じることを確認
        logger.info(f"Closing WebSocket connection for calculation {calculation_id}")
        try:
            if hasattr(ws, 'close'):
                ws.close()
        except Exception as close_error:
            logger.debug(f"Error closing WebSocket: {close_error}")

#
# ===== WebSocketエンドポイントの追加はここまで =====
#


@app.errorhandler(ValidationError)
def validation_error_handler(error):
    """Handle Pydantic validation errors with consistent error format."""
    errors = error.errors()
    error_messages = []
    for err in errors:
        field = '.'.join(str(x) for x in err['loc'])
        message = err['msg']
        error_messages.append(f"{field}: {message}")
    
    combined_message = "Validation failed: " + "; ".join(error_messages)
    logger.warning(f"Validation error: {combined_message}")
    return jsonify({'success': False, 'error': combined_message}), 400

@app.errorhandler(404)
def not_found(error):
    return jsonify({'success': False, 'error': 'Endpoint not found.'}), 404

@app.errorhandler(405)
def method_not_allowed(error):
    return jsonify({'success': False, 'error': 'Method not allowed for this endpoint.'}), 405

if __name__ == '__main__':
    host = os.environ.get('FLASK_RUN_HOST', '127.0.0.1')
    # If a port is passed as a command-line argument, use it. Otherwise, let OS assign one.
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    debug = os.environ.get('FLASK_ENV') == 'development'
    
    # Create server with OS-assigned port (port=0)
    server = make_server(host, port, app)
    actual_port = server.server_port
    
    # Print the actual port to stdout so the Electron process can capture it.
    print(f"FLASK_SERVER_PORT:{actual_port}", file=sys.stdout, flush=True)
    
    logger.info(f"Starting API server on http://{host}:{actual_port} (Debug: {debug})")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        logger.info("Server stopped by user.")
        server.shutdown()