"""
WebSocket handlers for real-time calculation monitoring.
Handles client connections for calculation status updates and progress monitoring.
"""

import logging
import os
import threading
from datetime import datetime
from typing import Dict
from flask import session
from flask_socketio import emit, join_room, leave_room

from quantum_calc.file_manager import CalculationFileManager
from quantum_calc import get_websocket_watcher

# Set up logging
logger = logging.getLogger(__name__)

# Global WebSocket connection registry for immediate notifications  
# calculation_id -> set of session IDs
active_websockets = {}
websocket_lock = threading.Lock()


def send_immediate_websocket_notification(socketio, calculation_id: str, status: str, error_message: str = None):
    """Send immediate WebSocket notification to all connected clients for a calculation."""
    # Build complete calculation instance
    try:
        from quantum_calc import get_current_settings
        settings = get_current_settings()
        file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
        calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        if os.path.exists(calc_dir):
            # Read current data from files
            parameters = file_manager.read_calculation_parameters(calc_dir) or {}
            results = file_manager.read_calculation_results(calc_dir)
            display_name = file_manager._get_display_name(calculation_id, parameters)
            
            # Build calculation instance
            calculation_instance = {
                'id': calculation_id,
                'name': display_name,
                'status': status,
                'createdAt': parameters.get('created_at', datetime.now().isoformat()),
                'updatedAt': datetime.now().isoformat(),
                'parameters': parameters,
                'results': results,
                'workingDirectory': calc_dir,
            }
            
            # Add error if provided (両方のフィールドに設定してフロントエンドとの整合性を確保)
            if error_message:
                calculation_instance['error'] = error_message
                calculation_instance['errorMessage'] = error_message
            
            # Send to all clients in the calculation room
            socketio.emit('calculation_update', calculation_instance, room=f'calculation_{calculation_id}')
            
            # Also send to global updates room for non-active calculations monitoring
            socketio.emit('calculation_update', calculation_instance, room='global_updates')
            logger.debug(f"Sent immediate notification for calculation {calculation_id} with status {status} to both specific and global rooms")
        else:
            logger.warning(f"Calculation directory not found for immediate notification: {calc_dir}")
            
    except Exception as e:
        logger.error(f"Error sending immediate WebSocket notification for {calculation_id}: {e}")


def build_calculation_instance(calc_id: str, calc_path: str, file_manager: CalculationFileManager) -> Dict:
    """Build a complete calculation instance from file system data."""
    try:
        parameters = file_manager.read_calculation_parameters(calc_path) or {}
        results = file_manager.read_calculation_results(calc_path)
        status = file_manager.read_calculation_status(calc_path) or 'pending'
        display_name = file_manager._get_display_name(calc_id, parameters)
        
        # Safe date retrieval
        try:
            creation_date = parameters.get('created_at', 
                datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat())
            updated_date = datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat()
        except OSError:
            current_time = datetime.now().isoformat()
            creation_date = current_time
            updated_date = current_time

        return {
            'id': calc_id,
            'name': display_name,
            'status': status,
            'createdAt': creation_date,
            'updatedAt': updated_date,
            'parameters': parameters,
            'results': results,
            'workingDirectory': calc_path,
        }
    except Exception as e:
        logger.error(f"Error building calculation instance for {calc_id}: {e}")
        raise


def register_websocket_handlers(socketio):
    """Register all WebSocket event handlers with the SocketIO instance."""
    
    @socketio.on('join_calculation')
    def on_join_calculation(data):
        """Join a calculation room to receive real-time updates."""
        calculation_id = data.get('calculation_id')
        if not calculation_id:
            emit('error', {'error': 'calculation_id is required'})
            return
        
        # 一時的IDの場合は特別なログメッセージを出力
        if calculation_id.startswith('new-calculation-'):
            logger.info(f"SocketIO connection attempt for temporary calculation ID: {calculation_id}")
        else:
            logger.info(f"SocketIO connection established for calculation {calculation_id}")

        from quantum_calc import get_current_settings
        settings = get_current_settings()
        file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
        calc_path = os.path.join(file_manager.get_base_directory(), calculation_id)
        
        # Calculation directory existence check
        if not os.path.isdir(calc_path):
            # 一時的IDの場合は特別なログメッセージを出力
            if calculation_id.startswith('new-calculation-'):
                logger.info(f"SocketIO connection attempted for temporary calculation ID: {calculation_id}. Sending error.")
                error_message = f'Temporary calculation ID "{calculation_id}" does not exist on server.'
            else:
                logger.warning(f"Calculation directory not found: {calc_path}")
                error_message = f'Calculation "{calculation_id}" not found.'
            
            emit('error', {
                'error': error_message,
                'id': calculation_id,
                'is_temporary': calculation_id.startswith('new-calculation-')
            })
            return

        # Join the calculation room
        room = f'calculation_{calculation_id}'
        join_room(room)
        
        # Store calculation_id in session for cleanup
        session['calculation_id'] = calculation_id
        
        # Define file change callback for this connection
        def on_file_change(file_data: Dict):
            """Callback for file system events - sends updated data to SocketIO client."""
            try:
                # Build complete calculation instance
                calculation_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
                
                # Send update to all clients in the room
                socketio.emit('calculation_update', calculation_instance, room=room)
                
                # Disconnect clients if calculation is finished
                if calculation_instance['status'] in ['completed', 'error']:
                    logger.info(f"Calculation {calculation_id} finished with status '{calculation_instance['status']}'.")
                    # Note: Don't force disconnect - let client handle completion
                    
            except Exception as e:
                logger.error(f"Error in file change callback for {calculation_id}: {e}")
                socketio.emit('error', {
                    'error': 'Failed to read calculation data',
                    'id': calculation_id
                }, room=room)

        try:
            # Initialize file watcher and register this connection
            watcher = get_websocket_watcher(file_manager.get_base_directory())
            watcher.add_connection(calculation_id, on_file_change)
            
            # Send initial state immediately
            initial_instance = build_calculation_instance(calculation_id, calc_path, file_manager)
            emit('calculation_update', initial_instance)
            logger.info(f"Sent initial state for calculation {calculation_id} (status: {initial_instance['status']})")
            
            # Store callback reference for cleanup
            session['file_change_callback'] = on_file_change
            
        except Exception as e:
            logger.error(f"Error setting up SocketIO monitoring for {calculation_id}: {e}")
            emit('error', {
                'error': 'Failed to set up calculation monitoring',
                'id': calculation_id
            })

    @socketio.on('leave_calculation')
    def on_leave_calculation(data):
        """Leave a calculation room."""
        calculation_id = data.get('calculation_id')
        if not calculation_id:
            return

        room = f'calculation_{calculation_id}'
        leave_room(room)

        # Clean up file watcher connection
        if 'file_change_callback' in session:
            try:
                from quantum_calc import get_current_settings
                settings = get_current_settings()
                file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
                watcher = get_websocket_watcher(file_manager.get_base_directory())
                watcher.remove_connection(calculation_id, session['file_change_callback'])
            except Exception as e:
                logger.debug(f"Error cleaning up file watcher for {calculation_id}: {e}")
        
        logger.info(f"Client left calculation {calculation_id}")

    @socketio.on('join_global_updates')
    def on_join_global_updates():
        """Join global updates room to receive all calculation updates."""
        join_room('global_updates')
        logger.info("Client joined global_updates room for real-time monitoring of all calculations")

    @socketio.on('leave_global_updates')
    def on_leave_global_updates():
        """Leave global updates room."""
        leave_room('global_updates')
        logger.info("Client left global_updates room")

    @socketio.on('disconnect')
    def on_disconnect():
        """Handle client disconnection."""
        calculation_id = session.get('calculation_id')
        if calculation_id and 'file_change_callback' in session:
            try:
                from quantum_calc import get_current_settings
                settings = get_current_settings()
                file_manager = CalculationFileManager(base_dir=settings.calculations_directory)
                watcher = get_websocket_watcher(file_manager.get_base_directory())
                watcher.remove_connection(calculation_id, session['file_change_callback'])
                logger.info(f"Cleaned up file watcher for disconnected client (calculation {calculation_id})")
            except Exception as e:
                logger.debug(f"Error cleaning up file watcher on disconnect for {calculation_id}: {e}")

    return {
        'send_immediate_websocket_notification': lambda calc_id, status, error_msg=None: send_immediate_websocket_notification(socketio, calc_id, status, error_msg)
    }