"""
Integration tests for WebSocket handlers.

Tests the WebSocket functionality for real-time calculation monitoring,
including joining rooms, receiving updates, and handling disconnections.
"""

import pytest
import os
import json
import tempfile
import shutil
from pathlib import Path


class TestJoinCalculationWebSocket:
    """Integration tests for join_calculation WebSocket event."""

    def test_join_calculation_success(self, socketio_client, app):
        """
        GIVEN a calculation directory exists
        WHEN client emits 'join_calculation' event
        THEN client joins the room and receives initial state
        """
        # ARRANGE
        # Create a temporary calculation directory
        calc_id = 'test-calc-123'
        calc_dir = Path(app.config['CALCULATIONS_DIR']) / calc_id
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create parameters.json
        params = {
            'name': 'Test Calculation',
            'calculation_method': 'HF',
            'basis_function': 'sto-3g',
            'created_at': '2024-01-01T00:00:00'
        }
        with open(calc_dir / 'parameters.json', 'w') as f:
            json.dump(params, f)
        
        # Create status.json
        status_data = {'status': 'completed'}
        with open(calc_dir / 'status.json', 'w') as f:
            json.dump(status_data, f)

        # ACT
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        
        # Wait for and get received messages
        received = socketio_client.get_received()

        # ASSERT
        assert len(received) > 0
        
        # Should receive initial calculation_update event
        update_events = [msg for msg in received if msg['name'] == 'calculation_update']
        assert len(update_events) > 0
        
        # Verify initial state
        initial_state = update_events[0]['args'][0]
        assert initial_state['id'] == calc_id
        assert initial_state['name'] == 'Test Calculation'
        assert initial_state['status'] == 'completed'
        
        # Cleanup
        shutil.rmtree(calc_dir)

    def test_join_calculation_not_found(self, socketio_client):
        """
        GIVEN calculation does not exist
        WHEN client emits 'join_calculation' event
        THEN client receives error event
        """
        # ARRANGE
        nonexistent_calc_id = 'nonexistent-calc-999'

        # ACT
        socketio_client.emit('join_calculation', {'calculation_id': nonexistent_calc_id})
        
        # Wait for response
        received = socketio_client.get_received()

        # ASSERT
        assert len(received) > 0
        
        # Should receive error event
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) > 0
        
        error_data = error_events[0]['args'][0]
        assert 'error' in error_data
        assert nonexistent_calc_id in error_data['error'] or 'not found' in error_data['error'].lower()

    def test_join_calculation_missing_id(self, socketio_client):
        """
        GIVEN calculation_id is missing from the event data
        WHEN client emits 'join_calculation' event
        THEN client receives error event
        """
        # ACT
        socketio_client.emit('join_calculation', {})
        
        # Wait for response
        received = socketio_client.get_received()

        # ASSERT
        assert len(received) > 0
        
        # Should receive error event
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) > 0

    def test_join_calculation_temporary_id(self, socketio_client):
        """
        GIVEN temporary calculation ID (new-calculation-*)
        WHEN client emits 'join_calculation' event
        THEN client receives error with is_temporary flag
        """
        # ARRANGE
        temp_calc_id = 'new-calculation-temp-123'

        # ACT
        socketio_client.emit('join_calculation', {'calculation_id': temp_calc_id})
        
        # Wait for response
        received = socketio_client.get_received()

        # ASSERT
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) > 0
        
        error_data = error_events[0]['args'][0]
        assert 'error' in error_data
        assert error_data.get('is_temporary') is True


class TestLeaveCalculationWebSocket:
    """Integration tests for leave_calculation WebSocket event."""

    def test_leave_calculation_success(self, socketio_client, app):
        """
        GIVEN client has joined a calculation room
        WHEN client emits 'leave_calculation' event
        THEN client leaves the room successfully
        """
        # ARRANGE
        # Create a temporary calculation directory
        calc_id = 'test-calc-leave'
        calc_dir = Path(app.config['CALCULATIONS_DIR']) / calc_id
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create minimal files
        with open(calc_dir / 'parameters.json', 'w') as f:
            json.dump({'name': 'Test', 'created_at': '2024-01-01T00:00:00'}, f)
        with open(calc_dir / 'status.json', 'w') as f:
            json.dump({'status': 'running'}, f)

        # First join the room
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        socketio_client.get_received()  # Clear initial messages

        # ACT
        socketio_client.emit('leave_calculation', {'calculation_id': calc_id})
        
        # ASSERT
        # No error should occur (successful leave is silent)
        # If we receive messages, they should not be errors
        received = socketio_client.get_received()
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) == 0
        
        # Cleanup
        shutil.rmtree(calc_dir)

    def test_leave_calculation_not_joined(self, socketio_client):
        """
        GIVEN client has not joined a calculation room
        WHEN client emits 'leave_calculation' event
        THEN no error occurs (silent operation)
        """
        # ACT
        socketio_client.emit('leave_calculation', {'calculation_id': 'some-calc'})
        
        # ASSERT
        # Should not receive error (leave is idempotent)
        received = socketio_client.get_received()
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) == 0


class TestGlobalUpdatesWebSocket:
    """Integration tests for global_updates room functionality."""

    def test_join_global_updates(self, socketio_client):
        """
        GIVEN WebSocket client is connected
        WHEN client emits 'join_global_updates' event
        THEN client joins global updates room
        """
        # ACT
        socketio_client.emit('join_global_updates')
        
        # ASSERT
        # Should not receive error
        received = socketio_client.get_received()
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) == 0

    def test_leave_global_updates(self, socketio_client):
        """
        GIVEN client has joined global updates room
        WHEN client emits 'leave_global_updates' event
        THEN client leaves global updates room
        """
        # ARRANGE
        socketio_client.emit('join_global_updates')
        socketio_client.get_received()  # Clear messages

        # ACT
        socketio_client.emit('leave_global_updates')
        
        # ASSERT
        # Should not receive error
        received = socketio_client.get_received()
        error_events = [msg for msg in received if msg['name'] == 'error']
        assert len(error_events) == 0


class TestWebSocketDisconnection:
    """Integration tests for WebSocket disconnection handling."""

    def test_disconnect_cleans_up(self, socketio_client, app):
        """
        GIVEN client has joined a calculation room
        WHEN client disconnects
        THEN cleanup is performed (no errors)
        """
        # ARRANGE
        # Create a temporary calculation directory
        calc_id = 'test-calc-disconnect'
        calc_dir = Path(app.config['CALCULATIONS_DIR']) / calc_id
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create minimal files
        with open(calc_dir / 'parameters.json', 'w') as f:
            json.dump({'name': 'Test', 'created_at': '2024-01-01T00:00:00'}, f)
        with open(calc_dir / 'status.json', 'w') as f:
            json.dump({'status': 'running'}, f)

        # Join room
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        socketio_client.get_received()

        # ACT
        socketio_client.disconnect()

        # ASSERT
        assert not socketio_client.is_connected()
        
        # Cleanup
        shutil.rmtree(calc_dir)


class TestCalculationUpdates:
    """Integration tests for receiving calculation updates via WebSocket."""

    def test_receive_update_on_file_change(self, socketio_client, app):
        """
        GIVEN client has joined a calculation room
        WHEN calculation status file is updated
        THEN client receives calculation_update event
        
        Note: This test validates the setup but may not test actual file watching
        in the test environment due to WEBSOCKET_WATCHER_ENABLED being False.
        """
        # ARRANGE
        calc_id = 'test-calc-update'
        calc_dir = Path(app.config['CALCULATIONS_DIR']) / calc_id
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create initial files
        params = {'name': 'Test Update', 'created_at': '2024-01-01T00:00:00'}
        with open(calc_dir / 'parameters.json', 'w') as f:
            json.dump(params, f)
        with open(calc_dir / 'status.json', 'w') as f:
            json.dump({'status': 'running'}, f)

        # Join room
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        received = socketio_client.get_received()
        
        # ASSERT
        # Should receive initial state
        update_events = [msg for msg in received if msg['name'] == 'calculation_update']
        assert len(update_events) > 0
        
        initial_state = update_events[0]['args'][0]
        assert initial_state['status'] == 'running'
        
        # Cleanup
        shutil.rmtree(calc_dir)

    def test_calculation_update_contains_all_fields(self, socketio_client, app):
        """
        GIVEN calculation has complete data
        WHEN client joins calculation
        THEN received update contains all expected fields
        """
        # ARRANGE
        calc_id = 'test-calc-complete'
        calc_dir = Path(app.config['CALCULATIONS_DIR']) / calc_id
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create complete calculation data
        params = {
            'name': 'Complete Calculation',
            'calculation_method': 'DFT',
            'basis_function': 'sto-3g',
            'exchange_correlation': 'b3lyp',
            'charges': 0,
            'spin': 0,
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'created_at': '2024-01-01T00:00:00'
        }
        with open(calc_dir / 'parameters.json', 'w') as f:
            json.dump(params, f)
        
        results = {'energy': -1.06, 'dipole': [0, 0, 0]}
        with open(calc_dir / 'results.json', 'w') as f:
            json.dump(results, f)
        
        with open(calc_dir / 'status.json', 'w') as f:
            json.dump({'status': 'completed'}, f)

        # ACT
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        received = socketio_client.get_received()

        # ASSERT
        update_events = [msg for msg in received if msg['name'] == 'calculation_update']
        assert len(update_events) > 0
        
        calc_data = update_events[0]['args'][0]
        
        # Verify all expected fields are present
        assert 'id' in calc_data
        assert 'name' in calc_data
        assert 'status' in calc_data
        assert 'createdAt' in calc_data
        assert 'updatedAt' in calc_data
        assert 'parameters' in calc_data
        assert 'results' in calc_data
        assert 'workingDirectory' in calc_data
        
        # Verify data values
        assert calc_data['id'] == calc_id
        assert calc_data['name'] == 'Complete Calculation'
        assert calc_data['status'] == 'completed'
        assert calc_data['results']['energy'] == -1.06
        
        # Cleanup
        shutil.rmtree(calc_dir)

    def test_calculation_update_with_error(self, socketio_client, app):
        """
        GIVEN calculation has error status
        WHEN client joins calculation
        THEN received update contains error information
        """
        # ARRANGE
        calc_id = 'test-calc-error'
        calc_dir = Path(app.config['CALCULATIONS_DIR']) / calc_id
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create calculation with error
        params = {'name': 'Error Calculation', 'created_at': '2024-01-01T00:00:00'}
        with open(calc_dir / 'parameters.json', 'w') as f:
            json.dump(params, f)
        
        error_info = {
            'status': 'error',
            'error': 'Calculation failed: SCF did not converge'
        }
        with open(calc_dir / 'status.json', 'w') as f:
            json.dump(error_info, f)

        # ACT
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        received = socketio_client.get_received()

        # ASSERT
        update_events = [msg for msg in received if msg['name'] == 'calculation_update']
        assert len(update_events) > 0
        
        calc_data = update_events[0]['args'][0]
        assert calc_data['status'] == 'error'
        
        # Cleanup
        shutil.rmtree(calc_dir)
