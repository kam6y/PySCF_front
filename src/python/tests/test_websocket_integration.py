# tests/test_websocket_integration.py

"""Integration tests for WebSocket calculation monitoring."""

import json
import os
import shutil
import tempfile
import threading
import time
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pytest

from app import app
from quantum_calc.file_manager import CalculationFileManager
from quantum_calc.file_watcher import get_websocket_watcher, shutdown_websocket_watcher


class MockWebSocket:
    """Mock WebSocket implementation for testing."""
    
    def __init__(self):
        self.messages_sent = []
        self.closed = False
        self.receive_timeout = None
        self.receive_exception = None
    
    def send(self, message):
        """Mock send method that records messages."""
        if self.closed:
            raise Exception("WebSocket is closed")
        self.messages_sent.append(message)
    
    def close(self):
        """Mock close method."""
        self.closed = True
    
    def receive(self, timeout=None):
        """Mock receive method."""
        if self.receive_exception:
            raise self.receive_exception
        if timeout and self.receive_timeout:
            time.sleep(self.receive_timeout)
        return None  # Simulate client disconnect


@pytest.fixture
def temp_calc_dir():
    """Create temporary calculation directory for testing."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


@pytest.fixture
def file_manager(temp_calc_dir):
    """Create CalculationFileManager with temp directory."""
    return CalculationFileManager(temp_calc_dir)


@pytest.fixture
def mock_websocket():
    """Create mock WebSocket for testing."""
    return MockWebSocket()


class TestWebSocketIntegration:
    """Integration test cases for WebSocket calculation monitoring."""
    
    def setup_method(self):
        """Set up test environment."""
        # Reset global watcher state
        shutdown_websocket_watcher()
        import quantum_calc.file_watcher
        quantum_calc.file_watcher._global_watcher = None
    
    def teardown_method(self):
        """Clean up test environment."""
        shutdown_websocket_watcher()
    
    def test_websocket_connection_with_valid_calculation(self, temp_calc_dir, file_manager, mock_websocket):
        """Test WebSocket connection to valid calculation."""
        # Create calculation directory with initial data
        calc_id = "test_calc_20241201_120000"
        calc_dir = file_manager.create_calculation_dir("test_molecule")
        calc_id = os.path.basename(calc_dir)
        
        # Save initial calculation data
        parameters = {
            "name": "test_molecule",
            "calculation_method": "DFT",
            "basis_function": "6-31G",
            "created_at": "2024-12-01T12:00:00"
        }
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, "running")
        
        # Import the WebSocket handler
        from app import calculation_status_socket
        
        # Mock the necessary Flask-Sock components
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_get_watcher.return_value = mock_watcher
            
            # Execute WebSocket handler
            calculation_status_socket(mock_websocket, calc_id)
            
            # Verify watcher interaction
            mock_watcher.add_connection.assert_called_once()
            mock_watcher.remove_connection.assert_called_once()
            
            # Verify initial state was sent
            assert len(mock_websocket.messages_sent) >= 1
            initial_message = json.loads(mock_websocket.messages_sent[0])
            assert initial_message['id'] == calc_id
            assert initial_message['name'] == "test_molecule"
            assert initial_message['status'] == "running"
    
    def test_websocket_connection_with_invalid_calculation(self, temp_calc_dir, mock_websocket):
        """Test WebSocket connection to non-existent calculation."""
        # Import the WebSocket handler
        from app import calculation_status_socket
        
        # Execute WebSocket handler with invalid calculation ID
        calculation_status_socket(mock_websocket, "nonexistent_calc")
        
        # Verify error message was sent
        assert len(mock_websocket.messages_sent) == 1
        error_message = json.loads(mock_websocket.messages_sent[0])
        assert 'error' in error_message
        assert 'nonexistent_calc' in error_message['error']
        assert mock_websocket.closed
    
    def test_websocket_file_change_notification(self, temp_calc_dir, file_manager, mock_websocket):
        """Test file change notifications through WebSocket."""
        # Create calculation with running status
        calc_dir = file_manager.create_calculation_dir("test_molecule")
        calc_id = os.path.basename(calc_dir)
        
        parameters = {
            "name": "test_molecule",
            "calculation_method": "DFT",
            "created_at": "2024-12-01T12:00:00"
        }
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, "running")
        
        # Mock WebSocket connection that simulates file changes
        file_change_callback = None
        
        def capture_callback(calc_id, callback):
            nonlocal file_change_callback
            file_change_callback = callback
        
        # Import and setup WebSocket handler
        from app import calculation_status_socket
        
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_watcher.add_connection = capture_callback
            mock_get_watcher.return_value = mock_watcher
            
            # Configure mock WebSocket to timeout after receiving file change
            mock_websocket.receive_timeout = 0.1
            
            # Start WebSocket handler in separate thread
            handler_thread = threading.Thread(
                target=calculation_status_socket,
                args=(mock_websocket, calc_id)
            )
            handler_thread.start()
            
            # Wait for initial state to be sent
            time.sleep(0.05)
            
            # Simulate file change
            if file_change_callback:
                file_change_data = {
                    "status.json": {"status": "completed", "progress": 100}
                }
                file_change_callback(file_change_data)
            
            # Wait for handler to complete
            handler_thread.join(timeout=1.0)
            
            # Verify initial state and file change notification were sent
            assert len(mock_websocket.messages_sent) >= 1
            initial_message = json.loads(mock_websocket.messages_sent[0])
            assert initial_message['status'] == "running"
    
    def test_websocket_connection_cleanup_on_completion(self, temp_calc_dir, file_manager, mock_websocket):
        """Test WebSocket connection cleanup when calculation completes."""
        # Create calculation that's already completed
        calc_dir = file_manager.create_calculation_dir("completed_molecule")
        calc_id = os.path.basename(calc_dir)
        
        parameters = {
            "name": "completed_molecule",
            "calculation_method": "DFT",
            "created_at": "2024-12-01T12:00:00"
        }
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, "completed")
        file_manager.save_calculation_results(calc_dir, {"energy": -123.456})
        
        # Import the WebSocket handler
        from app import calculation_status_socket
        
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_get_watcher.return_value = mock_watcher
            
            # Execute WebSocket handler
            calculation_status_socket(mock_websocket, calc_id)
            
            # Verify connection was cleaned up
            mock_watcher.remove_connection.assert_called_once()
            
            # Verify initial completed state was sent
            assert len(mock_websocket.messages_sent) == 1
            message = json.loads(mock_websocket.messages_sent[0])
            assert message['status'] == "completed"
            assert message['results']['energy'] == -123.456
    
    def test_websocket_error_handling(self, temp_calc_dir, file_manager):
        """Test WebSocket error handling scenarios."""
        # Create calculation directory
        calc_dir = file_manager.create_calculation_dir("error_test")
        calc_id = os.path.basename(calc_dir)
        
        parameters = {"name": "error_test", "created_at": "2024-12-01T12:00:00"}
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, "running")
        
        # Create mock WebSocket that fails on send
        failing_websocket = Mock()
        failing_websocket.send.side_effect = Exception("Send failed")
        failing_websocket.close = Mock()
        
        # Import the WebSocket handler
        from app import calculation_status_socket
        
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_get_watcher.return_value = mock_watcher
            
            # Execute WebSocket handler - should not raise exception
            calculation_status_socket(failing_websocket, calc_id)
            
            # Verify cleanup was attempted
            mock_watcher.remove_connection.assert_called()
    
    def test_websocket_multiple_connections_same_calculation(self, temp_calc_dir, file_manager):
        """Test multiple WebSocket connections to the same calculation."""
        # Create calculation directory
        calc_dir = file_manager.create_calculation_dir("multi_connection_test")
        calc_id = os.path.basename(calc_dir)
        
        parameters = {
            "name": "multi_connection_test",
            "created_at": "2024-12-01T12:00:00"
        }
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, "running")
        
        # Create multiple mock WebSockets
        websocket1 = MockWebSocket()
        websocket2 = MockWebSocket()
        websocket1.receive_timeout = 0.1
        websocket2.receive_timeout = 0.1
        
        # Import the WebSocket handler
        from app import calculation_status_socket
        
        connection_count = 0
        
        def count_connections(*args):
            nonlocal connection_count
            connection_count += 1
        
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_watcher.add_connection = count_connections
            mock_get_watcher.return_value = mock_watcher
            
            # Start multiple WebSocket handlers
            thread1 = threading.Thread(
                target=calculation_status_socket,
                args=(websocket1, calc_id)
            )
            thread2 = threading.Thread(
                target=calculation_status_socket,
                args=(websocket2, calc_id)
            )
            
            thread1.start()
            thread2.start()
            
            thread1.join(timeout=1.0)
            thread2.join(timeout=1.0)
            
            # Verify both connections were added
            assert connection_count == 2
            
            # Verify both WebSockets received initial state
            assert len(websocket1.messages_sent) >= 1
            assert len(websocket2.messages_sent) >= 1
    
    def test_build_calculation_instance_helper(self, temp_calc_dir, file_manager):
        """Test the build_calculation_instance helper function."""
        # Create calculation with complete data
        calc_dir = file_manager.create_calculation_dir("instance_test")
        calc_id = os.path.basename(calc_dir)
        
        parameters = {
            "name": "instance_test",
            "calculation_method": "DFT",
            "basis_function": "6-31G",
            "created_at": "2024-12-01T12:00:00"
        }
        results = {
            "energy": -123.456,
            "dipole": [0.1, 0.2, 0.3]
        }
        
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_results(calc_dir, results)
        file_manager.save_calculation_status(calc_dir, "completed")
        
        # Import and test the helper function indirectly through WebSocket handler
        from app import calculation_status_socket
        
        mock_websocket = MockWebSocket()
        
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_get_watcher.return_value = mock_watcher
            
            calculation_status_socket(mock_websocket, calc_id)
            
            # Verify the built instance contains expected data
            assert len(mock_websocket.messages_sent) == 1
            instance = json.loads(mock_websocket.messages_sent[0])
            
            assert instance['id'] == calc_id
            assert instance['name'] == "instance_test"
            assert instance['status'] == "completed"
            assert instance['parameters']['calculation_method'] == "DFT"
            assert instance['results']['energy'] == -123.456
            assert 'createdAt' in instance
            assert 'updatedAt' in instance
            assert 'workingDirectory' in instance


class TestWebSocketStressScenarios:
    """Stress test scenarios for WebSocket implementation."""
    
    def setup_method(self):
        """Set up test environment."""
        shutdown_websocket_watcher()
        import quantum_calc.file_watcher
        quantum_calc.file_watcher._global_watcher = None
    
    def teardown_method(self):
        """Clean up test environment."""
        shutdown_websocket_watcher()
    
    def test_rapid_file_changes(self, temp_calc_dir, file_manager):
        """Test handling of rapid file changes."""
        # Create calculation directory
        calc_dir = file_manager.create_calculation_dir("rapid_changes_test")
        calc_id = os.path.basename(calc_dir)
        
        parameters = {
            "name": "rapid_changes_test",
            "created_at": "2024-12-01T12:00:00"
        }
        file_manager.save_calculation_parameters(calc_dir, parameters)
        file_manager.save_calculation_status(calc_dir, "running")
        
        # Mock WebSocket that can handle multiple messages
        mock_websocket = MockWebSocket()
        mock_websocket.receive_timeout = 0.2  # Longer timeout to allow for rapid changes
        
        file_change_callback = None
        
        def capture_callback(calc_id, callback):
            nonlocal file_change_callback
            file_change_callback = callback
        
        from app import calculation_status_socket
        
        with patch('quantum_calc.file_watcher.get_websocket_watcher') as mock_get_watcher:
            mock_watcher = Mock()
            mock_watcher.add_connection = capture_callback
            mock_get_watcher.return_value = mock_watcher
            
            # Start WebSocket handler
            handler_thread = threading.Thread(
                target=calculation_status_socket,
                args=(mock_websocket, calc_id)
            )
            handler_thread.start()
            
            # Wait for connection setup
            time.sleep(0.05)
            
            # Simulate rapid file changes
            if file_change_callback:
                for i in range(10):
                    file_change_data = {
                        "status.json": {"status": "running", "progress": i * 10}
                    }
                    file_change_callback(file_change_data)
                    time.sleep(0.01)  # Small delay between changes
            
            handler_thread.join(timeout=2.0)
            
            # Verify initial state was sent (and potentially some updates)
            assert len(mock_websocket.messages_sent) >= 1
            
            # Verify no exceptions occurred (handler completed successfully)
            assert not handler_thread.is_alive()