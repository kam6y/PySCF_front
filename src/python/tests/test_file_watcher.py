# tests/test_file_watcher.py

"""Unit tests for file_watcher module."""

import json
import os
import shutil
import tempfile
import time
import threading
from pathlib import Path
from unittest.mock import Mock, patch, call
import pytest

from quantum_calc.file_watcher import (
    CalculationFileWatcher,
    WebSocketCalculationWatcher,
    get_websocket_watcher,
    shutdown_websocket_watcher
)
from quantum_calc.exceptions import FileManagerError, WebSocketError


class TestCalculationFileWatcher:
    """Test cases for CalculationFileWatcher class."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.callback_mock = Mock()
        self.watcher = CalculationFileWatcher(self.callback_mock)
    
    def teardown_method(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test proper initialization of CalculationFileWatcher."""
        assert self.watcher.callback == self.callback_mock
        assert self.watcher.monitored_files == {'status.json', 'results.json', 'parameters.json'}
        assert hasattr(self.watcher, '_lock')
    
    def test_handles_monitored_files_only(self):
        """Test that only monitored files trigger callbacks."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        # Create a non-monitored file
        non_monitored_file = calc_dir / "other_file.txt"
        non_monitored_file.write_text("test")
        
        # Simulate file modification event
        from watchdog.events import FileModifiedEvent
        event = FileModifiedEvent(str(non_monitored_file))
        self.watcher.on_modified(event)
        
        # Callback should not be called
        self.callback_mock.assert_not_called()
    
    def test_handles_status_json_change(self):
        """Test handling of status.json file changes."""
        # Create calculation directory and status file
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_data = {"status": "running", "progress": 50}
        status_file.write_text(json.dumps(status_data))
        
        # Simulate file modification event
        from watchdog.events import FileModifiedEvent
        event = FileModifiedEvent(str(status_file))
        self.watcher.on_modified(event)
        
        # Verify callback was called with correct data
        self.callback_mock.assert_called_once_with("test_calc_001", {"status.json": status_data})
    
    def test_handles_results_json_change(self):
        """Test handling of results.json file changes."""
        # Create calculation directory and results file
        calc_dir = Path(self.temp_dir) / "test_calc_002"
        calc_dir.mkdir()
        
        results_file = calc_dir / "results.json"
        results_data = {"energy": -123.456, "dipole": [0.1, 0.2, 0.3]}
        results_file.write_text(json.dumps(results_data))
        
        # Simulate file creation event
        from watchdog.events import FileCreatedEvent
        event = FileCreatedEvent(str(results_file))
        self.watcher.on_created(event)
        
        # Verify callback was called with correct data
        self.callback_mock.assert_called_once_with("test_calc_002", {"results.json": results_data})
    
    def test_handles_invalid_json_gracefully(self):
        """Test graceful handling of invalid JSON files."""
        # Create calculation directory and invalid JSON file
        calc_dir = Path(self.temp_dir) / "test_calc_003"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_file.write_text("invalid json content")
        
        # Simulate file modification event
        from watchdog.events import FileModifiedEvent
        event = FileModifiedEvent(str(status_file))
        self.watcher.on_modified(event)
        
        # Callback should not be called for invalid JSON
        self.callback_mock.assert_not_called()
    
    def test_handles_missing_file_gracefully(self):
        """Test graceful handling when file doesn't exist."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_004"
        calc_dir.mkdir()
        
        # Reference non-existent status file
        status_file = calc_dir / "status.json"
        
        # Simulate file modification event for non-existent file
        from watchdog.events import FileModifiedEvent
        event = FileModifiedEvent(str(status_file))
        self.watcher.on_modified(event)
        
        # Callback should not be called for missing file
        self.callback_mock.assert_not_called()
    
    def test_ignores_directory_events(self):
        """Test that directory events are ignored."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_005"
        calc_dir.mkdir()
        
        # Simulate directory modification event
        from watchdog.events import FileModifiedEvent
        event = FileModifiedEvent(str(calc_dir))
        event.is_directory = True
        self.watcher.on_modified(event)
        
        # Callback should not be called for directory events
        self.callback_mock.assert_not_called()


class TestWebSocketCalculationWatcher:
    """Test cases for WebSocketCalculationWatcher class."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.watcher = WebSocketCalculationWatcher(self.temp_dir)
    
    def teardown_method(self):
        """Clean up test environment."""
        self.watcher.stop()
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test proper initialization of WebSocketCalculationWatcher."""
        assert self.watcher.base_directory == Path(self.temp_dir)
        assert len(self.watcher.connections) == 0
        assert len(self.watcher.watched_dirs) == 0
        assert not self.watcher._started
    
    def test_start_and_stop(self):
        """Test starting and stopping the file observer."""
        # Start the watcher
        self.watcher.start()
        assert self.watcher._started
        
        # Stop the watcher
        self.watcher.stop()
        assert not self.watcher._started
        assert len(self.watcher.connections) == 0
        assert len(self.watcher.watched_dirs) == 0
    
    def test_add_connection(self):
        """Test adding WebSocket connections."""
        self.watcher.start()
        
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        callback_mock = Mock()
        
        # Add connection
        self.watcher.add_connection("test_calc_001", callback_mock)
        
        # Verify connection was added
        assert "test_calc_001" in self.watcher.connections
        assert callback_mock in self.watcher.connections["test_calc_001"]
        assert str(calc_dir) in self.watcher.watched_dirs
    
    def test_remove_connection(self):
        """Test removing WebSocket connections."""
        self.watcher.start()
        
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        callback_mock = Mock()
        
        # Add and then remove connection
        self.watcher.add_connection("test_calc_001", callback_mock)
        self.watcher.remove_connection("test_calc_001", callback_mock)
        
        # Verify connection was removed
        assert "test_calc_001" not in self.watcher.connections
        assert str(calc_dir) not in self.watcher.watched_dirs
    
    def test_multiple_connections_same_calculation(self):
        """Test multiple WebSocket connections for the same calculation."""
        self.watcher.start()
        
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        callback1 = Mock()
        callback2 = Mock()
        
        # Add multiple connections
        self.watcher.add_connection("test_calc_001", callback1)
        self.watcher.add_connection("test_calc_001", callback2)
        
        # Verify both connections are registered
        assert len(self.watcher.connections["test_calc_001"]) == 2
        assert callback1 in self.watcher.connections["test_calc_001"]
        assert callback2 in self.watcher.connections["test_calc_001"]
        
        # Remove one connection
        self.watcher.remove_connection("test_calc_001", callback1)
        
        # Verify only one connection remains and directory is still watched
        assert len(self.watcher.connections["test_calc_001"]) == 1
        assert callback2 in self.watcher.connections["test_calc_001"]
        assert str(calc_dir) in self.watcher.watched_dirs
        
        # Remove last connection
        self.watcher.remove_connection("test_calc_001", callback2)
        
        # Verify calculation is no longer watched
        assert "test_calc_001" not in self.watcher.connections
        assert str(calc_dir) not in self.watcher.watched_dirs
    
    def test_file_change_notification(self):
        """Test file change notifications to WebSocket clients."""
        self.watcher.start()
        
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        callback_mock = Mock()
        self.watcher.add_connection("test_calc_001", callback_mock)
        
        # Simulate file change event
        file_data = {"status.json": {"status": "completed"}}
        self.watcher._on_file_change("test_calc_001", file_data)
        
        # Verify callback was called
        callback_mock.assert_called_once_with(file_data)
    
    def test_callback_error_handling(self):
        """Test error handling when callback raises exception."""
        self.watcher.start()
        
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        # Create a callback that raises an exception
        failing_callback = Mock(side_effect=Exception("Callback error"))
        working_callback = Mock()
        
        self.watcher.add_connection("test_calc_001", failing_callback)
        self.watcher.add_connection("test_calc_001", working_callback)
        
        # Simulate file change event
        file_data = {"status.json": {"status": "completed"}}
        self.watcher._on_file_change("test_calc_001", file_data)
        
        # Verify working callback was still called
        working_callback.assert_called_once_with(file_data)
        
        # Verify failing callback was removed
        assert failing_callback not in self.watcher.connections["test_calc_001"]
        assert len(self.watcher.connections["test_calc_001"]) == 1
    
    def test_get_active_connections(self):
        """Test getting active connection summary."""
        self.watcher.start()
        
        # Create calculation directories
        calc_dir1 = Path(self.temp_dir) / "test_calc_001"
        calc_dir2 = Path(self.temp_dir) / "test_calc_002"
        calc_dir1.mkdir()
        calc_dir2.mkdir()
        
        # Add connections
        callback1 = Mock()
        callback2 = Mock()
        callback3 = Mock()
        
        self.watcher.add_connection("test_calc_001", callback1)
        self.watcher.add_connection("test_calc_001", callback2)
        self.watcher.add_connection("test_calc_002", callback3)
        
        # Check active connections
        active_connections = self.watcher.get_active_connections()
        
        assert active_connections == {
            "test_calc_001": 2,
            "test_calc_002": 1
        }
    
    def test_is_watching(self):
        """Test is_watching method."""
        self.watcher.start()
        
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        # Initially not watching
        assert not self.watcher.is_watching("test_calc_001")
        
        # Add connection and verify watching
        callback_mock = Mock()
        self.watcher.add_connection("test_calc_001", callback_mock)
        assert self.watcher.is_watching("test_calc_001")
        
        # Remove connection and verify not watching
        self.watcher.remove_connection("test_calc_001", callback_mock)
        assert not self.watcher.is_watching("test_calc_001")


class TestGlobalWatcherFunctions:
    """Test cases for global watcher management functions."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        # Ensure global watcher is None before each test
        import quantum_calc.file_watcher
        quantum_calc.file_watcher._global_watcher = None
    
    def teardown_method(self):
        """Clean up test environment."""
        shutdown_websocket_watcher()
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_get_websocket_watcher_first_call(self):
        """Test first call to get_websocket_watcher."""
        watcher = get_websocket_watcher(self.temp_dir)
        
        assert isinstance(watcher, WebSocketCalculationWatcher)
        assert watcher.base_directory == Path(self.temp_dir)
        assert watcher._started
    
    def test_get_websocket_watcher_subsequent_calls(self):
        """Test subsequent calls return same instance."""
        watcher1 = get_websocket_watcher(self.temp_dir)
        watcher2 = get_websocket_watcher()  # No base_directory needed
        
        assert watcher1 is watcher2
    
    def test_get_websocket_watcher_no_base_directory_error(self):
        """Test error when base_directory not provided on first call."""
        with pytest.raises(ValueError, match="base_directory required"):
            get_websocket_watcher()
    
    def test_shutdown_websocket_watcher(self):
        """Test shutting down the global watcher."""
        # Create global watcher
        watcher = get_websocket_watcher(self.temp_dir)
        assert watcher._started
        
        # Shut it down
        shutdown_websocket_watcher()
        assert not watcher._started
        
        # Verify global instance is cleared
        import quantum_calc.file_watcher
        assert quantum_calc.file_watcher._global_watcher is None


class TestIntegrationScenarios:
    """Integration test scenarios combining multiple components."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.watcher = WebSocketCalculationWatcher(self.temp_dir)
        self.watcher.start()
    
    def teardown_method(self):
        """Clean up test environment."""
        self.watcher.stop()
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_real_file_modification_detection(self):
        """Test detection of actual file system modifications."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        callback_mock = Mock()
        self.watcher.add_connection("test_calc_001", callback_mock)
        
        # Create and modify status file
        status_file = calc_dir / "status.json"
        status_data = {"status": "running", "progress": 0}
        
        # Initial write
        status_file.write_text(json.dumps(status_data))
        time.sleep(0.1)  # Allow file system event to propagate
        
        # Update status
        status_data["progress"] = 50
        status_file.write_text(json.dumps(status_data))
        time.sleep(0.1)  # Allow file system event to propagate
        
        # Final status
        status_data["status"] = "completed"
        status_data["progress"] = 100
        status_file.write_text(json.dumps(status_data))
        time.sleep(0.1)  # Allow file system event to propagate
        
        # Verify callbacks were made (should be at least one, possibly more depending on timing)
        assert callback_mock.call_count >= 1
        
        # Verify last call contains expected data
        last_call_args = callback_mock.call_args_list[-1][0][0]
        assert "status.json" in last_call_args
        assert last_call_args["status.json"]["status"] == "completed"