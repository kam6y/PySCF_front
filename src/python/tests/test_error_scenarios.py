# tests/test_error_scenarios.py

"""Error scenario tests for file watcher and WebSocket implementation."""

import json
import os
import shutil
import tempfile
import threading
import time
import stat
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pytest

from quantum_calc.file_watcher import (
    CalculationFileWatcher,
    WebSocketCalculationWatcher,
    get_websocket_watcher,
    shutdown_websocket_watcher
)
from quantum_calc.file_manager import CalculationFileManager
from quantum_calc.exceptions import FileManagerError, WebSocketError


class TestFileSystemErrorScenarios:
    """Test error scenarios related to file system operations."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.file_manager = CalculationFileManager(self.temp_dir)
        self.watcher = WebSocketCalculationWatcher(self.temp_dir)
        self.watcher.start()
    
    def teardown_method(self):
        """Clean up test environment."""
        self.watcher.stop()
        if os.path.exists(self.temp_dir):
            # Remove read-only permissions before cleanup
            for root, dirs, files in os.walk(self.temp_dir):
                for d in dirs:
                    os.chmod(os.path.join(root, d), stat.S_IRWXU)
                for f in files:
                    os.chmod(os.path.join(root, f), stat.S_IRWXU)
            shutil.rmtree(self.temp_dir)
    
    def test_file_deletion_during_monitoring(self):
        """Test behavior when monitored files are deleted."""
        # Create calculation directory and files
        calc_dir = Path(self.temp_dir) / "test_calc_001"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_data = {"status": "running", "progress": 50}
        status_file.write_text(json.dumps(status_data))
        
        callback_mock = Mock()
        
        # Set up file watcher
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification then deletion
        from watchdog.events import FileModifiedEvent, FileDeletedEvent
        
        # First modification should work
        mod_event = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event)
        callback_mock.assert_called_once()
        
        # Delete the file
        status_file.unlink()
        
        # Subsequent modification event should be handled gracefully
        callback_mock.reset_mock()
        mod_event2 = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event2)
        
        # Should not call callback for deleted file
        callback_mock.assert_not_called()
    
    def test_permission_denied_error(self):
        """Test behavior when file permissions prevent access."""
        # Create calculation directory and files
        calc_dir = Path(self.temp_dir) / "test_calc_002"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_data = {"status": "running", "progress": 50}
        status_file.write_text(json.dumps(status_data))
        
        # Remove read permissions
        os.chmod(str(status_file), 0o000)
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event)
        
        # Should handle permission error gracefully
        callback_mock.assert_not_called()
        
        # Restore permissions for cleanup
        os.chmod(str(status_file), stat.S_IRWXU)
    
    def test_corrupted_json_file(self):
        """Test behavior when JSON files are corrupted."""
        # Create calculation directory with corrupted JSON
        calc_dir = Path(self.temp_dir) / "test_calc_003"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_file.write_text("{ invalid json content")
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event)
        
        # Should handle JSON decode error gracefully
        callback_mock.assert_not_called()
    
    def test_directory_deletion_during_monitoring(self):
        """Test behavior when entire calculation directory is deleted."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_004"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_data = {"status": "running"}
        status_file.write_text(json.dumps(status_data))
        
        callback_mock = Mock()
        self.watcher.add_connection("test_calc_004", callback_mock)
        
        # Verify directory is being watched
        assert self.watcher.is_watching("test_calc_004")
        
        # Delete entire directory
        shutil.rmtree(str(calc_dir))
        
        # Simulate file events after deletion - should be handled gracefully
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(status_file))
        self.watcher.event_handler.on_modified(mod_event)
        
        # Should not crash or call callback
        callback_mock.assert_not_called()
    
    def test_disk_full_scenario(self):
        """Test behavior when disk is full (simulated by mocking)."""
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Mock file operations to raise OSError (disk full)
        with patch('builtins.open', side_effect=OSError("No space left on device")):
            # Create calculation directory
            calc_dir = Path(self.temp_dir) / "test_calc_005"
            calc_dir.mkdir()
            
            status_file = calc_dir / "status.json"
            
            # Simulate file modification
            from watchdog.events import FileModifiedEvent
            mod_event = FileModifiedEvent(str(status_file))
            file_watcher.on_modified(mod_event)
            
            # Should handle OSError gracefully
            callback_mock.assert_not_called()
    
    def test_file_locked_by_another_process(self):
        """Test behavior when file is locked by another process."""
        # Create calculation directory and file
        calc_dir = Path(self.temp_dir) / "test_calc_006"
        calc_dir.mkdir()
        
        status_file = calc_dir / "status.json"
        status_data = {"status": "running"}
        status_file.write_text(json.dumps(status_data))
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Mock file opening to raise PermissionError
        with patch('builtins.open', side_effect=PermissionError("File is locked")):
            from watchdog.events import FileModifiedEvent
            mod_event = FileModifiedEvent(str(status_file))
            file_watcher.on_modified(mod_event)
            
            # Should handle permission error gracefully
            callback_mock.assert_not_called()


class TestWebSocketErrorScenarios:
    """Test error scenarios specific to WebSocket operations."""
    
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
    
    def test_websocket_send_failure(self):
        """Test handling of WebSocket send failures."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_007"
        calc_dir.mkdir()
        
        # Create callback that raises exception
        failing_callback = Mock(side_effect=Exception("WebSocket send failed"))
        working_callback = Mock()
        
        # Add both callbacks
        self.watcher.add_connection("test_calc_007", failing_callback)
        self.watcher.add_connection("test_calc_007", working_callback)
        
        # Simulate file change
        file_data = {"status.json": {"status": "completed"}}
        self.watcher._on_file_change("test_calc_007", file_data)
        
        # Working callback should still be called
        working_callback.assert_called_once_with(file_data)
        
        # Failing callback should be removed
        assert failing_callback not in self.watcher.connections["test_calc_007"]
        assert len(self.watcher.connections["test_calc_007"]) == 1
    
    def test_multiple_callback_failures(self):
        """Test behavior when multiple callbacks fail."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_008"
        calc_dir.mkdir()
        
        # Create multiple failing callbacks
        failing_callback1 = Mock(side_effect=Exception("Callback 1 failed"))
        failing_callback2 = Mock(side_effect=Exception("Callback 2 failed"))
        
        # Add failing callbacks
        self.watcher.add_connection("test_calc_008", failing_callback1)
        self.watcher.add_connection("test_calc_008", failing_callback2)
        
        # Simulate file change
        file_data = {"status.json": {"status": "error"}}
        self.watcher._on_file_change("test_calc_008", file_data)
        
        # Both callbacks should be removed
        assert "test_calc_008" not in self.watcher.connections
    
    def test_observer_start_failure(self):
        """Test behavior when file observer fails to start."""
        watcher = WebSocketCalculationWatcher(self.temp_dir)
        
        # Mock observer start to raise exception
        with patch.object(watcher.observer, 'start', side_effect=Exception("Observer start failed")):
            # Should handle exception gracefully
            try:
                watcher.start()
            except Exception:
                pytest.fail("Watcher.start() should handle observer start failure gracefully")
        
        watcher.stop()
    
    def test_observer_stop_failure(self):
        """Test behavior when file observer fails to stop."""
        watcher = WebSocketCalculationWatcher(self.temp_dir)
        watcher.start()
        
        # Mock observer stop to raise exception
        with patch.object(watcher.observer, 'stop', side_effect=Exception("Observer stop failed")):
            # Should handle exception gracefully
            try:
                watcher.stop()
            except Exception:
                pytest.fail("Watcher.stop() should handle observer stop failure gracefully")
    
    def test_connection_removal_with_nonexistent_callback(self):
        """Test removing a callback that doesn't exist."""
        # Create calculation directory
        calc_dir = Path(self.temp_dir) / "test_calc_009"
        calc_dir.mkdir()
        
        callback_mock = Mock()
        nonexistent_callback = Mock()
        
        # Add one callback
        self.watcher.add_connection("test_calc_009", callback_mock)
        
        # Try to remove callback that was never added
        # Should not raise exception
        try:
            self.watcher.remove_connection("test_calc_009", nonexistent_callback)
        except Exception:
            pytest.fail("remove_connection should handle nonexistent callback gracefully")
        
        # Original callback should still be present
        assert callback_mock in self.watcher.connections["test_calc_009"]


class TestGlobalWatcherErrorScenarios:
    """Test error scenarios for global watcher management."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        shutdown_websocket_watcher()
        import quantum_calc.file_watcher
        quantum_calc.file_watcher._global_watcher = None
    
    def teardown_method(self):
        """Clean up test environment."""
        shutdown_websocket_watcher()
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_get_watcher_with_invalid_directory(self):
        """Test get_websocket_watcher with invalid base directory."""
        # Try to create watcher with non-existent directory
        nonexistent_dir = "/nonexistent/directory"
        
        # Should handle invalid directory gracefully
        try:
            watcher = get_websocket_watcher(nonexistent_dir)
            # If it doesn't raise an exception, ensure it's properly cleaned up
            watcher.stop()
        except Exception as e:
            # Some exceptions are expected for invalid directories
            assert "nonexistent" in str(e).lower() or "permission" in str(e).lower()
    
    def test_shutdown_already_shutdown_watcher(self):
        """Test shutting down an already shutdown watcher."""
        # Create and shutdown watcher
        watcher = get_websocket_watcher(self.temp_dir)
        shutdown_websocket_watcher()
        
        # Try to shutdown again - should not raise exception
        try:
            shutdown_websocket_watcher()
        except Exception:
            pytest.fail("shutdown_websocket_watcher should handle already shutdown watcher gracefully")
    
    def test_concurrent_watcher_operations(self):
        """Test concurrent operations on global watcher."""
        import threading
        
        results = []
        exceptions = []
        
        def create_watcher():
            try:
                watcher = get_websocket_watcher(self.temp_dir)
                results.append(watcher)
            except Exception as e:
                exceptions.append(e)
        
        def shutdown_watcher():
            try:
                shutdown_websocket_watcher()
            except Exception as e:
                exceptions.append(e)
        
        # Create multiple threads that try to access global watcher
        threads = []
        for i in range(5):
            if i < 3:
                thread = threading.Thread(target=create_watcher)
            else:
                thread = threading.Thread(target=shutdown_watcher)
            threads.append(thread)
        
        # Start all threads
        for thread in threads:
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join(timeout=1.0)
        
        # Should not have any unexpected exceptions
        assert len(exceptions) == 0, f"Unexpected exceptions: {exceptions}"


class TestFileSystemIntegrityScenarios:
    """Test scenarios involving file system integrity issues."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.file_manager = CalculationFileManager(self.temp_dir)
    
    def teardown_method(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_partial_file_write(self):
        """Test behavior with partially written files."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("partial_write_test")
        calc_id = os.path.basename(calc_dir)
        
        status_file = Path(calc_dir) / "status.json"
        
        # Write partial JSON (missing closing brace)
        status_file.write_text('{"status": "running", "progress": 50')
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event)
        
        # Should handle partial JSON gracefully
        callback_mock.assert_not_called()
    
    def test_file_encoding_issues(self):
        """Test behavior with files having encoding issues."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("encoding_test")
        
        status_file = Path(calc_dir) / "status.json"
        
        # Write file with invalid UTF-8 encoding
        with open(status_file, 'wb') as f:
            f.write(b'{"status": "running", "message": "\xff\xfe invalid utf-8"}')
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event)
        
        # Should handle encoding errors gracefully
        callback_mock.assert_not_called()
    
    def test_extremely_large_json_file(self):
        """Test behavior with extremely large JSON files."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("large_file_test")
        
        status_file = Path(calc_dir) / "status.json"
        
        # Create large JSON content (but not too large to cause memory issues in tests)
        large_data = {
            "status": "running",
            "large_array": [i for i in range(10000)],  # 10k integers
            "large_string": "x" * 100000  # 100k characters
        }
        
        status_file.write_text(json.dumps(large_data))
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(status_file))
        file_watcher.on_modified(mod_event)
        
        # Should handle large files
        callback_mock.assert_called_once()
        call_args = callback_mock.call_args[0]
        assert call_args[0] == os.path.basename(calc_dir)
        assert "status.json" in call_args[1]
        assert call_args[1]["status.json"]["status"] == "running"
    
    def test_symlink_handling(self):
        """Test behavior with symbolic links."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("symlink_test")
        
        # Create actual status file
        actual_status_file = Path(calc_dir) / "actual_status.json"
        status_data = {"status": "running", "progress": 75}
        actual_status_file.write_text(json.dumps(status_data))
        
        # Create symlink
        symlink_status_file = Path(calc_dir) / "status.json"
        try:
            os.symlink(str(actual_status_file), str(symlink_status_file))
        except OSError:
            # Skip test if symlinks not supported on this platform
            pytest.skip("Symbolic links not supported on this platform")
        
        callback_mock = Mock()
        file_watcher = CalculationFileWatcher(callback_mock)
        
        # Simulate file modification on symlink
        from watchdog.events import FileModifiedEvent
        mod_event = FileModifiedEvent(str(symlink_status_file))
        file_watcher.on_modified(mod_event)
        
        # Should handle symlinks properly
        callback_mock.assert_called_once()
        call_args = callback_mock.call_args[0]
        assert call_args[1]["status.json"]["status"] == "running"