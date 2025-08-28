# tests/test_concurrent_websockets.py

"""Concurrent WebSocket connection tests for file watcher implementation."""

import json
import os
import shutil
import tempfile
import threading
import time
import queue
from pathlib import Path
from unittest.mock import Mock, patch
import pytest

from quantum_calc.file_watcher import (
    WebSocketCalculationWatcher,
    get_websocket_watcher,
    shutdown_websocket_watcher
)
from quantum_calc.file_manager import CalculationFileManager


class ThreadSafeCallback:
    """Thread-safe callback for testing concurrent operations."""
    
    def __init__(self, name):
        self.name = name
        self.messages = queue.Queue()
        self.call_count = 0
        self.lock = threading.Lock()
    
    def __call__(self, data):
        with self.lock:
            self.call_count += 1
            self.messages.put((self.name, data, time.time()))
    
    def get_messages(self):
        messages = []
        while True:
            try:
                messages.append(self.messages.get_nowait())
            except queue.Empty:
                break
        return messages


class TestConcurrentWebSocketConnections:
    """Test concurrent WebSocket connections to the same calculation."""
    
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
            shutil.rmtree(self.temp_dir)
    
    def test_multiple_connections_same_calculation(self):
        """Test multiple WebSocket connections to the same calculation."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("multi_conn_test")
        calc_id = os.path.basename(calc_dir)
        
        # Create multiple callbacks
        num_connections = 5
        callbacks = [ThreadSafeCallback(f"client_{i}") for i in range(num_connections)]
        
        # Add all connections
        for callback in callbacks:
            self.watcher.add_connection(calc_id, callback)
        
        # Verify all connections are registered
        assert len(self.watcher.connections[calc_id]) == num_connections
        
        # Simulate file change
        file_data = {"status.json": {"status": "completed", "progress": 100}}
        self.watcher._on_file_change(calc_id, file_data)
        
        # Give some time for all callbacks to be processed
        time.sleep(0.1)
        
        # Verify all callbacks received the message
        for callback in callbacks:
            assert callback.call_count == 1
            messages = callback.get_messages()
            assert len(messages) == 1
            assert messages[0][1] == file_data
    
    def test_concurrent_connection_addition_removal(self):
        """Test adding and removing connections concurrently."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("concurrent_ops_test")
        calc_id = os.path.basename(calc_dir)
        
        callbacks = []
        threads = []
        
        def add_connection(callback_name):
            callback = ThreadSafeCallback(callback_name)
            callbacks.append(callback)
            self.watcher.add_connection(calc_id, callback)
            time.sleep(0.01)  # Small delay to allow interleaving
        
        def remove_connection(callback):
            time.sleep(0.05)  # Let some additions happen first
            self.watcher.remove_connection(calc_id, callback)
        
        # Start multiple threads adding connections
        for i in range(10):
            thread = threading.Thread(target=add_connection, args=(f"concurrent_client_{i}",))
            threads.append(thread)
            thread.start()
        
        # Wait for all additions to complete
        for thread in threads:
            thread.join()
        
        # Verify connections were added
        assert len(self.watcher.connections[calc_id]) == 10
        
        # Start threads removing half the connections
        remove_threads = []
        for i in range(0, 5):  # Remove first 5 callbacks
            thread = threading.Thread(target=remove_connection, args=(callbacks[i],))
            remove_threads.append(thread)
            thread.start()
        
        # Wait for removals to complete
        for thread in remove_threads:
            thread.join()
        
        # Verify remaining connections
        assert len(self.watcher.connections[calc_id]) == 5
    
    def test_high_frequency_file_changes(self):
        """Test behavior with high-frequency file changes to multiple calculations."""
        # Create multiple calculation directories
        num_calculations = 3
        calc_ids = []
        
        for i in range(num_calculations):
            calc_dir = self.file_manager.create_calculation_dir(f"high_freq_test_{i}")
            calc_ids.append(os.path.basename(calc_dir))
        
        # Add callbacks for each calculation
        all_callbacks = {}
        for calc_id in calc_ids:
            callbacks = [ThreadSafeCallback(f"{calc_id}_client_{i}") for i in range(3)]
            all_callbacks[calc_id] = callbacks
            for callback in callbacks:
                self.watcher.add_connection(calc_id, callback)
        
        # Function to generate file changes
        def generate_file_changes(calc_id, num_changes):
            for i in range(num_changes):
                file_data = {
                    "status.json": {
                        "status": "running",
                        "progress": i * 10,
                        "iteration": i
                    }
                }
                self.watcher._on_file_change(calc_id, file_data)
                time.sleep(0.01)  # Small delay between changes
        
        # Start threads generating changes for each calculation
        change_threads = []
        for calc_id in calc_ids:
            thread = threading.Thread(target=generate_file_changes, args=(calc_id, 20))
            change_threads.append(thread)
            thread.start()
        
        # Wait for all changes to complete
        for thread in change_threads:
            thread.join()
        
        # Give time for all callbacks to be processed
        time.sleep(0.2)
        
        # Verify all callbacks received messages
        for calc_id, callbacks in all_callbacks.items():
            for callback in callbacks:
                assert callback.call_count > 0, f"Callback {callback.name} didn't receive any messages"
                messages = callback.get_messages()
                assert len(messages) > 0
    
    def test_connection_isolation(self):
        """Test that connections to different calculations are isolated."""
        # Create two calculation directories
        calc_dir1 = self.file_manager.create_calculation_dir("isolation_test_1")
        calc_dir2 = self.file_manager.create_calculation_dir("isolation_test_2")
        calc_id1 = os.path.basename(calc_dir1)
        calc_id2 = os.path.basename(calc_dir2)
        
        # Create callbacks for each calculation
        callback1 = ThreadSafeCallback("calc1_client")
        callback2 = ThreadSafeCallback("calc2_client")
        
        self.watcher.add_connection(calc_id1, callback1)
        self.watcher.add_connection(calc_id2, callback2)
        
        # Send change to calculation 1
        file_data1 = {"status.json": {"status": "completed", "calc": "calc1"}}
        self.watcher._on_file_change(calc_id1, file_data1)
        
        time.sleep(0.1)
        
        # Send change to calculation 2
        file_data2 = {"status.json": {"status": "error", "calc": "calc2"}}
        self.watcher._on_file_change(calc_id2, file_data2)
        
        time.sleep(0.1)
        
        # Verify each callback only received its respective data
        messages1 = callback1.get_messages()
        messages2 = callback2.get_messages()
        
        assert len(messages1) == 1
        assert len(messages2) == 1
        assert messages1[0][1]["status.json"]["calc"] == "calc1"
        assert messages2[0][1]["status.json"]["calc"] == "calc2"
    
    def test_callback_exception_isolation(self):
        """Test that exceptions in one callback don't affect others."""
        # Create calculation directory
        calc_dir = self.file_manager.create_calculation_dir("exception_isolation_test")
        calc_id = os.path.basename(calc_dir)
        
        # Create callbacks - one that fails, others that work
        failing_callback = Mock(side_effect=Exception("Callback failed"))
        working_callback1 = ThreadSafeCallback("working_client_1")
        working_callback2 = ThreadSafeCallback("working_client_2")
        
        # Add all callbacks
        self.watcher.add_connection(calc_id, failing_callback)
        self.watcher.add_connection(calc_id, working_callback1)
        self.watcher.add_connection(calc_id, working_callback2)
        
        # Simulate file change
        file_data = {"status.json": {"status": "completed"}}
        self.watcher._on_file_change(calc_id, file_data)
        
        time.sleep(0.1)
        
        # Verify working callbacks received the message
        assert working_callback1.call_count == 1
        assert working_callback2.call_count == 1
        
        # Verify failing callback was removed
        assert failing_callback not in self.watcher.connections[calc_id]
        assert len(self.watcher.connections[calc_id]) == 2


class TestScalabilityScenarios:
    """Test scalability scenarios with many concurrent connections."""
    
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
    
    def test_many_concurrent_calculations(self):
        """Test handling many concurrent calculations."""
        # Create many calculation directories
        num_calculations = 20
        calc_ids = []
        
        file_manager = CalculationFileManager(self.temp_dir)
        
        for i in range(num_calculations):
            calc_dir = file_manager.create_calculation_dir(f"scale_test_{i:03d}")
            calc_ids.append(os.path.basename(calc_dir))
        
        # Add one callback per calculation
        callbacks = {}
        for calc_id in calc_ids:
            callback = ThreadSafeCallback(f"{calc_id}_client")
            callbacks[calc_id] = callback
            self.watcher.add_connection(calc_id, callback)
        
        # Verify all connections are registered
        assert len(self.watcher.connections) == num_calculations
        
        # Send changes to all calculations concurrently
        def send_change(calc_id):
            file_data = {"status.json": {"status": "completed", "calc_id": calc_id}}
            self.watcher._on_file_change(calc_id, file_data)
        
        threads = []
        for calc_id in calc_ids:
            thread = threading.Thread(target=send_change, args=(calc_id,))
            threads.append(thread)
            thread.start()
        
        # Wait for all changes to be sent
        for thread in threads:
            thread.join()
        
        time.sleep(0.2)  # Allow all callbacks to be processed
        
        # Verify all callbacks received their respective messages
        for calc_id, callback in callbacks.items():
            assert callback.call_count == 1
            messages = callback.get_messages()
            assert len(messages) == 1
            assert messages[0][1]["status.json"]["calc_id"] == calc_id
    
    def test_many_connections_per_calculation(self):
        """Test many connections to a single calculation."""
        # Create one calculation directory
        file_manager = CalculationFileManager(self.temp_dir)
        calc_dir = file_manager.create_calculation_dir("many_conn_test")
        calc_id = os.path.basename(calc_dir)
        
        # Add many callbacks to the same calculation
        num_connections = 50
        callbacks = []
        
        for i in range(num_connections):
            callback = ThreadSafeCallback(f"client_{i:03d}")
            callbacks.append(callback)
            self.watcher.add_connection(calc_id, callback)
        
        # Verify all connections are registered
        assert len(self.watcher.connections[calc_id]) == num_connections
        
        # Send a single file change
        file_data = {"status.json": {"status": "completed", "timestamp": time.time()}}
        self.watcher._on_file_change(calc_id, file_data)
        
        time.sleep(0.5)  # Allow all callbacks to be processed
        
        # Verify all callbacks received the message
        total_calls = sum(callback.call_count for callback in callbacks)
        assert total_calls == num_connections
        
        for callback in callbacks:
            assert callback.call_count == 1
            messages = callback.get_messages()
            assert len(messages) == 1
            assert messages[0][1] == file_data
    
    def test_rapid_connection_churn(self):
        """Test rapid addition and removal of connections."""
        # Create calculation directory
        file_manager = CalculationFileManager(self.temp_dir)
        calc_dir = file_manager.create_calculation_dir("churn_test")
        calc_id = os.path.basename(calc_dir)
        
        # Function to add and remove connections rapidly
        def connection_churn():
            for i in range(10):
                callback = ThreadSafeCallback(f"churn_client_{threading.current_thread().ident}_{i}")
                
                # Add connection
                self.watcher.add_connection(calc_id, callback)
                time.sleep(0.01)
                
                # Remove connection
                self.watcher.remove_connection(calc_id, callback)
                time.sleep(0.01)
        
        # Start multiple threads doing connection churn
        threads = []
        for i in range(5):
            thread = threading.Thread(target=connection_churn)
            threads.append(thread)
            thread.start()
        
        # Wait for all churn to complete
        for thread in threads:
            thread.join()
        
        # Verify no connections remain (or very few due to timing)
        remaining_connections = self.watcher.connections.get(calc_id, set())
        assert len(remaining_connections) == 0, f"Expected no connections, but found {len(remaining_connections)}"


class TestResourceManagement:
    """Test proper resource management in concurrent scenarios."""
    
    def setup_method(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
    
    def teardown_method(self):
        """Clean up test environment."""
        shutdown_websocket_watcher()
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_global_watcher_thread_safety(self):
        """Test thread safety of global watcher access."""
        watchers = []
        exceptions = []
        
        def get_global_watcher():
            try:
                watcher = get_websocket_watcher(self.temp_dir)
                watchers.append(watcher)
            except Exception as e:
                exceptions.append(e)
        
        # Start multiple threads accessing global watcher
        threads = []
        for i in range(10):
            thread = threading.Thread(target=get_global_watcher)
            threads.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join()
        
        # Verify no exceptions occurred
        assert len(exceptions) == 0, f"Exceptions occurred: {exceptions}"
        
        # Verify all threads got the same watcher instance
        assert len(set(id(w) for w in watchers)) == 1, "All threads should get the same watcher instance"
        assert len(watchers) == 10, "All threads should successfully get a watcher"
    
    def test_memory_usage_with_many_connections(self):
        """Test that memory usage remains reasonable with many connections."""
        import gc
        
        watcher = get_websocket_watcher(self.temp_dir)
        
        # Create many calculations with multiple connections each
        file_manager = CalculationFileManager(self.temp_dir)
        callbacks = []
        
        for calc_num in range(10):
            calc_dir = file_manager.create_calculation_dir(f"memory_test_{calc_num}")
            calc_id = os.path.basename(calc_dir)
            
            for conn_num in range(20):
                callback = ThreadSafeCallback(f"calc_{calc_num}_conn_{conn_num}")
                callbacks.append(callback)
                watcher.add_connection(calc_id, callback)
        
        # Force garbage collection to get accurate memory baseline
        gc.collect()
        initial_objects = len(gc.get_objects())
        
        # Remove all connections
        for calc_num in range(10):
            calc_id = f"memory_test_{calc_num}_" + "20241201_120000"  # Approximate calc_id format
            actual_calc_ids = [k for k in watcher.connections.keys() if k.startswith(f"memory_test_{calc_num}")]
            
            for actual_calc_id in actual_calc_ids:
                # Remove all callbacks for this calculation
                calc_callbacks = list(watcher.connections[actual_calc_id])
                for callback in calc_callbacks:
                    watcher.remove_connection(actual_calc_id, callback)
        
        # Force garbage collection
        gc.collect()
        final_objects = len(gc.get_objects())
        
        # Verify connections were cleaned up
        assert len(watcher.connections) == 0
        
        # Memory usage should not have grown significantly
        # (Allow some variance due to test framework overhead)
        object_growth = final_objects - initial_objects
        assert object_growth < 1000, f"Too many objects leaked: {object_growth}"