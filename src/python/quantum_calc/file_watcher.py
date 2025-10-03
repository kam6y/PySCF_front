# src/python/quantum_calc/file_watcher.py

"""File system monitoring utilities for quantum chemistry calculations."""

import os
import sys
import json
import logging
import threading
from pathlib import Path
from typing import Dict, Callable, Optional, Set, Any
from datetime import datetime
import queue
import threading

# Platform-specific observer selection
# FSEventsObserver (macOS default) fails after fork in Gunicorn workers
# Use PollingObserver on macOS for fork-safety, regular Observer elsewhere
if sys.platform == 'darwin':
    from watchdog.observers.polling import PollingObserver as Observer
    logging.getLogger(__name__).info("Using PollingObserver for macOS fork-safety")
else:
    from watchdog.observers import Observer

from watchdog.events import FileSystemEventHandler, FileModifiedEvent, FileCreatedEvent

from .exceptions import FileManagerError, WebSocketError


logger = logging.getLogger(__name__)


class CalculationFileWatcher(FileSystemEventHandler):
    """
    File system event handler for calculation directories.
    Monitors status.json, results.json, and parameters.json files.
    """
    
    def __init__(self, callback: Callable[[str, Dict], None]):
        """
        Initialize the file watcher.
        
        Args:
            callback: Function to call when a monitored file changes.
                      Takes calculation_id and file data as arguments.
        """
        super().__init__()
        self.callback = callback
        self.monitored_files = {'status.json', 'results.json', 'parameters.json'}
        self._lock = threading.Lock()
    
    def on_modified(self, event):
        """Handle file modification events."""
        if not event.is_directory:
            self._handle_file_event(event.src_path)
    
    def on_created(self, event):
        """Handle file creation events.""" 
        if not event.is_directory:
            self._handle_file_event(event.src_path)
    
    def _handle_file_event(self, file_path: str):
        """Process file system events for monitored files."""
        try:
            file_name = os.path.basename(file_path)
            if file_name not in self.monitored_files:
                return
                
            # Extract calculation_id from the directory path
            calc_dir = os.path.dirname(file_path)
            calculation_id = os.path.basename(calc_dir)
            
            logger.debug(f"File event detected: {file_name} in calculation {calculation_id}")
            
            # Read the file content
            file_data = self._read_file_safely(file_path, file_name)
            if file_data is not None:
                with self._lock:
                    self.callback(calculation_id, {file_name: file_data})
                    
        except Exception as e:
            logger.error(f"Error handling file event for {file_path}: {e}")
    
    def _read_file_safely(self, file_path: str, file_name: str) -> Optional[Dict]:
        """Safely read and parse JSON files."""
        try:
            if not os.path.exists(file_path):
                return None
                
            with open(file_path, 'r') as f:
                if file_name.endswith('.json'):
                    return json.load(f)
                else:
                    return {'content': f.read()}
                    
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Failed to read {file_name}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error reading {file_name}: {e}")
            return None


class WebSocketCalculationWatcher:
    """
    Manages file system watching for WebSocket connections to calculation directories.
    Provides efficient, event-driven monitoring instead of polling.
    Now with threading integration for thread-safe operation.
    """
    
    def __init__(self, base_directory: str):
        """
        Initialize the WebSocket calculation watcher.
        
        Args:
            base_directory: Base directory containing calculation folders
        """
        self.base_directory = Path(base_directory)
        self.observer = Observer()
        self.connections: Dict[str, Set[Callable]] = {}  # calculation_id -> set of callbacks
        self.watched_dirs: Dict[str, Any] = {}  # calculation_id -> ObservedWatch object from schedule()
        self._lock = threading.Lock()
        self._started = False
        
        # Standard queue for thread-safe communication
        self._event_queue = queue.Queue()
        self._queue_processor = None
        self._stop_processing = threading.Event()
        
        # Create event handler that puts events in the queue
        self.event_handler = CalculationFileWatcher(self._queue_file_change)
    
    def start(self):
        """Start the file system observer and event queue processor."""
        with self._lock:
            if not self._started:
                self.observer.start()
                # Start the threading queue processor
                self._queue_processor = threading.Thread(target=self._process_file_events, daemon=True)
                self._queue_processor.start()
                self._started = True
                logger.info("WebSocket file watcher started with threading integration")
    
    def stop(self):
        """Stop the file system observer and clean up resources."""
        with self._lock:
            if self._started:
                try:
                    # Stop the observer
                    self.observer.stop()
                    self.observer.join(timeout=5.0)
                    
                    # Stop the queue processor
                    if self._queue_processor:
                        self._stop_processing.set()
                        self._queue_processor.join(timeout=2.0)
                        self._queue_processor = None
                        
                except Exception as e:
                    logger.warning(f"Error stopping file observer: {e}")
                finally:
                    self._started = False
                    self.connections.clear()
                    self.watched_dirs.clear()
                    logger.info("WebSocket file watcher stopped")
    
    def add_connection(self, calculation_id: str, callback: Callable[[Dict], None]):
        """
        Add a WebSocket connection to monitor a calculation.
        
        Args:
            calculation_id: Unique calculation identifier
            callback: Function to call when calculation data changes
        """
        with self._lock:
            if calculation_id not in self.connections:
                self.connections[calculation_id] = set()
            
            self.connections[calculation_id].add(callback)
            
            # Add directory to watch list if not already watched
            calc_dir = self.base_directory / calculation_id
            calc_dir_str = str(calc_dir)
            
            if calc_dir_str not in self.watched_dirs and calc_dir.exists():
                # Store the ObservedWatch object returned by schedule() for later unscheduling
                watch = self.observer.schedule(self.event_handler, calc_dir_str, recursive=False)
                self.watched_dirs[calc_dir_str] = watch
                logger.debug(f"Started watching directory: {calc_dir_str}")
    
    def remove_connection(self, calculation_id: str, callback: Callable[[Dict], None]):
        """
        Remove a WebSocket connection from monitoring a calculation.
        
        Args:
            calculation_id: Unique calculation identifier
            callback: The callback function to remove
        """
        with self._lock:
            if calculation_id in self.connections:
                self.connections[calculation_id].discard(callback)
                
                # If no more connections for this calculation, stop watching
                if not self.connections[calculation_id]:
                    del self.connections[calculation_id]
                    self._stop_watching_calculation(calculation_id)
    
    def _stop_watching_calculation(self, calculation_id: str):
        """
        Stop watching a calculation directory when no connections remain.
        
        Uses the official watchdog API pattern: store the ObservedWatch object
        returned by schedule() and pass it to unschedule() later.
        This avoids relying on internal implementation details of the watchdog library.
        """
        calc_dir = self.base_directory / calculation_id
        calc_dir_str = str(calc_dir)
        
        if calc_dir_str in self.watched_dirs:
            try:
                # Use the stored ObservedWatch object for clean unscheduling (official pattern)
                watch = self.watched_dirs[calc_dir_str]
                self.observer.unschedule(watch)
                del self.watched_dirs[calc_dir_str]
                logger.debug(f"Stopped watching directory: {calc_dir_str}")
            except Exception as e:
                logger.error(f"Error unscheduling watch for {calc_dir_str}: {e}")
                # Clean up the entry even if unscheduling failed
                self.watched_dirs.pop(calc_dir_str, None)
    
    def _queue_file_change(self, calculation_id: str, file_data: Dict):
        """Queue file change events for processing in threading context."""
        # This runs in watchdog thread - just put event in queue
        try:
            self._event_queue.put((calculation_id, file_data), block=False)
        except Exception as e:
            logger.error(f"Error queuing file change event for {calculation_id}: {e}")
    
    def _process_file_events(self):
        """Process file events from the queue in threading context."""
        logger.debug("Started threading file event processor")
        try:
            while not self._stop_processing.is_set():
                try:
                    # Get event from queue with timeout to allow checking stop event
                    try:
                        calculation_id, file_data = self._event_queue.get(timeout=1.0)
                    except queue.Empty:
                        continue  # Check stop event and continue
                    
                    # Process the file change in threading context
                    self._handle_file_change_in_threading(calculation_id, file_data)
                    
                    # Mark task as done
                    self._event_queue.task_done()
                    
                except Exception as e:
                    logger.error(f"Error processing file event: {e}")
                    # Continue processing other events
                    continue
        except Exception as e:
            logger.error(f"Unexpected error in file event processor: {e}")
        finally:
            logger.debug("File event processor shutting down")
    
    def _handle_file_change_in_threading(self, calculation_id: str, file_data: Dict):
        """Handle file change events and notify connected WebSocket clients (threading context)."""
        with self._lock:
            if calculation_id not in self.connections:
                return
            
            # Get current connections for this calculation
            callbacks = self.connections[calculation_id].copy()
        
        # Notify all connected WebSocket clients (now in threading context)
        for callback in callbacks:
            try:
                callback(file_data)
            except Exception as e:
                logger.error(f"Error notifying WebSocket client for {calculation_id}: {e}")
                # Remove the failed callback to prevent future errors
                with self._lock:
                    if calculation_id in self.connections:
                        self.connections[calculation_id].discard(callback)
    
    def get_active_connections(self) -> Dict[str, int]:
        """Get a summary of active connections by calculation ID."""
        with self._lock:
            return {calc_id: len(callbacks) for calc_id, callbacks in self.connections.items()}
    
    def is_watching(self, calculation_id: str) -> bool:
        """Check if a calculation is currently being watched."""
        calc_dir = self.base_directory / calculation_id
        calc_dir_str = str(calc_dir)
        return calc_dir_str in self.watched_dirs


# Global instance for the application
_global_watcher: Optional[WebSocketCalculationWatcher] = None
_watcher_lock = threading.Lock()


def get_websocket_watcher(base_directory: Optional[str] = None) -> WebSocketCalculationWatcher:
    """
    Get or create the global WebSocket file watcher instance.
    
    Args:
        base_directory: Base directory for calculations (required on first call)
    
    Returns:
        WebSocketCalculationWatcher instance
    """
    global _global_watcher
    
    with _watcher_lock:
        if _global_watcher is None:
            if base_directory is None:
                raise ValueError("base_directory required for first call to get_websocket_watcher")
            _global_watcher = WebSocketCalculationWatcher(base_directory)
            _global_watcher.start()
        
        return _global_watcher


def shutdown_websocket_watcher():
    """Shutdown the global WebSocket file watcher."""
    global _global_watcher
    
    with _watcher_lock:
        if _global_watcher is not None:
            _global_watcher.stop()
            _global_watcher = None