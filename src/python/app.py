import logging
import os
import sys
import json
import signal
import atexit
import socket
from flask import Flask, jsonify
from flask_cors import CORS
from flask_socketio import SocketIO
from pydantic import ValidationError

# Import blueprint registration and WebSocket handlers
from api import register_blueprints
from websocket import register_websocket_handlers

# Import cleanup functions from quantum_calc
from quantum_calc import shutdown_process_manager, shutdown_websocket_watcher
from quantum_calc.exceptions import ProcessManagerError


def load_server_config():
    """Load server configuration from JSON file."""
    try:
        # Load from config directory (single source of truth)
        config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config', 'server-config.json')
        
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Server configuration file not found at {config_path}")
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = json.load(f)
        
        print(f"Loaded server configuration from: {config_path}")
        return config
    except Exception as e:
        print(f"WARNING: Failed to load server configuration: {e}. Using defaults.")
        # Return default configuration that matches config/server-config.json structure
        return {
            "server": {
                "host": "127.0.0.1",
                "port": {
                    "default": 5000,
                    "auto_detect": True,
                    "range": {"start": 5000, "end": 5100}
                },
                "debug": False
            },
            "gunicorn": {
                "workers": 1,
                "threads": 4,
                "worker_class": "sync",
                "timeout": 0,
                "keep_alive": 30,
                "preload_app": True,
                "access_logfile": "-",
                "log_level": "info"
            },
            "socketio": {
                "cors_allowed_origins": ["http://127.0.0.1:*", "ws://127.0.0.1:*", "file://"],
                "async_mode": "threading",
                "ping_timeout": 60,
                "ping_interval": 25,
                "allow_unsafe_werkzeug": True,
                "logger": True,
                "engineio_logger": False
            },
            "development": {"use_reloader": False, "debug": False, "enable_dev_tools": True},
            "production": {"use_gunicorn": True, "optimize_performance": True, "enable_logging": True},
            "quantum_calculations": {
                "max_concurrent_calculations": "auto",
                "process_pool_size": "auto",
                "timeout_calculation": 0,
                "memory_limit_mb": 0
            },
            "logging": {
                "level": "INFO",
                "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                "enable_access_log": True,
                "enable_error_log": True
            }
        }


# Global configuration
SERVER_CONFIG = load_server_config()

# Configure logging based on configuration
log_level = getattr(logging, SERVER_CONFIG.get('logging', {}).get('level', 'INFO').upper())
log_format = SERVER_CONFIG.get('logging', {}).get('format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logging.basicConfig(
    level=log_level,
    format=log_format
)
logger = logging.getLogger(__name__)


def create_app():
    """Application factory for Gunicorn compatibility."""
    # Initialize Flask app
    app = Flask(__name__)
    CORS(app)  # Enable CORS for cross-origin requests

    # Initialize SocketIO with configuration-based settings
    socketio_config = SERVER_CONFIG.get('socketio', {})
    socketio = SocketIO(
        app, 
        cors_allowed_origins=socketio_config.get('cors_allowed_origins', ["http://127.0.0.1:*", "ws://127.0.0.1:*", "file://"]),
        async_mode=socketio_config.get('async_mode', 'threading'),
        ping_timeout=socketio_config.get('ping_timeout', 60),
        ping_interval=socketio_config.get('ping_interval', 25),
        allow_unsafe_werkzeug=socketio_config.get('allow_unsafe_werkzeug', True),
        logger=socketio_config.get('logger', True),
        engineio_logger=socketio_config.get('engineio_logger', False)
    )

    # Register all API blueprints
    register_blueprints(app)
    logger.info("Registered all API blueprints")

    # Register WebSocket handlers and get immediate notification function
    websocket_funcs = register_websocket_handlers(socketio)
    logger.info("Registered all WebSocket handlers")
    
    # Store the immediate notification function globally for other modules to use
    app.send_immediate_websocket_notification = websocket_funcs['send_immediate_websocket_notification']

    # Error handlers
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

    @app.errorhandler(400)
    def bad_request_handler(error):
        """Handle bad requests, filtering out WebSocket close frame errors."""
        import re
        
        # Check if this is a WebSocket close frame being misinterpreted as HTTP
        error_description = str(error).lower()
        
        # Expanded list of WebSocket protocol indicators
        websocket_indicators = [
            'invalid http method', 'expected get method', 'invalid method', 
            'websocket', 'connection upgrade', 'upgrade required', 
            'bad request line', 'malformed request', 'protocol error',
            'connection reset', 'connection closed', 'connection aborted',
            'invalid utf-8', 'decode error', 'unicode error'
        ]
        
        # Check simple string indicators
        is_websocket_related = any(indicator in error_description for indicator in websocket_indicators)
        
        # Additional check: if error description contains mostly non-printable characters
        if not is_websocket_related:
            non_printable_count = sum(1 for c in str(error) if ord(c) < 32 and c not in '\n\r\t')
            if non_printable_count > 3:  # Threshold for binary data detection
                is_websocket_related = True
        
        if is_websocket_related:
            # This is likely a WebSocket close frame or protocol error, log as debug
            logger.debug(f"WebSocket protocol frame misinterpreted as HTTP: {error}")
            return '', 400  # Return empty response for WebSocket frames
        
        # For genuine bad requests, return proper error response
        logger.warning(f"Genuine bad request: {error}")
        return jsonify({'success': False, 'error': 'Bad request.'}), 400

    @app.errorhandler(404)
    def not_found(error):
        return jsonify({'success': False, 'error': 'Endpoint not found.'}), 404

    @app.errorhandler(405)
    def method_not_allowed(error):
        return jsonify({'success': False, 'error': 'Method not allowed for this endpoint.'}), 405

    # Store socketio instance for access by other modules
    app.socketio = socketio
    
    return app, socketio


def cleanup_resources():
    """Clean up resources including the process pool and file watcher."""
    logger.info("Cleaning up resources...")
    try:
        shutdown_process_manager()
        logger.info("Process manager shut down successfully")
    except ProcessManagerError as e:
        logger.warning(f"Process manager already shut down or unavailable: {e}")
    except Exception as e:
        logger.error(f"Unexpected error shutting down process manager: {e}", exc_info=True)
    
    try:
        shutdown_websocket_watcher()
        logger.info("WebSocket file watcher shut down successfully")
    except Exception as e:
        logger.error(f"Unexpected error shutting down WebSocket file watcher: {e}", exc_info=True)


def signal_handler(signum, frame):
    """Handle shutdown signals."""
    logger.info(f"Received signal {signum}, shutting down gracefully...")
    cleanup_resources()
    sys.exit(0)


def find_available_port(host, start_port, end_port):
    """Find an available port within the specified range."""
    for port in range(start_port, end_port + 1):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind((host, port))
                return port
        except OSError:
            continue
    raise RuntimeError(f"No available port found in range {start_port}-{end_port}")


def start_development_server():
    """Start the server in development mode using Flask-SocketIO."""
    server_config = SERVER_CONFIG.get('server', {})
    dev_config = SERVER_CONFIG.get('development', {})
    socketio_config = SERVER_CONFIG.get('socketio', {})
    
    host = server_config.get('host', '127.0.0.1')
    port_config = server_config.get('port', {})
    debug = dev_config.get('debug', False)
    
    # Command line port argument takes precedence
    port_arg = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    
    if port_arg > 0:
        actual_port = port_arg
    elif port_config.get('auto_detect', True):
        port_range = port_config.get('range', {'start': 5000, 'end': 5100})
        start_port = port_config.get('default', 5000)
        actual_port = find_available_port(host, start_port, port_range['end'])
    else:
        actual_port = port_config.get('default', 5000)
    
    # Notify Electron process of the port
    print(f"FLASK_SERVER_PORT:{actual_port}", file=sys.stdout, flush=True)
    
    logger.info(f"Starting API server with Flask-SocketIO on http://{host}:{actual_port}")
    logger.info(f"Configuration: Debug={debug}, Async_mode={socketio_config.get('async_mode', 'threading')}")
    
    # Create app and socketio instances
    app, socketio = create_app()
    
    # Register cleanup functions
    atexit.register(cleanup_resources)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    # Start Flask-SocketIO server with configuration-based settings
    try:
        socketio.run(
            app, 
            host=host, 
            port=actual_port, 
            debug=debug, 
            use_reloader=False, 
            allow_unsafe_werkzeug=socketio_config.get('allow_unsafe_werkzeug', True)
        )
    except KeyboardInterrupt:
        logger.info("Server stopped by user.")
    finally:
        logger.info("Shutting down server...")
        cleanup_resources()


def create_socketio():
    """SocketIO factory for Gunicorn compatibility."""
    app, socketio = create_app()
    return socketio


# Create global app instance for import by other modules
app, socketio = create_app()


if __name__ == '__main__':
    start_development_server()