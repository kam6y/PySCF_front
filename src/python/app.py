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
    
    # Initialize process manager with WebSocket notification callback
    from quantum_calc import initialize_process_manager_with_callback
    try:
        initialize_process_manager_with_callback(
            notification_callback=websocket_funcs['send_immediate_websocket_notification']
        )
        logger.info("Process manager initialized with WebSocket notification callback")
    except Exception as e:
        logger.error(f"Failed to initialize process manager with callback: {e}")
        # Continue anyway - process manager will work without notifications

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
        """
        Handle bad requests, filtering out WebSocket close frame errors.
        
        KNOWN ISSUE & WORKAROUND:
        When WebSocket connections close, the close frame (binary data) is sometimes
        misinterpreted as an HTTP request by Werkzeug/Flask-SocketIO, resulting in
        400 Bad Request errors with messages like "Invalid HTTP method" or garbled
        binary data in the error description.
        
        This is a known issue in the Flask-SocketIO/Werkzeug stack when using
        threading async_mode. The issue has been observed across multiple versions
        and environments (see Flask-SocketIO issues #287, #466, #1417, #1811).
        
        ROOT CAUSE:
        - WebSocket close frames contain binary protocol data
        - When connection teardown occurs, these frames may be read by HTTP handlers
        - The binary data fails to parse as valid HTTP, triggering 400 errors
        
        CURRENT WORKAROUND:
        This handler detects WebSocket-related errors by examining the error message
        for known patterns (protocol keywords, binary data indicators) and silently
        logs them as debug messages rather than warnings to avoid log pollution.
        
        FUTURE MONITORING:
        - Monitor Flask-SocketIO and Werkzeug changelogs for protocol handling fixes
        - Consider upgrading to newer async modes (eventlet/gevent) if issues persist
        - Review this workaround when upgrading major versions of dependencies
        
        VALIDATION:
        This approach has been validated as the recommended workaround by the
        Flask-SocketIO community and is safe as it only affects cosmetic logging,
        not functionality.
        """
        error_description = str(error).lower()
        
        # Comprehensive list of WebSocket protocol error indicators
        # Based on observed patterns from Flask-SocketIO issues and WebSocket RFC 6455
        websocket_indicators = [
            # HTTP method errors (most common)
            'invalid http method', 'expected get method', 'invalid method',
            # WebSocket protocol keywords
            'websocket', 'connection upgrade', 'upgrade required',
            # Request parsing errors
            'bad request line', 'malformed request', 'protocol error',
            # Connection state errors
            'connection reset', 'connection closed', 'connection aborted',
            # Encoding errors (binary data misinterpreted as text)
            'invalid utf-8', 'decode error', 'unicode error'
        ]
        
        # Primary detection: Check for known error message patterns
        is_websocket_related = any(indicator in error_description for indicator in websocket_indicators)
        
        # Secondary detection: Binary data heuristic
        # WebSocket close frames contain non-printable bytes that appear in error messages
        if not is_websocket_related:
            error_str = str(error)
            non_printable_count = sum(1 for c in error_str if ord(c) < 32 and c not in '\n\r\t')
            # If error message contains significant binary data, likely a WebSocket frame
            if non_printable_count > 3:  # Empirically determined threshold
                is_websocket_related = True
        
        if is_websocket_related:
            # WebSocket protocol frame misinterpreted as HTTP - expected behavior
            # Log at debug level to avoid polluting logs with normal connection teardown
            logger.debug(f"WebSocket close frame detected (expected): {error}")
            return '', 400  # Return minimal response
        
        # Genuine HTTP 400 error - log and return proper error response
        logger.warning(f"Genuine bad HTTP request: {error}")
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


# Create global app and socketio instances for import by other modules
# These are created at module import time for Gunicorn compatibility
app, socketio = create_app()


if __name__ == '__main__':
    start_development_server()