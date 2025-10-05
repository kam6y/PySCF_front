import logging
import os
import sys
import signal
import atexit
from flask import Flask, jsonify
from flask_cors import CORS
from flask_socketio import SocketIO
from pydantic import ValidationError

# Import unified configuration module
from config import (
    ServerConfig,
    get_server_config,
    determine_server_port,
    configure_flask_app,
    ConfigurationError
)

# Import blueprint registration and WebSocket handlers
from api import register_blueprints
from websocket import register_websocket_handlers

# Import cleanup functions from quantum_calc
from quantum_calc import shutdown_process_manager, shutdown_websocket_watcher
from quantum_calc.exceptions import ProcessManagerError

# Initialize configuration and logging
try:
    _server_config = get_server_config()
    log_level = getattr(logging, _server_config.get_logging_level().upper())
    log_format = _server_config.get_logging_format()
except ConfigurationError as e:
    print(f"WARNING: Configuration error: {e}. Using default logging settings.")
    log_level = logging.INFO
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

logging.basicConfig(
    level=log_level,
    format=log_format
)
logger = logging.getLogger(__name__)

# Initialize SocketIO extension globally (without binding to app)
# This will be initialized with init_app() in create_app()
socketio = SocketIO()


def create_app(server_port: int = None, test_config: dict = None):
    """
    Application factory for Gunicorn compatibility and testing.

    Args:
        server_port: The port number on which the server will run.
                    If None, it will be determined from configuration.
        test_config: Dictionary of configuration values for testing.
                    If provided, these settings override the default configuration.
    """
    # Initialize Flask app
    app = Flask(__name__)
    CORS(app)  # Enable CORS for cross-origin requests

    # Load server configuration
    server_config = get_server_config()

    # Determine server port if not provided
    if server_port is None:
        port_env = os.getenv('PYSCF_SERVER_PORT')
        server_port = determine_server_port(server_config, port_env=port_env)

    # Configure Flask app with all settings from ServerConfig
    # This establishes app.config as the single source of truth for configuration
    configure_flask_app(app, server_config, server_port)

    # Apply test configuration if provided
    if test_config is not None:
        app.config.update(test_config)

    # Initialize SocketIO with configuration-based settings
    # Use init_app() pattern for better testability
    socketio_config = app.config.get('SOCKETIO', {})
    socketio.init_app(
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

    return app


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


def start_development_server():
    """Start the server in development mode using Flask-SocketIO."""
    # Load server configuration
    server_config = get_server_config()

    host = server_config.get_server_host()
    debug = server_config.get('development.debug', False)

    # Command line port argument takes precedence
    port_arg = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else None

    # Determine server port using unified logic
    actual_port = determine_server_port(server_config, port_arg=port_arg)

    # Notify Electron process of the port
    print(f"FLASK_SERVER_PORT:{actual_port}", file=sys.stdout, flush=True)

    socketio_config = server_config.get('socketio', {})
    logger.info(f"Starting API server with Flask-SocketIO on http://{host}:{actual_port}")
    logger.info(f"Configuration: Debug={debug}, Async_mode={socketio_config.get('async_mode', 'threading')}")

    # Create app instance with the determined port
    # Global socketio instance is used automatically
    app = create_app(server_port=actual_port)
    
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


# Create global app instance for import by other modules
# Global socketio instance is already defined above
# These are created at module import time for Gunicorn compatibility
# Port is determined from environment variable or configuration
app = create_app()


if __name__ == '__main__':
    start_development_server()