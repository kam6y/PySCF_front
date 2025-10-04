"""
Unified configuration management module for PySCF Native App.

This module serves as the single source of truth for all application configuration,
loading settings from config/server-config.json and providing a unified interface
for accessing configuration values throughout the application.

Configuration is loaded once at application startup and stored in Flask's app.config,
ensuring consistent configuration access across all modules.
"""

import json
import logging
import os
import socket
from pathlib import Path
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Raised when configuration loading or validation fails."""
    pass


class ServerConfig:
    """
    Centralized configuration manager for server settings.

    Loads configuration from config/server-config.json and provides
    typed access to configuration values with fallback defaults.
    """

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the configuration manager.

        Args:
            config_path: Path to server-config.json. If None, auto-detects location.
        """
        self.config_path = config_path or self._find_config_file()
        self._config: Dict[str, Any] = {}
        self._load_config()

    def _find_config_file(self) -> str:
        """
        Find the server-config.json file in the project hierarchy.

        Returns:
            str: Path to the configuration file.

        Raises:
            ConfigurationError: If config file cannot be found.
        """
        # Start from current file's directory and walk up
        current_dir = Path(__file__).parent

        # Walk up the directory tree to find config/server-config.json
        for _ in range(5):  # Limit search to 5 levels up
            config_path = current_dir / ".." / "config" / "server-config.json"
            config_path = config_path.resolve()
            if config_path.exists():
                logger.info(f"Found config file at: {config_path}")
                return str(config_path)
            current_dir = current_dir.parent

        # Fallback paths
        fallback_paths = [
            Path(__file__).parent.parent.parent / "config" / "server-config.json",
            Path.cwd() / "config" / "server-config.json",
        ]

        for path in fallback_paths:
            if path.exists():
                logger.info(f"Found config file at: {path}")
                return str(path)

        raise ConfigurationError(
            "Could not find config/server-config.json. "
            "Please ensure the configuration file exists in the project root."
        )

    def _load_config(self) -> None:
        """
        Load configuration from JSON file.

        Raises:
            ConfigurationError: If config file cannot be loaded or parsed.
        """
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                self._config = json.load(f)
            logger.info(f"Successfully loaded configuration from: {self.config_path}")
        except json.JSONDecodeError as e:
            raise ConfigurationError(f"Invalid JSON in config file {self.config_path}: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error loading config file {self.config_path}: {e}")

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a configuration value by key.

        Args:
            key: Configuration key (can use dot notation for nested keys, e.g., 'server.host')
            default: Default value if key is not found.

        Returns:
            Configuration value or default.
        """
        keys = key.split('.')
        value = self._config

        for k in keys:
            if isinstance(value, dict):
                value = value.get(k)
                if value is None:
                    return default
            else:
                return default

        return value

    def get_server_host(self) -> str:
        """Get server host address."""
        return self.get('server.host', '127.0.0.1')

    def get_default_port(self) -> int:
        """Get default server port."""
        port_config = self.get('server.port', {})
        if isinstance(port_config, dict):
            return port_config.get('default', 5000)
        return port_config if isinstance(port_config, int) else 5000

    def is_port_auto_detect_enabled(self) -> bool:
        """Check if automatic port detection is enabled."""
        port_config = self.get('server.port', {})
        if isinstance(port_config, dict):
            return port_config.get('auto_detect', True)
        return True

    def get_port_range(self) -> tuple:
        """Get port range for auto-detection."""
        port_config = self.get('server.port', {})
        if isinstance(port_config, dict):
            port_range = port_config.get('range', {'start': 5000, 'end': 5100})
            return (port_range.get('start', 5000), port_range.get('end', 5100))
        return (5000, 5100)

    def get_logging_level(self) -> str:
        """Get logging level."""
        return self.get('logging.level', 'INFO')

    def get_logging_format(self) -> str:
        """Get logging format string."""
        return self.get('logging.format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    def to_dict(self) -> Dict[str, Any]:
        """
        Get entire configuration as a dictionary.

        Returns:
            Copy of the configuration dictionary.
        """
        return self._config.copy()


def find_available_port(host: str, start_port: int, end_port: int) -> int:
    """
    Find an available port within the specified range.

    Args:
        host: Host address to bind to.
        start_port: Starting port number.
        end_port: Ending port number (inclusive).

    Returns:
        int: Available port number.

    Raises:
        RuntimeError: If no available port is found in the range.
    """
    for port in range(start_port, end_port + 1):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind((host, port))
                return port
        except OSError:
            continue
    raise RuntimeError(f"No available port found in range {start_port}-{end_port}")


def determine_server_port(config: ServerConfig, port_arg: Optional[int] = None,
                         port_env: Optional[str] = None) -> int:
    """
    Determine the server port based on configuration, arguments, and environment.

    Priority order:
    1. Command line argument (port_arg)
    2. Environment variable (port_env)
    3. Auto-detected port (if enabled)
    4. Default port from configuration

    Args:
        config: ServerConfig instance.
        port_arg: Port from command line argument.
        port_env: Port from environment variable.

    Returns:
        int: Determined port number.
    """
    # Priority 1: Command line argument
    if port_arg and port_arg > 0:
        logger.info(f"Using port from command line argument: {port_arg}")
        return port_arg

    # Priority 2: Environment variable
    if port_env:
        try:
            port = int(port_env)
            logger.info(f"Using port from environment variable: {port}")
            return port
        except (ValueError, TypeError):
            logger.warning(f"Invalid port in environment variable: {port_env}")

    # Priority 3: Auto-detection
    if config.is_port_auto_detect_enabled():
        start_port, end_port = config.get_port_range()
        default_port = config.get_default_port()
        host = config.get_server_host()

        # Start searching from the default port
        port = find_available_port(host, default_port, end_port)
        logger.info(f"Auto-detected available port: {port}")
        return port

    # Priority 4: Default port
    default_port = config.get_default_port()
    logger.info(f"Using default port: {default_port}")
    return default_port


def configure_flask_app(app, config: ServerConfig, server_port: int) -> None:
    """
    Configure Flask application with settings from ServerConfig.

    This function populates Flask's app.config with all necessary configuration
    values, establishing app.config as the single source of truth for configuration
    access throughout the application.

    Args:
        app: Flask application instance.
        config: ServerConfig instance.
        server_port: Determined server port number.
    """
    # Store entire configuration in app.config
    app.config['SERVER_CONFIG'] = config.to_dict()

    # Store commonly used values for easy access
    app.config['SERVER_HOST'] = config.get_server_host()
    app.config['SERVER_PORT'] = server_port
    app.config['DEBUG'] = config.get('server.debug', False)

    # Gunicorn settings
    app.config['GUNICORN'] = config.get('gunicorn', {})

    # SocketIO settings
    app.config['SOCKETIO'] = config.get('socketio', {})

    # Development settings
    app.config['DEVELOPMENT'] = config.get('development', {})

    # Production settings
    app.config['PRODUCTION'] = config.get('production', {})

    # Quantum calculation settings
    app.config['QUANTUM_CALCULATIONS'] = config.get('quantum_calculations', {})
    app.config['QUANTUM_CALCULATION_DEFAULTS'] = config.get('quantum_calculation_defaults', {})

    # App info
    app.config['APP_INFO'] = config.get('app_info', {})
    app.config['APP_VERSION'] = config.get('app_info.version', 'unknown')

    # External services
    app.config['EXTERNAL_SERVICES'] = config.get('external_services', {})

    # Logging settings
    app.config['LOGGING'] = config.get('logging', {})

    # AI Agent settings
    app.config['AI_AGENT'] = config.get('ai_agent', {})

    logger.info(f"Flask app configured with server port: {server_port}")


def get_fallback_config() -> Dict[str, Any]:
    """
    Get fallback configuration when config file is not available.

    Returns:
        Dictionary with default configuration values.
    """
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
            "preload_app": False,
            "access_logfile": None,
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
        "development": {
            "use_reloader": False,
            "debug": False,
            "enable_dev_tools": True
        },
        "production": {
            "use_gunicorn": True,
            "optimize_performance": True,
            "enable_logging": True
        },
        "quantum_calculations": {
            "max_concurrent_calculations": "auto",
            "process_pool_size": "auto",
            "timeout_calculation": 0,
            "memory_limit_mb": 0
        },
        "quantum_calculation_defaults": {
            "memory_settings": {
                "CASCI": 6000,
                "CASSCF": 6000,
                "CCSD": 4000,
                "CCSD_T": 4000,
                "TDDFT": 4000,
                "DFT": 2000,
                "HF": 2000,
                "MP2": 2000,
                "default": 2000
            },
            "cycle_settings": {
                "max_cycle": 150,
                "max_cycle_macro": 50,
                "max_cycle_micro": 4,
                "ah_max_cycle": 30
            },
            "spectrum_settings": {
                "ir_frequency_range": [400, 4000],
                "frequency_threshold": 80.0
            }
        },
        "app_info": {
            "version": "0.3.0"
        },
        "external_services": {
            "pubchem_timeout": 30,
            "api_timeout": 30
        },
        "logging": {
            "level": "INFO",
            "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            "enable_access_log": True,
            "enable_error_log": True
        },
        "ai_agent": {
            "model_name": "gemini-2.5-flash",
            "max_retries": 3,
            "retry_delays": [1, 2, 4],
            "retryable_status_codes": [429, 500, 502, 503, 504]
        }
    }


# Global configuration instance (lazy-loaded)
_server_config: Optional[ServerConfig] = None


def get_server_config() -> ServerConfig:
    """
    Get the global ServerConfig instance.

    Returns:
        ServerConfig: Global configuration instance.
    """
    global _server_config
    if _server_config is None:
        try:
            _server_config = ServerConfig()
        except ConfigurationError as e:
            logger.error(f"Failed to load configuration: {e}")
            logger.warning("Using fallback configuration")
            # Create a temporary config with fallback values
            _server_config = ServerConfig.__new__(ServerConfig)
            _server_config._config = get_fallback_config()
            _server_config.config_path = "fallback"
    return _server_config
