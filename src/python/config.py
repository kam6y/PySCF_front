"""
Unified configuration management module for PySCF_front.

This module serves as the single source of truth for all application configuration,
loading settings from config/server-config.json and providing a unified interface
for accessing configuration values throughout the application.

Configuration is loaded once at application startup and stored on FastAPI app.state,
ensuring consistent configuration access across all modules.
"""

import json
import logging
import os
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
        logger.debug("Starting config file search...")
        search_paths = []

        # Priority 1: Environment variable (for packaged apps)
        # main.ts sets PYSCF_RESOURCES_PATH to process.resourcesPath in packaged mode
        resources_path = os.getenv('PYSCF_RESOURCES_PATH')
        if resources_path:
            env_config_path = Path(resources_path) / "config" / "server-config.json"
            search_paths.append(("Environment variable (PYSCF_RESOURCES_PATH)", env_config_path))
            if env_config_path.exists():
                logger.info(f"Found config file at: {env_config_path} (from PYSCF_RESOURCES_PATH)")
                return str(env_config_path)

        # Priority 2: Walk up the directory tree (for development)
        current_dir = Path(__file__).parent
        for i in range(5):  # Limit search to 5 levels up
            config_path = current_dir / ".." / "config" / "server-config.json"
            config_path = config_path.resolve()
            search_paths.append((f"Parent directory level {i+1}", config_path))
            if config_path.exists():
                logger.info(f"Found config file at: {config_path} (development mode)")
                return str(config_path)
            current_dir = current_dir.parent

        # Priority 3: Additional fallback paths
        fallback_paths = [
            ("Relative from __file__", Path(__file__).parent.parent.parent / "config" / "server-config.json"),
            ("Current working directory", Path.cwd() / "config" / "server-config.json"),
        ]
        search_paths.extend(fallback_paths)

        for description, path in fallback_paths:
            if path.exists():
                logger.info(f"Found config file at: {path} ({description})")
                return str(path)

        # Config file not found - log all attempted paths for debugging
        logger.error("Config file not found. Searched the following locations:")
        for description, path in search_paths:
            logger.error(f"  - {description}: {path} (exists: {path.exists()})")

        raise ConfigurationError(
            "Could not find config/server-config.json. "
            "Please ensure the configuration file exists in the project root. "
            f"Searched {len(search_paths)} locations."
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


def determine_server_port(config: ServerConfig, port_arg: Optional[int] = None,
                         port_env: Optional[str] = None) -> int:
    """
    Determine the server port from environment or arguments.

    Port detection is handled by Electron in production mode.
    In standalone mode (python app.py), uses command line or default.

    Priority order:
    1. Command line argument (port_arg) - for `python app.py 5001`
    2. Environment variable (port_env) - set by Electron (PYSCF_SERVER_PORT)
    3. Default port from configuration - standalone fallback

    Args:
        config: ServerConfig instance.
        port_arg: Port from command line argument.
        port_env: Port from environment variable (PYSCF_SERVER_PORT).

    Returns:
        int: Determined port number.
    """
    # Priority 1: Command line argument
    if port_arg and port_arg > 0:
        logger.info(f"Using port from CLI: {port_arg}")
        return port_arg

    # Priority 2: Environment variable (set by Electron)
    if port_env:
        try:
            port = int(port_env)
            logger.info(f"Using port from Electron: {port}")
            return port
        except (ValueError, TypeError):
            logger.warning(f"Invalid port in env: {port_env}")

    # Priority 3: Default (standalone fallback)
    default_port = config.get_default_port()
    logger.info(f"Using default port: {default_port}")
    return default_port


def configure_fastapi_app(app, config: ServerConfig, server_port: int) -> None:
    """
    Configure FastAPI application state from ServerConfig.

    FastAPI has no app.config dict, so application settings are stored on
    app.state with the same logical keys the Flask runtime used.
    """
    app.state.SERVER_CONFIG = config.to_dict()
    app.state.SERVER_CONFIG_OBJECT = config
    app.state.SERVER_HOST = config.get_server_host()
    app.state.SERVER_PORT = server_port
    app.state.DEBUG = config.get('server.debug', False)
    app.state.TESTING = False

    app.state.GUNICORN = config.get('gunicorn', {})
    app.state.SOCKETIO = config.get('socketio', {})
    app.state.DEVELOPMENT = config.get('development', {})
    app.state.PRODUCTION = config.get('production', {})
    app.state.QUANTUM_CALCULATIONS = config.get('quantum_calculations', {})
    app.state.QUANTUM_CALCULATION_DEFAULTS = config.get('quantum_calculation_defaults', {})
    app.state.APP_INFO = config.get('app_info', {})
    app.state.APP_VERSION = config.get('app_info.version', 'unknown')
    app.state.EXTERNAL_SERVICES = config.get('external_services', {})
    app.state.LOGGING = config.get('logging', {})
    app.state.AI_AGENT = config.get('ai_agent', {})
    app.state.VERSION = config.get('app_info.version', 'unknown')

    logger.info(f"FastAPI app configured with server port: {server_port}")


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
            logger.info("Successfully initialized ServerConfig from file")
        except ConfigurationError as e:
            logger.critical(f"Failed to load configuration: {e}")
            # Re-raise the exception to prevent the application from starting with invalid state
            raise
    return _server_config
