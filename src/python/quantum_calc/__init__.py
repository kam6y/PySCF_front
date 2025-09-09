"""Quantum chemistry calculation module for PySCF Native App."""

from .base_calculator import BaseCalculator
from .dft_calculator import DFTCalculator
from .hf_calculator import HFCalculator
from .mp2_calculator import MP2Calculator
from .ccsd_calculator import CCSDCalculator
from .tddft_calculator import TDDFTCalculator
from .orbital_generator import MolecularOrbitalGenerator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError, FileManagerError, ProcessManagerError, WebSocketError, XYZValidationError
from .process_manager import CalculationProcessManager, get_process_manager, shutdown_process_manager, update_process_manager_settings
from .file_watcher import WebSocketCalculationWatcher, get_websocket_watcher, shutdown_websocket_watcher
from .supported_parameters import get_all_supported_parameters
from .settings_manager import SettingsManager, get_settings_manager, get_current_settings, update_app_settings

__all__ = [
    'BaseCalculator',
    'DFTCalculator',
    'HFCalculator',
    'MP2Calculator',
    'CCSDCalculator',
    'TDDFTCalculator',
    'MolecularOrbitalGenerator',
    'CalculationError',
    'ConvergenceError',
    'InputError',
    'GeometryError',
    'FileManagerError',
    'ProcessManagerError',
    'WebSocketError',
    'XYZValidationError',
    'CalculationProcessManager',
    'get_process_manager',
    'shutdown_process_manager',
    'update_process_manager_settings',
    'WebSocketCalculationWatcher',
    'get_websocket_watcher',
    'shutdown_websocket_watcher',
    'get_all_supported_parameters',
    'SettingsManager',
    'get_settings_manager',
    'get_current_settings',
    'update_app_settings'
]