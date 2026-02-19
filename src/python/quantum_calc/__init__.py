"""Quantum chemistry calculation module for PySCF_front."""

from .base_calculator import BaseCalculator
from .dft_calculator import DFTCalculator
from .hf_calculator import HFCalculator
from .mp2_calculator import MP2Calculator
from .ccsd_calculator import CCSDCalculator
from .tddft_calculator import TDDFTCalculator
from .casci_calculator import CASCICalculator
from .casscf_calculator import CASSCFCalculator
from .orbital_generator import MolecularOrbitalGenerator
from .exceptions import CalculationError, ConvergenceError, InputError, GeometryError, FileManagerError, ProcessManagerError, WebSocketError, XYZValidationError, PauseRequestedException
from ._status_transition import CalculationStatus
from .process_manager import CalculationProcessManager, get_process_manager, initialize_process_manager_with_callback, shutdown_process_manager, update_process_manager_settings
from .file_watcher import WebSocketCalculationWatcher, get_websocket_watcher, shutdown_websocket_watcher
from .supported_parameters import get_all_supported_parameters
from .settings_manager import SettingsManager, get_settings_manager, get_current_settings, update_app_settings, mask_settings
from ._calculation_repository import CalculationRepository
from ._cube_artifact_service import CubeArtifactService
from ._calculation_directory_migration import CalculationDirectoryMigration

__all__ = [
    'BaseCalculator',
    'DFTCalculator',
    'HFCalculator',
    'MP2Calculator',
    'CCSDCalculator',
    'TDDFTCalculator',
    'CASCICalculator',
    'CASSCFCalculator',
    'MolecularOrbitalGenerator',
    'CalculationError',
    'ConvergenceError',
    'InputError',
    'GeometryError',
    'FileManagerError',
    'ProcessManagerError',
    'WebSocketError',
    'XYZValidationError',
    'PauseRequestedException',
    'CalculationStatus',
    'CalculationProcessManager',
    'get_process_manager',
    'initialize_process_manager_with_callback',
    'shutdown_process_manager',
    'update_process_manager_settings',
    'WebSocketCalculationWatcher',
    'get_websocket_watcher',
    'shutdown_websocket_watcher',
    'get_all_supported_parameters',
    'SettingsManager',
    'get_settings_manager',
    'get_current_settings',
    'update_app_settings',
    'CalculationRepository',
    'CubeArtifactService',
    'CalculationDirectoryMigration',
    'mask_settings',
]
