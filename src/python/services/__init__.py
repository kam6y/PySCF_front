"""
Service layer module.

This module provides singleton instances of all service classes,
ensuring consistent state across API endpoints and AI agent tools.
"""

from .exceptions import (
    ServiceError,
    NotFoundError,
    ValidationError,
    ResourceUnavailableError,
    PermissionDeniedError,
    InsufficientResourcesError
)
from .quantum_service import QuantumService
from .pubchem_service import PubChemService
from .smiles_service import SMILESService
from .settings_service import SettingsService
from .system_service import SystemService

# Singleton instances
_quantum_service = None
_pubchem_service = None
_smiles_service = None
_settings_service = None
_system_service = None


def get_quantum_service() -> QuantumService:
    """
    Get the singleton QuantumService instance.
    
    Returns:
        QuantumService singleton instance
    """
    global _quantum_service
    if _quantum_service is None:
        _quantum_service = QuantumService()
    return _quantum_service


def get_pubchem_service() -> PubChemService:
    """
    Get the singleton PubChemService instance.
    
    Returns:
        PubChemService singleton instance
    """
    global _pubchem_service
    if _pubchem_service is None:
        _pubchem_service = PubChemService()
    return _pubchem_service


def get_smiles_service() -> SMILESService:
    """
    Get the singleton SMILESService instance.
    
    Returns:
        SMILESService singleton instance
    """
    global _smiles_service
    if _smiles_service is None:
        _smiles_service = SMILESService()
    return _smiles_service


def get_settings_service() -> SettingsService:
    """
    Get the singleton SettingsService instance.
    
    Returns:
        SettingsService singleton instance
    """
    global _settings_service
    if _settings_service is None:
        _settings_service = SettingsService()
    return _settings_service


def get_system_service() -> SystemService:
    """
    Get the singleton SystemService instance.
    
    Returns:
        SystemService singleton instance
    """
    global _system_service
    if _system_service is None:
        _system_service = SystemService()
    return _system_service


__all__ = [
    # Exceptions
    'ServiceError',
    'NotFoundError',
    'ValidationError',
    'ResourceUnavailableError',
    'PermissionDeniedError',
    'InsufficientResourcesError',
    # Service classes
    'QuantumService',
    'PubChemService',
    'SMILESService',
    'SettingsService',
    'SystemService',
    # Singleton getters
    'get_quantum_service',
    'get_pubchem_service',
    'get_smiles_service',
    'get_settings_service',
    'get_system_service'
]
