"""
Service layer exceptions.

This module defines common exceptions used across all service classes
to provide consistent error handling between the service layer and API endpoints.
"""


class ServiceError(Exception):
    """Base exception for all service layer errors."""
    
    def __init__(self, message: str, status_code: int = 500):
        """
        Initialize ServiceError.
        
        Args:
            message: Human-readable error message
            status_code: HTTP status code that should be returned (default: 500)
        """
        super().__init__(message)
        self.message = message
        self.status_code = status_code


class NotFoundError(ServiceError):
    """Exception raised when a requested resource is not found."""
    
    def __init__(self, message: str):
        """
        Initialize NotFoundError.
        
        Args:
            message: Human-readable error message
        """
        super().__init__(message, status_code=404)


class ValidationError(ServiceError):
    """Exception raised when input validation fails."""
    
    def __init__(self, message: str):
        """
        Initialize ValidationError.
        
        Args:
            message: Human-readable error message
        """
        super().__init__(message, status_code=400)


class ResourceUnavailableError(ServiceError):
    """Exception raised when a required resource is temporarily unavailable."""
    
    def __init__(self, message: str):
        """
        Initialize ResourceUnavailableError.
        
        Args:
            message: Human-readable error message
        """
        super().__init__(message, status_code=503)


class PermissionDeniedError(ServiceError):
    """Exception raised when access to a resource is denied."""
    
    def __init__(self, message: str):
        """
        Initialize PermissionDeniedError.
        
        Args:
            message: Human-readable error message
        """
        super().__init__(message, status_code=403)


class InsufficientResourcesError(ServiceError):
    """Exception raised when insufficient system resources are available."""
    
    def __init__(self, message: str):
        """
        Initialize InsufficientResourcesError.
        
        Args:
            message: Human-readable error message
        """
        super().__init__(message, status_code=507)
