"""Custom exceptions for quantum chemistry calculations."""


class CalculationError(Exception):
    """Base exception for quantum chemistry calculation errors."""
    pass


class ConvergenceError(CalculationError):
    """Exception raised when SCF or optimization fails to converge."""
    pass


class InputError(CalculationError):
    """Exception raised for invalid input parameters."""
    pass


class GeometryError(CalculationError):
    """Exception raised for invalid molecular geometry."""
    pass


class FileManagerError(Exception):
    """Exception raised for file management operations."""
    pass


class ProcessManagerError(Exception):
    """Exception raised for process pool management operations."""
    pass


class WebSocketError(Exception):
    """Exception raised for WebSocket communication errors."""
    pass


class XYZValidationError(Exception):
    """Exception raised for XYZ format validation errors."""
    pass