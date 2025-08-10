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