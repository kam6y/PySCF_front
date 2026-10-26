"""Client module for retrieving molecular structures from the PubChem API."""

# Re-export the main client and data class from the package root.
from .client import PubChemClient, PubChemError, CompoundData

# Make the parser module available as pubchem.parser.
from . import parser

__version__ = "0.2.0"

__all__ = [
    "PubChemClient",
    "PubChemError",
    "CompoundData",
    "parser"
]
