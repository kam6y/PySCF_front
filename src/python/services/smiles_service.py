"""
SMILES conversion service.

This service encapsulates SMILES to XYZ conversion logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
from typing import Dict, Any

from SMILES.smiles_converter import smiles_to_xyz, SMILESError
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)


class SMILESService:
    """Service for SMILES to XYZ conversion operations."""
    
    def convert_smiles(self, smiles: str, title: str = None) -> Dict[str, Any]:
        """
        Convert a SMILES string to XYZ format.
        
        Args:
            smiles: SMILES string representing the molecular structure
            title: Optional title for the XYZ file (defaults to SMILES string)
            
        Returns:
            Dict containing XYZ data
            
        Raises:
            ValidationError: If SMILES string is invalid
            ServiceError: For other errors
        """
        try:
            # Validate input
            if not smiles or not isinstance(smiles, str):
                raise ValidationError('SMILES string is required and must be a non-empty string.')
            
            if len(smiles.strip()) == 0:
                raise ValidationError('SMILES string cannot be empty or contain only whitespace.')
            
            # SMILES strings typically don't exceed 500 characters for reasonable molecules
            if len(smiles) > 500:
                raise ValidationError('SMILES string is too long. Maximum length is 500 characters.')
            
            logger.info(f"Converting SMILES: {smiles}")
            
            # Use provided title or default
            if title is None:
                title = f"Molecule from SMILES: {smiles}"
            
            xyz_string = smiles_to_xyz(smiles, title=title)
            
            logger.info(f"Successfully converted SMILES to XYZ structure")
            
            return {'xyz': xyz_string}
            
        except SMILESError as e:
            logger.error(f"SMILES conversion failed: {e}")
            raise ValidationError(str(e))
        except Exception as e:
            logger.error(f"An unexpected error occurred during SMILES conversion: {e}", exc_info=True)
            raise ServiceError('An internal server error occurred.')
