"""
PubChem search service.

This service encapsulates PubChem molecular data retrieval and XYZ validation logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
from typing import Dict, Any

from pubchem.client import PubChemClient, PubChemError, PubChemNotFoundError
from pubchem import parser as xyz_parser
from .exceptions import ServiceError, NotFoundError, ValidationError

logger = logging.getLogger(__name__)


class PubChemService:
    """Service for PubChem molecular data operations."""
    
    def __init__(self, timeout: int = 30):
        """
        Initialize PubChemService.
        
        Args:
            timeout: Request timeout in seconds
        """
        self.client = PubChemClient(timeout=timeout)
    
    def search_compound(self, query: str, search_type: str = 'name') -> Dict[str, Any]:
        """
        Search PubChem for a compound and return its 3D structure in XYZ format.
        
        Args:
            query: Search query (compound name, CID, or formula)
            search_type: Type of search ('name', 'cid', or 'formula')
            
        Returns:
            Dict containing XYZ data and compound information
            
        Raises:
            NotFoundError: If compound not found
            ValidationError: If search type is invalid
            ServiceError: For other errors
        """
        try:
            # Validate search type
            valid_types = ['name', 'cid', 'formula']
            if search_type not in valid_types:
                raise ValidationError(f"Invalid search type. Must be one of: {', '.join(valid_types)}")
            
            logger.info(f"Searching PubChem for '{query}' (type: {search_type})")
            
            compound_data = self.client.search_compound(query, search_type)
            if not compound_data or not compound_data.atoms:
                raise NotFoundError(f'No compound with a 3D structure found for query: {query}')
            
            logger.info(f"Found CID {compound_data.cid} with {len(compound_data.atoms)} atoms.")
            
            title = xyz_parser.format_compound_title(compound_data, query)
            xyz_string = xyz_parser.atoms_to_xyz(compound_data.atoms, title)
            
            logger.info(f"Successfully generated XYZ for CID {compound_data.cid}")
            
            return {
                'xyz': xyz_string,
                'compound_info': {
                    'cid': compound_data.cid,
                    'iupac_name': compound_data.iupac_name,
                    'molecular_formula': compound_data.molecular_formula,
                    'molecular_weight': compound_data.molecular_weight,
                    'synonyms': compound_data.synonyms,
                },
                'atom_count': len(compound_data.atoms)
            }
            
        except ValidationError:
            raise  # Re-raise ValidationError to propagate to caller
        except NotFoundError:
            raise  # Re-raise NotFoundError to propagate to caller
        except PubChemNotFoundError as e:
            logger.warning(f"PubChem search failed (Not Found): {e}")
            raise NotFoundError(str(e))
        except PubChemError as e:
            logger.error(f"A PubChem API error occurred: {e}", exc_info=True)
            # Map PubChem error status codes to service exceptions
            if hasattr(e, 'status_code'):
                if e.status_code == 404:
                    raise NotFoundError(str(e))
                elif e.status_code == 400:
                    raise ValidationError(str(e))
            raise ServiceError(str(e))
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}", exc_info=True)
            raise ServiceError('An internal server error occurred.')
    
    def validate_xyz(self, xyz_string: str) -> Dict[str, Any]:
        """
        Validate an XYZ format string.

        Args:
            xyz_string: XYZ format string to validate

        Returns:
            Dict containing validation results

        Raises:
            ValidationError: If xyz_string is empty or invalid format
            ServiceError: For other errors
        """
        try:
            # Check for empty or None input
            if not xyz_string or not xyz_string.strip():
                raise ValidationError('XYZ string cannot be empty.')

            validation_result = xyz_parser.validate_xyz(xyz_string)

            return validation_result

        except ValidationError:
            raise  # Re-raise ValidationError to propagate to caller
        except (TypeError, AttributeError) as e:
            logger.warning(f"Invalid input data for XYZ validation: {e}")
            raise ValidationError('Invalid input format for XYZ validation.')
        except MemoryError as e:
            logger.error(f"Memory error during XYZ validation: {e}")
            raise ServiceError('XYZ string too large to process.', status_code=413)
        except Exception as e:
            logger.error(f"Unexpected error during XYZ validation: {e}", exc_info=True)
            raise ServiceError('An internal server error occurred.')

    def get_compound_smiles(self, query: str, search_type: str = 'name') -> Dict[str, Any]:
        """
        Get the canonical SMILES and properties for a compound from PubChem.

        Args:
            query: Search query (compound name, CID, or formula)
            search_type: Type of search ('name', 'cid', or 'formula')

        Returns:
            Dict containing SMILES and compound information:
            {
                'smiles': str,
                'compound_info': {
                    'cid': int,
                    'iupac_name': str,
                    'molecular_formula': str,
                    'molecular_weight': float
                }
            }

        Raises:
            NotFoundError: If compound not found
            ValidationError: If search type is invalid
            ServiceError: For other errors
        """
        try:
            # Validate search type
            valid_types = ['name', 'cid', 'formula']
            if search_type not in valid_types:
                raise ValidationError(f"Invalid search type. Must be one of: {', '.join(valid_types)}")

            logger.info(f"Retrieving SMILES from PubChem for '{query}' (type: {search_type})")

            result = self.client.get_compound_smiles(query, search_type)

            if not result or not result.get('smiles'):
                raise NotFoundError(f'No SMILES found for query: {query}')

            logger.info(f"Successfully retrieved SMILES for '{query}': {result['smiles']}")

            return {
                'smiles': result['smiles'],
                'compound_info': {
                    'cid': result['cid'],
                    'iupac_name': result['iupac_name'],
                    'molecular_formula': result['molecular_formula'],
                    'molecular_weight': result['molecular_weight']
                }
            }

        except ValidationError:
            raise  # Re-raise ValidationError to propagate to caller
        except NotFoundError:
            raise  # Re-raise NotFoundError to propagate to caller
        except PubChemNotFoundError as e:
            logger.warning(f"PubChem SMILES search failed (Not Found): {e}")
            raise NotFoundError(str(e))
        except PubChemError as e:
            logger.error(f"A PubChem API error occurred: {e}", exc_info=True)
            # Map PubChem error status codes to service exceptions
            if hasattr(e, 'status_code'):
                if e.status_code == 404:
                    raise NotFoundError(str(e))
                elif e.status_code == 400:
                    raise ValidationError(str(e))
            raise ServiceError(str(e))
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}", exc_info=True)
            raise ServiceError('An internal server error occurred.')
