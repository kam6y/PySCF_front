# src/python/pubchem/client.py

"""PubChem API client for molecular structure retrieval."""

import requests
import logging
from typing import Optional, Dict, Any, List
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

@dataclass
class CompoundData:
    """A dataclass to hold structured compound information from PubChem."""
    cid: int
    molecular_formula: str = "Unknown"
    molecular_weight: Optional[float] = None
    iupac_name: str = "Unknown"
    synonyms: List[str] = field(default_factory=list)
    atoms: List[Dict[str, Any]] = field(default_factory=list)

class PubChemError(Exception):
    """Base exception for PubChem API errors."""
    def __init__(self, message, status_code=None):
        super().__init__(message)
        self.status_code = status_code

class PubChemNotFoundError(PubChemError):
    """Exception for 404 Not Found errors."""
    pass

class PubChemClient:
    """Client for interacting with the PubChem PUG REST API."""

    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self, timeout: int = 30):
        self.timeout = timeout
        self.session = requests.Session()

    def _make_request(self, url: str) -> requests.Response:
        """Makes a GET request and handles common errors."""
        try:
            response = self.session.get(url, timeout=self.timeout)
            
            if response.status_code == 404:
                raise PubChemNotFoundError(f"Resource not found at {url}", status_code=404)
            
            response.raise_for_status()
            return response
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Request to {url} failed: {e}")
            raise PubChemError(f"API request failed: {e}")

    def search_compound(self, query: str, search_type: str = "name") -> Optional[CompoundData]:
        try:
            cid = self._find_cid(query, search_type)
            if cid is None:
                return None

            compound_info = self._get_compound_properties(cid)
            atoms = self._get_3d_structure_from_sdf(cid)
            
            if not atoms:
                logger.warning(f"No 3D structure found for CID {cid}")
                raise PubChemNotFoundError(f"A 3D structure is not available for the compound with CID {cid}.")

            compound_info['atoms'] = atoms
            return CompoundData(**compound_info)

        except PubChemNotFoundError:
             # 見つからなかった場合はそのまま上位に伝播させる
            raise
        except (ValueError, KeyError, TypeError) as e:
            logger.error(f"Error processing data for query '{query}': {e}")
            raise PubChemError(f"Failed to process data: {e}")

    def _find_cid(self, query: str, search_type: str) -> Optional[int]:
        if search_type == "cid":
            try:
                return int(query)
            except ValueError:
                raise PubChemError(f"Invalid CID format: {query}")

        search_path = f"compound/{search_type}/{query.strip()}/cids/JSON"
        url = f"{self.BASE_URL}/{search_path}"
        
        try:
            response = self._make_request(url)
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID')
            if not cids:
                logger.info(f"No CID found for query: {query}")
                return None
            return cids[0]
        except PubChemNotFoundError:
            # 404の場合は単に見つからなかったケースとしてNoneを返す
            logger.warning(f"Could not find CID for '{query}' (type: {search_type})")
            return None

    def _get_compound_properties(self, cid: int) -> Dict[str, Any]:
        prop_url = f"{self.BASE_URL}/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight/JSON"
        prop_data = self._make_request(prop_url).json()
        
        props = prop_data.get('PropertyTable', {}).get('Properties', [{}])[0]
        
        info = {
            'cid': cid,
            'iupac_name': props.get('IUPACName'),
            'molecular_formula': props.get('MolecularFormula'),
            'molecular_weight': props.get('MolecularWeight'),
            'synonyms': []
        }

        try:
            syn_url = f"{self.BASE_URL}/compound/cid/{cid}/synonyms/JSON"
            syn_data = self._make_request(syn_url).json()
            synonyms = syn_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
            info['synonyms'] = synonyms[:5]
        except (PubChemError, KeyError, PubChemNotFoundError):
            logger.warning(f"Could not retrieve synonyms for CID {cid}. Continuing without them.")

        return info

    def _get_3d_structure_from_sdf(self, cid: int) -> Optional[List[Dict[str, Any]]]:
        url = f'{self.BASE_URL}/compound/cid/{cid}/SDF?record_type=3d'
        try:
            response = self._make_request(url)
            return self._parse_sdf_atoms(response.text)
        except PubChemNotFoundError:
            logger.info(f"No 3D SDF record found for CID {cid}.")
            return None

    def _parse_sdf_atoms(self, sdf_data: str) -> List[Dict[str, Any]]:
        atoms = []
        lines = sdf_data.strip().split('\n')
        if len(lines) < 4:
            return atoms

        try:
            counts_line = lines[3].strip()
            num_atoms = int(counts_line[:3].strip())
            
            atom_block = lines[4 : 4 + num_atoms]
            for line in atom_block:
                parts = line.split()
                if len(parts) < 4:
                    continue

                x = float(parts[0])
                y = float(parts[1])
                z = float(parts[2])
                element = parts[3]
                atoms.append({'element': element, 'x': x, 'y': y, 'z': z})
        except (ValueError, IndexError) as e:
            logger.error(f"Failed to parse SDF data: {e}")
            raise PubChemError("Could not parse SDF atom block.")
            
        return atoms