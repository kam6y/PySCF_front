"""PubChem API client for molecular structure retrieval."""

import pubchempy as pcp
from typing import Optional, Dict, Any, List, Tuple
import requests
import logging
import urllib3
import ssl

# Suppress SSL warnings for development
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Configure SSL for development - disable certificate verification
import certifi
import os
os.environ['REQUESTS_CA_BUNDLE'] = ''
os.environ['SSL_VERIFY'] = '0'

# Monkey patch SSL verification for pubchempy
try:
    # Patch urllib for pubchempy
    import urllib.request
    import urllib.error
    
    # Create unverified SSL context for urllib
    original_urlopen = urllib.request.urlopen
    
    def unverified_urlopen(url, data=None, timeout=None, *args, **kwargs):
        if isinstance(url, str) and url.startswith('https://'):
            # Create unverified context
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            kwargs['context'] = ctx
        elif hasattr(url, 'full_url') and url.full_url.startswith('https://'):
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            kwargs['context'] = ctx
        return original_urlopen(url, data, timeout, *args, **kwargs)
    
    urllib.request.urlopen = unverified_urlopen
    
    # Also patch requests session for good measure
    import requests.adapters
    from requests.packages.urllib3.util.ssl_ import create_urllib3_context
    
    class NoSSLVerifyHTTPAdapter(requests.adapters.HTTPAdapter):
        def init_poolmanager(self, *args, **pool_kwargs):
            pool_kwargs['ssl_context'] = create_urllib3_context()
            pool_kwargs['ssl_context'].check_hostname = False
            pool_kwargs['ssl_context'].verify_mode = ssl.CERT_NONE
            return super().init_poolmanager(*args, **pool_kwargs)
    
    # Create a session with SSL verification disabled
    pcp.session = requests.Session()
    pcp.session.mount('https://', NoSSLVerifyHTTPAdapter())
    pcp.session.verify = False
    
except Exception as e:
    # If patching fails, log but continue
    print(f"SSL patching failed: {e}")
    pass

logger = logging.getLogger(__name__)


class PubChemError(Exception):
    """Custom exception for PubChem API errors."""
    pass


class PubChemClient:
    """Client for interacting with PubChem API to retrieve molecular structures."""
    
    def __init__(self, timeout: int = 30):
        """Initialize the PubChem client.
        
        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
    
    def search_compound(self, query: str, search_type: str = "name") -> Optional[pcp.Compound]:
        """Search for a compound by name or CID.
        
        Args:
            query: Compound name or CID
            search_type: Search type ("name", "cid", or "formula")
            
        Returns:
            PubChem Compound object or None if not found
            
        Raises:
            PubChemError: If the search fails
        """
        try:
            # For SSL issues, try using PubChem REST API directly
            return self._direct_pubchem_search(query, search_type)
            
        except Exception as e:
            logger.error(f"Error searching for compound {query}: {e}")
            raise PubChemError(f"Failed to search compound: {e}")
    
    def _direct_pubchem_search(self, query: str, search_type: str) -> Optional[pcp.Compound]:
        """Direct PubChem REST API search to bypass SSL issues."""
        try:
            # First, get the CID
            if search_type == "cid":
                try:
                    cid = int(query)
                except ValueError:
                    raise PubChemError(f"Invalid CID format: {query}")
            else:
                # Search for CID by name using REST API
                search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/cids/JSON"
                response = requests.get(search_url, timeout=self.timeout, verify=False)
                
                if response.status_code != 200:
                    return None
                
                data = response.json()
                if 'IdentifierList' not in data or 'CID' not in data['IdentifierList']:
                    return None
                
                cid = data['IdentifierList']['CID'][0]
            
            # Get compound information including name
            compound_info = self._get_compound_info_by_cid(cid)
            
            # Create a mock compound object with the CID and name
            class MockCompound:
                def __init__(self, cid, info=None):
                    self.cid = cid
                    self.atoms = []
                    if info:
                        self.iupac_name = info.get('iupac_name', 'Unknown')
                        self.molecular_formula = info.get('molecular_formula', 'Unknown')
                        self.molecular_weight = info.get('molecular_weight', None)
                        self.synonyms = info.get('synonyms', [])
                    else:
                        self.iupac_name = 'Unknown'
                        self.molecular_formula = 'Unknown'
                        self.molecular_weight = None
                        self.synonyms = []
                    
            return MockCompound(cid, compound_info)
            
        except Exception as e:
            logger.error(f"Direct PubChem search failed: {e}")
            return None
    
    def get_3d_structure(self, compound: pcp.Compound) -> Optional[List[Dict[str, Any]]]:
        """Get 3D atomic coordinates from a compound.
        
        Args:
            compound: PubChem compound object
            
        Returns:
            List of atom dictionaries with coordinates, or None if no 3D data
        """
        try:
            if not hasattr(compound, 'atoms') or not compound.atoms:
                # Try to get 3D data via SDF
                return self._get_3d_from_sdf(compound.cid)
            
            atoms = []
            for atom in compound.atoms:
                atom_data = {
                    'element': atom.element,
                    'x': getattr(atom, 'x', None),
                    'y': getattr(atom, 'y', None),
                    'z': getattr(atom, 'z', None)
                }
                # Only include atoms with 3D coordinates
                if all(coord is not None for coord in [atom_data['x'], atom_data['y'], atom_data['z']]):
                    atoms.append(atom_data)
            
            return atoms if atoms else None
            
        except Exception as e:
            logger.error(f"Error getting 3D structure: {e}")
            return None
    
    def _get_3d_from_sdf(self, cid: int) -> Optional[List[Dict[str, Any]]]:
        """Get 3D structure from SDF format.
        
        Args:
            cid: Compound ID
            
        Returns:
            List of atom dictionaries with coordinates
        """
        try:
            url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d'
            # Disable SSL verification for development to avoid certificate issues
            response = requests.get(url, timeout=self.timeout, verify=False)
            response.raise_for_status()
            
            sdf_data = response.text
            return self._parse_sdf_atoms(sdf_data)
            
        except Exception as e:
            logger.error(f"Error getting 3D SDF for CID {cid}: {e}")
            return None
    
    def _parse_sdf_atoms(self, sdf_data: str) -> List[Dict[str, Any]]:
        """Parse atom coordinates from SDF data.
        
        Args:
            sdf_data: SDF file content
            
        Returns:
            List of atom dictionaries
        """
        atoms = []
        lines = sdf_data.split('\n')
        
        # Find the counts line (4th line in SDF, skip empty lines)
        counts_line_idx = None
        for i, line in enumerate(lines):
            if line.strip() and i >= 3:  # Find first non-empty line after header
                counts_line_idx = i
                break
        
        if counts_line_idx is None:
            return atoms
            
        try:
            counts_line = lines[counts_line_idx]
            if len(counts_line.strip()) < 3:
                return atoms
            num_atoms = int(counts_line[:3].strip())
            
            # Parse atom block (starts after counts line)
            for i in range(counts_line_idx + 1, counts_line_idx + 1 + num_atoms):
                if i < len(lines):
                    line = lines[i]
                    if len(line) >= 31:  # Minimum length for coordinate data
                        x = float(line[0:10].strip())
                        y = float(line[10:20].strip())
                        z = float(line[20:30].strip())
                        element = line[31:34].strip()
                        
                        atoms.append({
                            'element': element,
                            'x': x,
                            'y': y,
                            'z': z
                        })
        except (ValueError, IndexError) as e:
            logger.error(f"Error parsing SDF data: {e}")
        
        return atoms
    
    def _get_compound_info_by_cid(self, cid: int) -> Dict[str, Any]:
        """Get compound information directly by CID using REST API."""
        try:
            # Get compound properties
            prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,MolecularFormula,MolecularWeight/JSON"
            response = requests.get(prop_url, timeout=self.timeout, verify=False)
            
            info = {
                'cid': cid,
                'iupac_name': 'Unknown',
                'molecular_formula': 'Unknown',
                'molecular_weight': None,
                'synonyms': []
            }
            
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    props = data['PropertyTable']['Properties'][0]
                    info['iupac_name'] = props.get('IUPACName', 'Unknown')
                    info['molecular_formula'] = props.get('MolecularFormula', 'Unknown')
                    info['molecular_weight'] = props.get('MolecularWeight', None)
            
            # Try to get synonyms (common names)
            try:
                syn_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
                syn_response = requests.get(syn_url, timeout=self.timeout, verify=False)
                if syn_response.status_code == 200:
                    syn_data = syn_response.json()
                    if 'InformationList' in syn_data and 'Information' in syn_data['InformationList']:
                        synonyms = syn_data['InformationList']['Information'][0].get('Synonym', [])
                        info['synonyms'] = synonyms[:3]  # First 3 synonyms
            except:
                pass  # Synonyms are optional
            
            return info
            
        except Exception as e:
            logger.error(f"Error getting compound info for CID {cid}: {e}")
            return {
                'cid': cid,
                'iupac_name': 'Unknown',
                'molecular_formula': 'Unknown',
                'molecular_weight': None,
                'synonyms': []
            }

    def get_compound_info(self, compound: pcp.Compound) -> Dict[str, Any]:
        """Get basic information about a compound.
        
        Args:
            compound: PubChem compound object
            
        Returns:
            Dictionary with compound information
        """
        # Handle synonyms safely (may be Mock object in tests)
        synonyms = getattr(compound, 'synonyms', [])
        if hasattr(synonyms, '__getitem__') and hasattr(synonyms, '__len__'):
            synonyms = synonyms[:3]
        else:
            synonyms = []
            
        info = {
            'cid': compound.cid,
            'molecular_formula': getattr(compound, 'molecular_formula', 'Unknown'),
            'molecular_weight': getattr(compound, 'molecular_weight', None),
            'iupac_name': getattr(compound, 'iupac_name', 'Unknown'),
            'synonyms': synonyms
        }
        return info