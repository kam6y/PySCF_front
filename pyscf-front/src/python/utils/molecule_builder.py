"""
Molecule Builder Utility
Functions for building PySCF molecule objects from various input formats.
"""

import logging
from typing import Dict, Any, List, Tuple, Optional
import numpy as np
from pyscf import gto

logger = logging.getLogger(__name__)


class MoleculeBuilder:
    """Utility class for building PySCF molecule objects from various input formats."""
    
    def __init__(self):
        pass
    
    def build_from_data(self, data: Dict[str, Any]) -> gto.Mole:
        """
        Build a PySCF molecule from input data.
        
        Args:
            data: Dictionary containing molecule information
                  Supported formats:
                  - 'coordinates': List of [element, x, y, z] coordinates
                  - 'smiles': SMILES string (requires RDKit)
                  - 'xyz': XYZ format string
                  - 'pdb': PDB format string
        
        Returns:
            PySCF Mole object
        """
        input_type = data.get('type', 'coordinates')
        
        if input_type == 'coordinates':
            return self._build_from_coordinates(data)
        elif input_type == 'smiles':
            return self._build_from_smiles(data)
        elif input_type == 'xyz':
            return self._build_from_xyz(data)
        elif input_type == 'pdb':
            return self._build_from_pdb(data)
        else:
            raise ValueError(f"Unsupported input type: {input_type}")
    
    def _build_from_coordinates(self, data: Dict[str, Any]) -> gto.Mole:
        """Build molecule from coordinate list."""
        coordinates = data.get('coordinates', [])
        charge = data.get('charge', 0)
        spin = data.get('spin', 0)
        unit = data.get('unit', 'Angstrom')  # Default to Angstrom, allow override
        
        if not coordinates:
            raise ValueError("No coordinates provided")
        
        # Convert coordinates to PySCF format
        atom_string = ""
        for coord in coordinates:
            if len(coord) != 4:
                raise ValueError("Each coordinate must be [element, x, y, z]")
            element, x, y, z = coord
            atom_string += f"{element} {x:.6f} {y:.6f} {z:.6f}; "
        
        # Remove trailing semicolon and space
        atom_string = atom_string.rstrip('; ')
        
        mol = gto.Mole()
        mol.atom = atom_string
        mol.charge = charge
        mol.spin = spin
        mol.unit = unit  # Use specified unit
        mol.build()
        
        logger.info(f"Built molecule from coordinates: {len(coordinates)} atoms")
        logger.info(f"PySCF molecule unit: {mol.unit}")
        return mol
    
    def get_coordinates_in_angstrom(self, mol: gto.Mole) -> List[List[float]]:
        """
        Get molecule coordinates in Angstrom units.
        
        PySCF's atom_coords() always returns coordinates in Bohr units,
        regardless of the input unit setting. This method converts them back to Angstrom.
        """
        # PySCF conversion factor from Bohr to Angstrom
        BOHR_TO_ANGSTROM = 0.5291772083
        
        coords_bohr = mol.atom_coords()
        coords_angstrom = coords_bohr * BOHR_TO_ANGSTROM
        
        # Return as list of [element, x, y, z] format
        result = []
        for i in range(mol.natm):
            element = mol.atom_symbol(i)
            x, y, z = coords_angstrom[i]
            result.append([element, float(x), float(y), float(z)])
        
        return result
    
    def _build_from_smiles(self, data: Dict[str, Any]) -> gto.Mole:
        """Build molecule from SMILES string."""
        smiles = data.get('smiles', '')
        charge = data.get('charge', 0)
        spin = data.get('spin', 0)
        
        if not smiles:
            raise ValueError("No SMILES string provided")
        
        try:
            # Try to use RDKit if available
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Parse SMILES
            mol_rdkit = Chem.MolFromSmiles(smiles)
            if mol_rdkit is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
            
            # Add hydrogens
            mol_rdkit = Chem.AddHs(mol_rdkit)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol_rdkit)
            AllChem.UFFOptimizeMolecule(mol_rdkit)
            
            # Extract coordinates
            conf = mol_rdkit.GetConformer()
            coordinates = []
            
            for i, atom in enumerate(mol_rdkit.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                coordinates.append([
                    atom.GetSymbol(),
                    pos.x,
                    pos.y,
                    pos.z
                ])
            
            # Build PySCF molecule
            return self._build_from_coordinates({
                'coordinates': coordinates,
                'charge': charge,
                'spin': spin
            })
            
        except ImportError:
            raise RuntimeError("RDKit not available. Cannot process SMILES strings.")
        except Exception as e:
            logger.error(f"Error processing SMILES: {e}")
            raise
    
    def _build_from_xyz(self, data: Dict[str, Any]) -> gto.Mole:
        """Build molecule from XYZ format string."""
        xyz_string = data.get('xyz', '')
        charge = data.get('charge', 0)
        spin = data.get('spin', 0)
        
        if not xyz_string:
            raise ValueError("No XYZ string provided")
        
        lines = xyz_string.strip().split('\n')
        
        if len(lines) < 2:
            raise ValueError("Invalid XYZ format")
        
        try:
            num_atoms = int(lines[0])
            # Skip comment line (line 1)
            
            coordinates = []
            for i in range(2, 2 + num_atoms):
                if i >= len(lines):
                    raise ValueError("Insufficient coordinate lines in XYZ")
                
                parts = lines[i].split()
                if len(parts) < 4:
                    raise ValueError(f"Invalid coordinate line: {lines[i]}")
                
                element = parts[0]
                x, y, z = map(float, parts[1:4])
                coordinates.append([element, x, y, z])
            
            # Convert to internal coordinate format and set unit explicitly  
            coord_data = {
                'coordinates': coordinates,
                'charge': charge,
                'spin': spin,
                'unit': 'Angstrom'  # XYZ files typically use Angstrom units
            }
            return self._build_from_coordinates(coord_data)
            
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error parsing XYZ format: {e}")
    
    def _build_from_pdb(self, data: Dict[str, Any]) -> gto.Mole:
        """Build molecule from PDB format string."""
        pdb_string = data.get('pdb', '')
        charge = data.get('charge', 0)
        spin = data.get('spin', 0)
        
        if not pdb_string:
            raise ValueError("No PDB string provided")
        
        lines = pdb_string.strip().split('\n')
        coordinates = []
        
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    # Extract coordinates from PDB line
                    element = line[76:78].strip()
                    if not element:
                        # Fallback to atom name
                        element = line[12:16].strip()[0]
                    
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    coordinates.append([element, x, y, z])
                    
                except (ValueError, IndexError) as e:
                    logger.warning(f"Skipping invalid PDB line: {line[:20]}...")
                    continue
        
        if not coordinates:
            raise ValueError("No valid coordinates found in PDB")
        
        return self._build_from_coordinates({
            'coordinates': coordinates,
            'charge': charge,
            'spin': spin
        })
    
    def get_test_molecules(self) -> Dict[str, Dict[str, Any]]:
        """Get a dictionary of test molecules for validation."""
        return {
            'water': {
                'type': 'coordinates',
                'coordinates': [
                    ['O', 0.0, 0.0, 0.0],
                    ['H', 0.757, 0.586, 0.0],
                    ['H', -0.757, 0.586, 0.0]
                ],
                'charge': 0,
                'spin': 0
            },
            'methane': {
                'type': 'coordinates',
                'coordinates': [
                    ['C', 0.0, 0.0, 0.0],
                    ['H', 1.089, 0.0, 0.0],
                    ['H', -0.363, 1.027, 0.0],
                    ['H', -0.363, -0.513, 0.889],
                    ['H', -0.363, -0.513, -0.889]
                ],
                'charge': 0,
                'spin': 0
            },
            'hydrogen': {
                'type': 'coordinates',
                'coordinates': [
                    ['H', 0.0, 0.0, 0.0],
                    ['H', 0.74, 0.0, 0.0]
                ],
                'charge': 0,
                'spin': 0
            }
        }