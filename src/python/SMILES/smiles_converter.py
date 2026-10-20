"""SMILES to XYZ converter using RDKit."""
import logging
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

class SMILESError(Exception):
    """Custom exception for SMILES conversion errors."""
    pass


def _format_xyz_atom_line(element: str, x: float, y: float, z: float) -> str:
    """Format one atom line in XYZ coordinate format."""
    return f"{element:<3} {x:>7.4f} {y:>7.4f} {z:>7.4f}"


def smiles_to_xyz(smiles_string: str, title: str = "Molecule from SMILES") -> str:
    """
    Converts a SMILES string to a 3D structure in XYZ format.
    
    Args:
        smiles_string: The SMILES string of the molecule.
        title: A title for the second line of the XYZ string.

    Returns:
        A string in XYZ format.
        
    Raises:
        SMILESError: If the SMILES is invalid or 3D embedding fails.
    """
    try:
        # Create the RDKit molecule from SMILES.
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            raise SMILESError("Invalid SMILES string")

        # Add explicit hydrogens before 3D embedding.
        mol = Chem.AddHs(mol)

        # Use a fixed seed so generated coordinates are reproducible.
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
            raise SMILESError("Failed to generate 3D coordinates. The structure may be too complex or constrained.")

        # Optimize with the MMFF94 force field.
        AllChem.MMFFOptimizeMolecule(mol)
        
        conformer = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        xyz_lines = [str(num_atoms), title]
        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            element = atom.GetSymbol()
            xyz_lines.append(_format_xyz_atom_line(element, pos.x, pos.y, pos.z))

        return "\n".join(xyz_lines)

    except SMILESError:
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred during SMILES conversion: {e}")
        raise SMILESError("An internal error occurred during SMILES conversion.") from None
