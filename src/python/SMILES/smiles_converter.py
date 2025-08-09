"""SMILES to XYZ converter using RDKit."""
import logging
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

class SMILESError(Exception):
    """Custom exception for SMILES conversion errors."""
    pass

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
        # 1. SMILESからMoleculeオブジェクトを生成
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            raise SMILESError(f"Invalid SMILES string: {smiles_string}")

        # 2. 水素原子を付加
        mol = Chem.AddHs(mol)

        # 3. 3D構造を生成 (EmbedMolecule)
        # randomSeedを指定して再現性を確保
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
             raise SMILESError("Failed to generate 3D coordinates. The structure may be too complex or constrained.")

        # 4. 構造最適化 (MMFF94力場を使用)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # 5. XYZ形式の文字列を構築
        conformer = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        xyz_lines = [str(num_atoms), title]
        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            element = atom.GetSymbol()
            xyz_lines.append(f"{element:<3} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f}")
            
        return "\n".join(xyz_lines)

    except SMILESError:
        # SMILESErrorはそのまま上位に投げる
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred during SMILES conversion: {e}")
        raise SMILESError("An internal error occurred during SMILES conversion.")