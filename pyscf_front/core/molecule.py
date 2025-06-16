"""
分子データ管理モジュール
"""
import uuid
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import numpy as np
from loguru import logger

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available. Some molecular features will be limited.")
    RDKIT_AVAILABLE = False
    Chem = None
    AllChem = None


@dataclass
class Atom:
    """原子データクラス"""
    symbol: str
    x: float
    y: float
    z: float
    charge: int = 0
    
    @property
    def coordinates(self) -> Tuple[float, float, float]:
        """座標タプルを返す"""
        return (self.x, self.y, self.z)
    
    @property
    def atomic_number(self) -> int:
        """原子番号を返す"""
        atomic_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
            'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
            'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20
        }
        return atomic_numbers.get(self.symbol, 0)


class Molecule:
    """分子クラス"""
    
    def __init__(self, name: str = "Untitled", charge: int = 0, multiplicity: int = 1):
        self.id = str(uuid.uuid4())
        self.name = name
        self.atoms: List[Atom] = []
        self.charge = charge
        self.multiplicity = multiplicity
        self.bonds: List[Tuple[int, int, int]] = []  # (atom1_idx, atom2_idx, bond_order)
        self.properties: Dict[str, Any] = {}
        
    def add_atom(self, atom: Atom) -> int:
        """原子を追加し、インデックスを返す"""
        self.atoms.append(atom)
        return len(self.atoms) - 1
    
    def add_bond(self, atom1_idx: int, atom2_idx: int, bond_order: int = 1):
        """結合を追加"""
        if 0 <= atom1_idx < len(self.atoms) and 0 <= atom2_idx < len(self.atoms):
            self.bonds.append((atom1_idx, atom2_idx, bond_order))
        else:
            raise ValueError(f"Invalid atom indices: {atom1_idx}, {atom2_idx}")
    
    def get_coordinates(self) -> np.ndarray:
        """座標行列を取得 (N x 3)"""
        if not self.atoms:
            return np.empty((0, 3))
        return np.array([atom.coordinates for atom in self.atoms])
    
    def get_atomic_numbers(self) -> np.ndarray:
        """原子番号配列を取得"""
        return np.array([atom.atomic_number for atom in self.atoms])
    
    def get_symbols(self) -> List[str]:
        """原子記号リストを取得"""
        return [atom.symbol for atom in self.atoms]
    
    def to_xyz(self) -> str:
        """XYZ形式の文字列に変換"""
        lines = [str(len(self.atoms))]
        lines.append(f"{self.name} (charge={self.charge}, mult={self.multiplicity})")
        
        for atom in self.atoms:
            lines.append(f"{atom.symbol:2s} {atom.x:12.6f} {atom.y:12.6f} {atom.z:12.6f}")
        
        return "\n".join(lines)
    
    def from_xyz_string(self, xyz_string: str):
        """XYZ形式の文字列から分子を構築"""
        lines = xyz_string.strip().split('\n')
        if len(lines) < 2:
            raise ValueError("Invalid XYZ format")
        
        try:
            num_atoms = int(lines[0])
            self.name = lines[1] if len(lines) > 1 else "Unknown"
            
            self.atoms.clear()
            for i in range(2, 2 + num_atoms):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 4:
                    symbol = parts[0]
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    self.add_atom(Atom(symbol, x, y, z))
                    
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error parsing XYZ string: {e}")
    
    @property
    def formula(self) -> str:
        """分子式プロパティ"""
        return self.get_molecular_formula()
    
    def get_molecular_formula(self) -> str:
        """分子式を取得"""
        element_count = {}
        for atom in self.atoms:
            element_count[atom.symbol] = element_count.get(atom.symbol, 0) + 1
        
        # 炭素、水素、その他の順で並べる
        formula_parts = []
        for element in ['C', 'H']:
            if element in element_count:
                count = element_count[element]
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
                del element_count[element]
        
        # 残りの元素をアルファベット順で追加
        for element in sorted(element_count.keys()):
            count = element_count[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
        
        return "".join(formula_parts) if formula_parts else ""
    
    def calculate_center_of_mass(self) -> np.ndarray:
        """重心を計算"""
        if not self.atoms:
            return np.array([0.0, 0.0, 0.0])
        
        atomic_masses = {
            'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
            'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
            'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
            'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078
        }
        
        total_mass = 0.0
        center = np.array([0.0, 0.0, 0.0])
        
        for atom in self.atoms:
            mass = atomic_masses.get(atom.symbol, 1.0)
            total_mass += mass
            center += mass * np.array(atom.coordinates)
        
        return center / total_mass if total_mass > 0 else center
    
    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """分子の境界ボックスを取得"""
        if not self.atoms:
            return np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])
        
        coords = self.get_coordinates()
        return np.min(coords, axis=0), np.max(coords, axis=0)
    
    @staticmethod
    def create_from_smiles(smiles: str, name: Optional[str] = None) -> Optional['Molecule']:
        """SMILES文字列から分子を作成"""
        return MoleculeBuilder.from_smiles(smiles, name)
    
    @staticmethod
    def create_from_xyz(file_path: str, name: Optional[str] = None) -> Optional['Molecule']:
        """XYZファイルから分子を作成"""
        try:
            with open(file_path, 'r') as f:
                xyz_content = f.read()
            
            molecule = Molecule(name or f"Molecule_from_{file_path}")
            molecule.from_xyz_string(xyz_content)
            return molecule
        except Exception as e:
            logger.error(f"Error creating molecule from XYZ file {file_path}: {e}")
            return None


class MoleculeBuilder:
    """分子構築ヘルパークラス"""
    
    @staticmethod
    def from_smiles(smiles: str, name: Optional[str] = None) -> Optional[Molecule]:
        """SMILES文字列から分子を構築"""
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is required for SMILES parsing")
            return None
        
        try:
            if not RDKIT_AVAILABLE or Chem is None:
                logger.error("RDKit is not available")
                return None
                
            mol = Chem.MolFromSmiles(smiles)  # type: ignore
            if mol is None:
                logger.error(f"Invalid SMILES string: {smiles}")
                return None
            
            # 3D座標を生成
            mol = Chem.AddHs(mol)  # type: ignore
            if AllChem is not None:
                try:
                    AllChem.EmbedMolecule(mol, randomSeed=42)  # type: ignore
                    AllChem.UFFOptimizeMolecule(mol)  # type: ignore
                except Exception as e:
                    logger.warning(f"Could not optimize molecule geometry: {e}")
            
            # PySCF_Front分子オブジェクトに変換
            molecule = Molecule(name or f"Molecule_from_{smiles}")
            
            conf = mol.GetConformer()
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                symbol = atom.GetSymbol()
                molecule.add_atom(Atom(symbol, pos.x, pos.y, pos.z))
            
            # 結合情報を追加
            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                bond_order = int(bond.GetBondType())
                molecule.add_bond(atom1_idx, atom2_idx, bond_order)
            
            return molecule
            
        except Exception as e:
            logger.error(f"Error creating molecule from SMILES {smiles}: {e}")
            return None
    
    @staticmethod
    def create_water() -> Molecule:
        """水分子を作成"""
        mol = Molecule("Water", charge=0, multiplicity=1)
        mol.add_atom(Atom("O", 0.0, 0.0, 0.0))
        mol.add_atom(Atom("H", 0.757, 0.586, 0.0))
        mol.add_atom(Atom("H", -0.757, 0.586, 0.0))
        mol.add_bond(0, 1, 1)  # O-H
        mol.add_bond(0, 2, 1)  # O-H
        return mol
    
    @staticmethod
    def create_methane() -> Molecule:
        """メタン分子を作成"""
        mol = Molecule("Methane", charge=0, multiplicity=1)
        mol.add_atom(Atom("C", 0.0, 0.0, 0.0))
        mol.add_atom(Atom("H", 1.089, 0.0, 0.0))
        mol.add_atom(Atom("H", -0.363, 1.027, 0.0))
        mol.add_atom(Atom("H", -0.363, -0.513, 0.889))
        mol.add_atom(Atom("H", -0.363, -0.513, -0.889))
        
        # C-H結合を追加
        for i in range(1, 5):
            mol.add_bond(0, i, 1)
        
        return mol
    
    @staticmethod
    def create_benzene() -> Molecule:
        """ベンゼン分子を作成"""
        mol = Molecule("Benzene", charge=0, multiplicity=1)
        
        # 炭素原子の配置（正六角形）
        import math
        radius = 1.397  # C-C結合長 (Å)
        
        for i in range(6):
            angle = i * math.pi / 3
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            mol.add_atom(Atom("C", x, y, 0.0))
        
        # 水素原子の配置
        h_radius = radius + 1.084  # C-H結合長
        for i in range(6):
            angle = i * math.pi / 3
            x = h_radius * math.cos(angle)
            y = h_radius * math.sin(angle)
            mol.add_atom(Atom("H", x, y, 0.0))
        
        # C-C結合（芳香環）
        for i in range(6):
            mol.add_bond(i, (i + 1) % 6, 1)  # 単結合として簡略化
        
        # C-H結合
        for i in range(6):
            mol.add_bond(i, i + 6, 1)
        
        return mol


class MoleculeManager:
    """分子管理クラス"""
    
    def __init__(self):
        self.molecules: Dict[str, Molecule] = {}
        
    def add_molecule(self, molecule: Molecule) -> str:
        """分子を追加し、IDを返す"""
        self.molecules[molecule.id] = molecule
        logger.info(f"Added molecule: {molecule.name} (ID: {molecule.id})")
        return molecule.id
    
    def get_molecule(self, molecule_id: str) -> Optional[Molecule]:
        """分子を取得"""
        return self.molecules.get(molecule_id)
    
    def remove_molecule(self, molecule_id: str) -> bool:
        """分子を削除"""
        if molecule_id in self.molecules:
            del self.molecules[molecule_id]
            logger.info(f"Removed molecule ID: {molecule_id}")
            return True
        return False
    
    def list_molecules(self) -> List[Dict[str, Any]]:
        """分子リストを取得"""
        return [
            {
                'id': mol.id,
                'name': mol.name,
                'formula': mol.get_molecular_formula(),
                'atoms': len(mol.atoms),
                'charge': mol.charge,
                'multiplicity': mol.multiplicity
            }
            for mol in self.molecules.values()
        ]
    
    def create_from_smiles(self, smiles: str, name: Optional[str] = None) -> Optional[str]:
        """SMILES文字列から分子を作成"""
        molecule = MoleculeBuilder.from_smiles(smiles, name)
        if molecule:
            return self.add_molecule(molecule)
        return None
    
    def create_from_xyz(self, xyz_string: str, name: Optional[str] = None) -> Optional[str]:
        """XYZ文字列から分子を作成"""
        try:
            molecule = Molecule(name or "XYZ_Molecule")
            molecule.from_xyz_string(xyz_string)
            return self.add_molecule(molecule)
        except ValueError as e:
            logger.error(f"Error creating molecule from XYZ: {e}")
            return None