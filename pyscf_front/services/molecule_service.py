"""
分子管理サービス
既存のMoleculeクラスとデータベースを統合
"""
from typing import List, Dict, Any, Optional
from loguru import logger

from pyscf_front.core.molecule import Molecule as CoreMolecule, Atom
from pyscf_front.database import (
    MoleculeRepository, InstanceRepository,
    Molecule as DBMolecule
)


class MoleculeService:
    """分子管理サービスクラス"""
    
    def __init__(self):
        self.molecule_repo = MoleculeRepository()
        self.instance_repo = InstanceRepository()
    
    def create_molecule_instance(self, name: str, molecule: CoreMolecule,
                               description: Optional[str] = None) -> str:
        """分子インスタンスを作成してデータベースに保存"""
        try:
            # インスタンス作成
            with self.instance_repo as repo:
                instance = repo.create(
                    name=name,
                    description=description or f"Molecule: {molecule.name}"
                )
                
                # 分子データをデータベース形式に変換
                geometry_data = self._convert_molecule_to_db_format(molecule)
                
                # 分子をデータベースに保存
                with MoleculeRepository(repo.session) as mol_repo:
                    db_molecule = mol_repo.create(
                        instance_id=str(instance.id),  # type: ignore
                        name=molecule.name,
                        formula=molecule.formula if hasattr(molecule, 'formula') else None,
                        geometry_data=geometry_data,
                        charge=molecule.charge,
                        multiplicity=molecule.multiplicity
                    )
                
                logger.info(f"Created molecule instance: {instance.id}")
                return str(instance.id)
                
        except Exception as e:
            logger.error(f"Failed to create molecule instance: {e}")
            raise
    
    def load_molecule_from_instance(self, instance_id: str) -> Optional[CoreMolecule]:
        """インスタンスIDから分子オブジェクトを復元"""
        try:
            with self.molecule_repo as repo:
                molecules = repo.get_by_instance(instance_id)
                if not molecules:
                    return None
                
                # 最初の分子を取得（通常は1つのインスタンスに1つの分子）
                db_molecule = molecules[0]
                
                # データベース形式からコア分子クラスに変換
                core_molecule = self._convert_db_to_molecule_format(db_molecule)
                
                logger.info(f"Loaded molecule from instance: {instance_id}")
                return core_molecule
                
        except Exception as e:
            logger.error(f"Failed to load molecule from instance {instance_id}: {e}")
            return None
    
    def update_molecule_geometry(self, instance_id: str, molecule: CoreMolecule) -> bool:
        """分子の形状データを更新"""
        try:
            with self.molecule_repo as repo:
                molecules = repo.get_by_instance(instance_id)
                if not molecules:
                    return False
                
                db_molecule = molecules[0]
                geometry_data = self._convert_molecule_to_db_format(molecule)
                
                success = repo.update_geometry(str(db_molecule.id), geometry_data)
                if success:
                    logger.info(f"Updated molecule geometry for instance: {instance_id}")
                
                return success
                
        except Exception as e:
            logger.error(f"Failed to update molecule geometry: {e}")
            return False
    
    def get_molecule_instances(self, user_id: Optional[str] = None) -> List[Dict[str, Any]]:
        """分子インスタンスの一覧を取得"""
        try:
            with self.instance_repo as repo:
                instances = repo.get_all(user_id=user_id)
                
                result = []
                for instance in instances:
                    # 分子情報を含むインスタンス詳細を取得
                    molecules = self.molecule_repo.get_by_instance(str(instance.id))
                    molecule_info = molecules[0] if molecules else None
                    
                    result.append({
                        'instance_id': instance.id,
                        'name': instance.name,
                        'description': instance.description,
                        'status': instance.status.value,
                        'created_at': instance.created_at,
                        'molecule': {
                            'name': molecule_info.name if molecule_info else None,
                            'formula': molecule_info.formula if molecule_info else None,
                            'charge': molecule_info.charge if molecule_info else 0,
                            'multiplicity': molecule_info.multiplicity if molecule_info else 1
                        } if molecule_info else None
                    })
                
                return result
                
        except Exception as e:
            logger.error(f"Failed to get molecule instances: {e}")
            return []
    
    def delete_molecule_instance(self, instance_id: str) -> bool:
        """分子インスタンスを削除"""
        try:
            with self.instance_repo as repo:
                success = repo.delete(instance_id)
                if success:
                    logger.info(f"Deleted molecule instance: {instance_id}")
                return success
                
        except Exception as e:
            logger.error(f"Failed to delete molecule instance {instance_id}: {e}")
            return False
    
    def _convert_molecule_to_db_format(self, molecule: CoreMolecule) -> Dict[str, Any]:
        """CoreMoleculeオブジェクトをデータベース形式に変換"""
        atoms_data = []
        for atom in molecule.atoms:
            atoms_data.append({
                'symbol': atom.symbol,
                'x': atom.x,
                'y': atom.y,
                'z': atom.z,
                'charge': atom.charge if hasattr(atom, 'charge') else 0
            })
        
        return {
            'geometry_type': 'xyz',
            'atoms': atoms_data,
            'bonds': getattr(molecule, 'bonds', []),
            'properties': getattr(molecule, 'properties', {})
        }
    
    def _convert_db_to_molecule_format(self, db_molecule: DBMolecule) -> CoreMolecule:
        """データベース分子をCoreMoleculeオブジェクトに変換"""
        molecule = CoreMolecule(
            name=str(db_molecule.name) if db_molecule.name else "Untitled",  # type: ignore
            charge=int(db_molecule.charge) if db_molecule.charge else 0,  # type: ignore
            multiplicity=int(db_molecule.multiplicity) if db_molecule.multiplicity else 1  # type: ignore
        )
        
        # 形状データから原子を復元
        if db_molecule.geometry_data and isinstance(db_molecule.geometry_data, dict) and 'atoms' in db_molecule.geometry_data:  # type: ignore
            for atom_data in db_molecule.geometry_data['atoms']:
                atom = Atom(
                    symbol=str(atom_data['symbol']),
                    x=float(atom_data['x']),  # type: ignore
                    y=float(atom_data['y']),  # type: ignore
                    z=float(atom_data['z']),  # type: ignore
                    charge=int(atom_data.get('charge', 0))
                )
                molecule.atoms.append(atom)
        
        return molecule
    
    def create_from_smiles(self, smiles: str, name: Optional[str] = None) -> Optional[str]:
        """SMILES文字列から分子インスタンスを作成"""
        try:
            # CoreMoleculeクラスのcreate_from_smilesメソッドを使用
            molecule = CoreMolecule.create_from_smiles(smiles, name or f"SMILES_{smiles[:20]}")
            if molecule:
                instance_name = name or f"SMILES_{smiles[:20]}"
                return self.create_molecule_instance(instance_name, molecule)
            return None
            
        except Exception as e:
            logger.error(f"Failed to create molecule from SMILES {smiles}: {e}")
            return None
    
    def create_from_xyz_file(self, file_path: str, name: Optional[str] = None) -> Optional[str]:
        """XYZファイルから分子インスタンスを作成"""
        try:
            # CoreMoleculeクラスのcreate_from_xyzメソッドを使用
            default_name = f"XYZ_{file_path.split('/')[-1]}"
            molecule = CoreMolecule.create_from_xyz(file_path, name or default_name)
            if molecule:
                instance_name = name or default_name
                return self.create_molecule_instance(instance_name, molecule)
            return None
            
        except Exception as e:
            logger.error(f"Failed to create molecule from XYZ file {file_path}: {e}")
            return None