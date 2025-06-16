"""
インスタンス管理サービス
分子と計算を統合管理
"""
from typing import List, Dict, Any, Optional
from loguru import logger

from pyscf_front.database import (
    InstanceRepository, MoleculeRepository, InstanceStatus
)
from .molecule_service import MoleculeService
from .calculation_service import CalculationService


class InstanceService:
    """インスタンス統合管理サービス"""
    
    def __init__(self):
        self.instance_repo = InstanceRepository()
        self.molecule_service = MoleculeService()
        self.calculation_service = CalculationService()
    
    def create_complete_instance(self, name: str, description: Optional[str] = None,
                                molecule_data: Optional[Dict[str, Any]] = None,
                                calculation_config: Optional[Dict[str, Any]] = None) -> str:
        """分子と計算設定を含む完全なインスタンスを作成"""
        try:
            # 分子インスタンスを作成
            if molecule_data:
                # 分子オブジェクトから作成
                if 'molecule_object' in molecule_data:
                    instance_id = self.molecule_service.create_molecule_instance(
                        name=name,
                        molecule=molecule_data['molecule_object'],
                        description=description
                    )
                # SMILES文字列から作成
                elif 'smiles' in molecule_data:
                    instance_id = self.molecule_service.create_from_smiles(
                        smiles=molecule_data['smiles'],
                        name=name
                    )
                # XYZファイルから作成
                elif 'xyz_file' in molecule_data:
                    instance_id = self.molecule_service.create_from_xyz_file(
                        file_path=molecule_data['xyz_file'],
                        name=name
                    )
                else:
                    raise ValueError("Invalid molecule_data format")
            else:
                # 分子なしでインスタンスのみ作成
                with self.instance_repo as repo:
                    instance = repo.create(name=name, description=description)
                    instance_id = instance.id
            
            # 計算設定を追加
            if calculation_config and instance_id:  # type: ignore
                self.calculation_service.create_calculation(
                    instance_id=str(instance_id),
                    method=calculation_config['method'],
                    basis_set=calculation_config['basis_set'],
                    parameters=calculation_config.get('parameters'),
                    priority=calculation_config.get('priority', 5)
                )
            
            logger.info(f"Created complete instance: {instance_id}")
            return str(instance_id)
            
        except Exception as e:
            logger.error(f"Failed to create complete instance: {e}")
            raise
    
    def get_instance_details(self, instance_id: str) -> Optional[Dict[str, Any]]:
        """インスタンスの詳細情報を取得"""
        try:
            with self.instance_repo as repo:
                instance = repo.get_by_id(instance_id)
                if not instance:
                    return None
                
                # 基本情報
                details = {
                    'instance_id': instance.id,
                    'name': instance.name,
                    'description': instance.description,
                    'status': instance.status.value,
                    'created_at': instance.created_at,
                    'updated_at': instance.updated_at,
                    'user_id': instance.user_id,
                    'project_id': instance.project_id
                }
                
                # 分子情報を追加
                molecule = self.molecule_service.load_molecule_from_instance(instance_id)
                if molecule:
                    details['molecule'] = {
                        'name': molecule.name,
                        'formula': getattr(molecule, 'formula', None),
                        'charge': molecule.charge,
                        'multiplicity': molecule.multiplicity,
                        'atom_count': len(molecule.atoms),
                        'atoms': [
                            {
                                'symbol': atom.symbol,
                                'coordinates': [atom.x, atom.y, atom.z]
                            } for atom in molecule.atoms
                        ]
                    }
                
                # 計算情報を追加
                calculations = self.calculation_service.get_calculations_by_instance(instance_id)
                details['calculations'] = calculations
                
                # 統計情報を追加
                details['statistics'] = {
                    'total_calculations': len(calculations),
                    'completed_calculations': len([c for c in calculations if c['status'] == 'completed']),
                    'running_calculations': len([c for c in calculations if c['status'] == 'running']),
                    'failed_calculations': len([c for c in calculations if c['status'] == 'failed'])
                }
                
                return details
                
        except Exception as e:
            logger.error(f"Failed to get instance details {instance_id}: {e}")
            return None
    
    def get_all_instances(self, user_id: Optional[str] = None,
                         include_details: bool = False) -> List[Dict[str, Any]]:
        """すべてのインスタンスを取得"""
        try:
            with self.instance_repo as repo:
                instances = repo.get_all(user_id=user_id)
                
                result = []
                for instance in instances:
                    if include_details:
                        instance_data = self.get_instance_details(str(instance.id))
                    else:
                        # 基本情報のみ
                        instance_data = {
                            'instance_id': instance.id,
                            'name': instance.name,
                            'description': instance.description,
                            'status': instance.status.value,
                            'created_at': instance.created_at,
                            'updated_at': instance.updated_at
                        }
                        
                        # 分子の基本情報を追加
                        with MoleculeRepository(repo.session) as mol_repo:
                            molecules = mol_repo.get_by_instance(str(instance.id))
                        if molecules:
                            mol = molecules[0]
                            instance_data['molecule_name'] = mol.name
                            instance_data['molecule_formula'] = mol.formula
                    
                    if instance_data:
                        result.append(instance_data)
                
                return result
                
        except Exception as e:
            logger.error(f"Failed to get all instances: {e}")
            return []
    
    def update_instance_status(self, instance_id: str, status: str) -> bool:
        """インスタンススステータスを更新"""
        try:
            status_enum = InstanceStatus(status)
            with self.instance_repo as repo:
                success = repo.update_status(instance_id, status_enum)
                if success:
                    logger.info(f"Updated instance {instance_id} status to {status}")
                return success
                
        except ValueError:
            logger.error(f"Invalid status: {status}")
            return False
        except Exception as e:
            logger.error(f"Failed to update instance status {instance_id}: {e}")
            return False
    
    def delete_instance(self, instance_id: str) -> bool:
        """インスタンスを削除（分子と計算も含む）"""
        try:
            with self.instance_repo as repo:
                success = repo.delete(instance_id)
                if success:
                    logger.info(f"Deleted instance: {instance_id}")
                return success
                
        except Exception as e:
            logger.error(f"Failed to delete instance {instance_id}: {e}")
            return False
    
    def run_calculation(self, instance_id: str, method: Optional[str] = None,
                       basis_set: Optional[str] = None) -> Optional[str]:
        """インスタンスの計算を実行"""
        try:
            # インスタンスの詳細を取得
            details = self.get_instance_details(instance_id)
            if not details:
                logger.error(f"Instance not found: {instance_id}")
                return None
            
            # 分子が存在するかチェック
            if not details.get('molecule'):
                logger.error(f"No molecule found for instance: {instance_id}")
                return None
            
            # 計算設定を決定
            if not method or not basis_set:
                # デフォルト値を使用
                method = method or "B3LYP"
                basis_set = basis_set or "6-31G"
            
            # 計算を作成
            calculation_id = self.calculation_service.create_calculation(
                instance_id=instance_id,
                method=str(method),
                basis_set=str(basis_set)
            )
            
            # インスタンスステータスを更新
            self.update_instance_status(instance_id, "running")
            
            logger.info(f"Started calculation for instance {instance_id}: {calculation_id}")
            return calculation_id
            
        except Exception as e:
            logger.error(f"Failed to run calculation for instance {instance_id}: {e}")
            return None
    
    def get_instance_results(self, instance_id: str) -> Dict[str, Any]:
        """インスタンスのすべての計算結果を取得"""
        try:
            results = {}
            calculations = self.calculation_service.get_calculations_by_instance(instance_id)
            
            for calc in calculations:
                calc_results = self.calculation_service.get_calculation_results(calc['id'])
                if calc_results:
                    results[calc['id']] = {
                        'calculation_info': calc,
                        'results': calc_results
                    }
            
            return results
            
        except Exception as e:
            logger.error(f"Failed to get instance results {instance_id}: {e}")
            return {}
    
    def search_instances(self, query: str, user_id: Optional[str] = None) -> List[Dict[str, Any]]:
        """インスタンスを検索"""
        try:
            # 簡単な検索実装（名前と説明で検索）
            all_instances = self.get_all_instances(user_id=user_id)
            
            filtered_instances = []
            query_lower = query.lower()
            
            for instance in all_instances:
                # 名前、説明、分子名で検索（None値対応）
                name = instance.get('name') or ''
                description = instance.get('description') or ''
                molecule_name = instance.get('molecule_name') or ''
                
                if (query_lower in name.lower() or
                    query_lower in description.lower() or
                    query_lower in molecule_name.lower()):
                    filtered_instances.append(instance)
            
            return filtered_instances
            
        except Exception as e:
            logger.error(f"Failed to search instances with query '{query}': {e}")
            return []