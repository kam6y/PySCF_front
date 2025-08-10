# src/python/quantum_calc/file_manager.py

"""File management utilities for quantum chemistry calculations."""

import os
import shutil
import json
import re
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime


class CalculationFileManager:
    """Manages files generated during quantum chemistry calculations."""
    
    def __init__(self, base_dir: Optional[str] = None):
        if base_dir is None:
            home = Path.home()
            self.base_dir = home / "PySCF_Calculations"
        else:
            self.base_dir = Path(base_dir)
        
        self.base_dir.mkdir(exist_ok=True)
    
    def create_calculation_dir(self, molecule_name: Optional[str] = None) -> str:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if molecule_name:
            clean_name = "".join(c for c in molecule_name if c.isalnum() or c in "._-").strip() or "unnamed"
            dir_name = f"{clean_name}_{timestamp}"
        else:
            dir_name = f"calculation_{timestamp}"
        
        calc_dir = self.base_dir / dir_name
        calc_dir.mkdir(exist_ok=True)
        
        return str(calc_dir)

    # --- 新規追加: 計算の名前変更とJSON更新を行う関数 ---
    def rename_calculation(self, old_id: str, new_name: str) -> Optional[str]:
        """Renames a calculation directory and updates its metadata."""
        old_path = self.base_dir / old_id
        if not old_path.is_dir():
            return None

        # タイムスタンプ部分を維持しつつ、新しいディレクトリ名を生成
        timestamp_match = re.search(r'_(\d{8}_\d{6})$', old_id)
        timestamp = timestamp_match.group(1) if timestamp_match else datetime.now().strftime("%Y%m%d_%H%M%S")
        
        clean_new_name = "".join(c for c in new_name if c.isalnum() or c in "._-").strip()
        new_id = f"{clean_new_name}_{timestamp}"
        new_path = self.base_dir / new_id

        if old_path == new_path:
            return old_id # 名前が実質的に変わらない場合は何もしない

        if new_path.exists():
            raise FileExistsError(f"A calculation with the name '{new_id}' already exists.")

        # ディレクトリ名を変更
        os.rename(old_path, new_path)

        # parameters.json内のmolecule_nameを更新
        params = self.read_calculation_parameters(str(new_path))
        if params:
            params['molecule_name'] = new_name
            self.save_calculation_parameters(str(new_path), params)
            
        return new_id

    def list_calculations(self) -> list:
        """List all calculation directories."""
        if not self.base_dir.exists():
            return []
        
        calculations = []
        for item in self.base_dir.iterdir():
            if item.is_dir():
                # --- 変更点: JSONから名前を取得し、なければディレクトリ名から推測 ---
                params = self.read_calculation_parameters(str(item))
                display_name = params.get('molecule_name', item.name.replace(r'_(\d{8}_\d{6})$', '')) if params else item.name.replace(r'_(\d{8}_\d{6})$', '')

                status = self.read_calculation_status(str(item))
                
                calculations.append({
                    'name': item.name, # `id`としてディレクトリ名を使用
                    'displayName': display_name, # UI表示用の名前
                    'path': str(item),
                    'date': datetime.fromtimestamp(item.stat().st_mtime).isoformat(),
                    'has_checkpoint': (item / "calculation.chk").exists(),
                    'status': status
                })
        
        return sorted(calculations, key=lambda x: x['date'], reverse=True)

    # ... (以降の関数は変更なし) ...
    def get_base_directory(self) -> str:
        """Get the base directory path."""
        return str(self.base_dir)

    def file_exists(self, calc_dir: str, filename: str) -> bool:
        """Check if a specific file exists in the calculation directory."""
        return (Path(calc_dir) / filename).exists()

    def save_calculation_parameters(self, calc_dir: str, parameters: Dict[str, Any]) -> None:
        """Save calculation parameters to JSON file."""
        params_file = Path(calc_dir) / "parameters.json"
        with open(params_file, 'w') as f:
            json.dump(parameters, f, indent=2, default=str)

    def read_calculation_parameters(self, calc_dir: str) -> Optional[Dict[str, Any]]:
        """Read calculation parameters from JSON file."""
        params_file = Path(calc_dir) / "parameters.json"
        if not params_file.exists():
            return None
        try:
            with open(params_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            return None

    def save_calculation_results(self, calc_dir: str, results: Dict[str, Any]) -> None:
        """Save calculation results to JSON file."""
        results_file = Path(calc_dir) / "results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

    def read_calculation_results(self, calc_dir: str) -> Optional[Dict[str, Any]]:
        """Read calculation results from JSON file."""
        results_file = Path(calc_dir) / "results.json"
        if not results_file.exists():
            return None
        try:
            with open(results_file, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            return None
            
    def save_calculation_status(self, calc_dir: str, status: str) -> None:
        """Save calculation status to JSON file."""
        status_file = Path(calc_dir) / "status.json"
        status_data = {
            'status': status,
            'updated_at': datetime.now().isoformat()
        }
        with open(status_file, 'w') as f:
            json.dump(status_data, f, indent=2)

    def read_calculation_status(self, calc_dir: str) -> str:
        """Read calculation status from JSON file."""
        status_file = Path(calc_dir) / "status.json"
        if not status_file.exists():
            return 'completed' if self.file_exists(calc_dir, 'calculation.chk') else 'pending'
        try:
            with open(status_file, 'r') as f:
                return json.load(f).get('status', 'pending')
        except (json.JSONDecodeError, OSError):
            return 'pending'