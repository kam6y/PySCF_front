# src/python/quantum_calc/file_manager.py

"""File management utilities for quantum chemistry calculations."""

import os
import shutil
import json
import re
from pathlib import Path
from typing import Optional, Dict, Any, List
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

    def rename_calculation(self, calculation_id: str, new_name: str) -> Optional[str]:
        """Updates the display name of a calculation without changing the directory name."""
        calc_path = self.base_dir / calculation_id
        if not calc_path.is_dir():
            return None

        # Read existing parameters
        params = self.read_calculation_parameters(str(calc_path))
        if params is None:
            params = {}

        # Update the display name
        params['name'] = new_name
        
        # Save updated parameters
        self.save_calculation_parameters(str(calc_path), params)
        
        # Return the same ID
        return calculation_id

    def _get_display_name(self, directory_name: str, parameters: Optional[Dict[str, Any]]) -> str:
        """Get display name from parameters."""
        if parameters and parameters.get('name'):
            return parameters['name'].strip()
        
        return "Unnamed Calculation"

    def list_calculations(self) -> list:
        """List all calculation directories."""
        if not self.base_dir.exists():
            return []
        
        calculations = []
        for item in self.base_dir.iterdir():
            if item.is_dir():
                params = self.read_calculation_parameters(str(item))
                # Use robust fallback logic for display names
                display_name = self._get_display_name(item.name, params)
                status = self.read_calculation_status(str(item))
                
                calculations.append({
                    'id': item.name, # Use directory name as the unique ID
                    'name': display_name, # Use this for display purposes
                    'path': str(item),
                    'date': datetime.fromtimestamp(item.stat().st_mtime).isoformat(),
                    'has_checkpoint': (item / "calculation.chk").exists(),
                    'status': status
                })
        
        return sorted(calculations, key=lambda x: x['date'], reverse=True)

    def save_geometry(self, calc_dir: str, xyz_string: str, filename: str = "optimized_geometry.xyz"):
        """Save geometry to an XYZ file."""
        if not xyz_string:
            return
        geom_file = Path(calc_dir) / filename
        with open(geom_file, 'w') as f:
            f.write(xyz_string)

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
            
    def save_calculation_status(self, calc_dir: str, status: str, waiting_reason: Optional[str] = None) -> None:
        """Save calculation status to JSON file."""
        status_file = Path(calc_dir) / "status.json"
        status_data = {
            'status': status,
            'updated_at': datetime.now().isoformat()
        }
        if waiting_reason is not None:
            status_data['waiting_reason'] = waiting_reason
        
        with open(status_file, 'w') as f:
            json.dump(status_data, f, indent=2)

    def read_calculation_status(self, calc_dir: str) -> str:
        """Read calculation status from JSON file."""
        status_file = Path(calc_dir) / "status.json"
        if not status_file.exists():
            if self.file_exists(calc_dir, 'results.json'):
                return 'completed'
            if self.file_exists(calc_dir, 'parameters.json'):
                 return 'pending'
            return 'error' 
        try:
            with open(status_file, 'r') as f:
                return json.load(f).get('status', 'pending')
        except (json.JSONDecodeError, OSError):
            return 'pending'
    
    def read_calculation_status_details(self, calc_dir: str) -> tuple[str, Optional[str]]:
        """Read calculation status and waiting reason from JSON file."""
        status_file = Path(calc_dir) / "status.json"
        if not status_file.exists():
            if self.file_exists(calc_dir, 'results.json'):
                return 'completed', None
            if self.file_exists(calc_dir, 'parameters.json'):
                return 'pending', None
            return 'error', None
        
        try:
            with open(status_file, 'r') as f:
                status_data = json.load(f)
                status = status_data.get('status', 'pending')
                waiting_reason = status_data.get('waiting_reason')
                return status, waiting_reason
        except (json.JSONDecodeError, OSError):
            return 'pending', None

    def get_cube_files_info(self, calc_dir: str) -> List[Dict[str, Any]]:
        """Get information about CUBE files in a calculation directory."""
        cube_files = []
        calc_path = Path(calc_dir)
        orbital_dir = calc_path / "orbital"
        
        if not orbital_dir.exists():
            return cube_files
        
        # Find all cube files matching our naming pattern
        import glob
        pattern = str(orbital_dir / "orbital_*_grid*.cube")
        
        for file_path in glob.glob(pattern):
            filename = os.path.basename(file_path)
            
            try:
                # Parse filename: orbital_{index}_grid{size}.cube
                parts = filename.replace('.cube', '').split('_')
                if len(parts) >= 3 and parts[0] == 'orbital' and parts[2].startswith('grid'):
                    orbital_index = int(parts[1])
                    grid_size = int(parts[2].replace('grid', ''))
                    
                    file_size_kb = os.path.getsize(file_path) / 1024.0
                    modified_time = datetime.fromtimestamp(os.path.getmtime(file_path))
                    
                    cube_files.append({
                        "filename": filename,
                        "file_path": file_path,
                        "orbital_index": orbital_index,
                        "grid_size": grid_size,
                        "file_size_kb": file_size_kb,
                        "modified": modified_time.isoformat()
                    })
            except (ValueError, IndexError):
                continue
        
        return sorted(cube_files, key=lambda x: x["orbital_index"])

    def delete_cube_files(self, calc_dir: str, orbital_index: int = None) -> int:
        """
        Delete CUBE files from a calculation directory.
        
        Args:
            calc_dir: Calculation directory path
            orbital_index: Specific orbital index to delete, or None to delete all
        
        Returns:
            Number of files deleted
        """
        deleted_count = 0
        calc_path = Path(calc_dir)
        orbital_dir = calc_path / "orbital"
        
        if not orbital_dir.exists():
            return deleted_count
        
        import glob
        if orbital_index is not None:
            pattern = str(orbital_dir / f"orbital_{orbital_index}_grid*.cube")
        else:
            pattern = str(orbital_dir / "orbital_*_grid*.cube")
        
        for file_path in glob.glob(pattern):
            try:
                os.unlink(file_path)
                deleted_count += 1
            except Exception:
                continue
        
        return deleted_count

    def cleanup_calculation_directory(self, calc_dir: str, include_cube_files: bool = True) -> Dict[str, int]:
        """
        Clean up files in a calculation directory.
        
        Args:
            calc_dir: Calculation directory path
            include_cube_files: Whether to delete CUBE files
        
        Returns:
            Dictionary with counts of deleted file types
        """
        cleanup_stats = {
            "cube_files": 0,
            "temp_files": 0,
            "log_files": 0
        }
        
        calc_path = Path(calc_dir)
        if not calc_path.exists():
            return cleanup_stats
        
        # Delete CUBE files
        if include_cube_files:
            cleanup_stats["cube_files"] = self.delete_cube_files(calc_dir)
        
        # Delete temporary files
        import glob
        temp_patterns = [
            str(calc_path / "*.tmp"),
            str(calc_path / "temp_*"),
            str(calc_path / "*.temp")
        ]
        
        for pattern in temp_patterns:
            for file_path in glob.glob(pattern):
                try:
                    os.unlink(file_path)
                    cleanup_stats["temp_files"] += 1
                except Exception:
                    continue
        
        # Delete log files (optional, keep recent ones)
        log_pattern = str(calc_path / "*.log")
        log_files = glob.glob(log_pattern)
        
        # Keep only the most recent 3 log files
        if len(log_files) > 3:
            log_files.sort(key=lambda f: os.path.getmtime(f), reverse=True)
            for old_log in log_files[3:]:
                try:
                    os.unlink(old_log)
                    cleanup_stats["log_files"] += 1
                except Exception:
                    continue
        
        return cleanup_stats