# src/python/quantum_calc/file_manager.py

"""File management utilities for quantum chemistry calculations."""

import os
import shutil
import json
import re
import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

logger = logging.getLogger(__name__)


class CalculationFileManager:
    """Manages files generated during quantum chemistry calculations."""

    def __init__(self, base_dir: Optional[str] = None):
        if base_dir is None:
            # Default to ~/PySCF_calculations if no base_dir is provided
            home = Path.home()
            self.base_dir = home / "PySCF_calculations"
        else:
            # Use the provided path as-is (should already include PySCF_calculations)
            self.base_dir = Path(base_dir)

        self.base_dir.mkdir(parents=True, exist_ok=True)
    
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

    def list_calculations(
        self,
        name_query: Optional[str] = None,
        status: Optional[str] = None,
        calculation_method: Optional[str] = None,
        basis_function: Optional[str] = None,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None
    ) -> list:
        """
        List calculation directories with optional filtering.

        Args:
            name_query: Partial match search in calculation name (case-insensitive)
            status: Filter by status ("completed", "running", "error", etc.)
            calculation_method: Filter by calculation method ("DFT", "HF", "MP2", etc.)
            basis_function: Filter by basis set (case-insensitive)
            date_from: Start date for date range filtering (ISO format: YYYY-MM-DD)
            date_to: End date for date range filtering (ISO format: YYYY-MM-DD)

        Returns:
            List of calculation dictionaries matching the filter criteria, sorted by date (newest first)
        """
        if not self.base_dir.exists():
            return []

        # Helper function to parse ISO date strings
        def parse_date(date_str: str) -> Optional[datetime]:
            """Parse ISO date string to datetime object."""
            try:
                return datetime.fromisoformat(date_str.replace('Z', '+00:00'))
            except (ValueError, AttributeError):
                return None

        # Parse date filters once
        date_from_dt = parse_date(date_from) if date_from else None
        date_to_dt = parse_date(date_to) if date_to else None

        calculations = []
        for item in self.base_dir.iterdir():
            if not item.is_dir():
                continue

            params = self.read_calculation_parameters(str(item))
            display_name = self._get_display_name(item.name, params)
            calc_status = self.read_calculation_status(str(item))

            # Apply filters - skip this calculation if it doesn't match

            # Filter by name (case-insensitive partial match)
            if name_query:
                if name_query.lower() not in display_name.lower():
                    continue

            # Filter by status (exact match, case-insensitive)
            if status:
                if calc_status.lower() != status.lower():
                    continue

            # Filter by calculation_method (requires params)
            if calculation_method and params:
                calc_method = params.get('calculation_method', '')
                if calc_method.lower() != calculation_method.lower():
                    continue
            elif calculation_method:
                # No params available, can't match
                continue

            # Filter by basis_function (requires params)
            if basis_function and params:
                calc_basis = params.get('basis_function', '')
                if calc_basis.lower() != basis_function.lower():
                    continue
            elif basis_function:
                # No params available, can't match
                continue

            # Filter by date range (requires params with created_at)
            if (date_from_dt or date_to_dt) and params:
                created_at = params.get('created_at')
                if not created_at:
                    continue

                calc_date = parse_date(created_at)
                if not calc_date:
                    continue

                # Check date range
                if date_from_dt and calc_date < date_from_dt:
                    continue
                if date_to_dt and calc_date > date_to_dt:
                    continue
            elif (date_from_dt or date_to_dt):
                # No params or created_at available, can't match
                continue

            # All filters passed - add to results
            calculations.append({
                'id': item.name,  # Use directory name as the unique ID
                'name': display_name,  # Use this for display purposes
                'path': str(item),
                'date': datetime.fromtimestamp(item.stat().st_mtime).isoformat(),
                'has_checkpoint': (item / "calculation.chk").exists(),
                'status': calc_status,
                # Include additional fields for better search results
                'calculation_method': params.get('calculation_method') if params else None,
                'basis_function': params.get('basis_function') if params else None,
                'created_at': params.get('created_at') if params else None
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

    def _parse_cube_filename(self, filename: str) -> Optional[tuple[int, int]]:
        """
        Parse CUBE filename to extract orbital index and grid size.
        
        Args:
            filename: CUBE filename (e.g., "orbital_5_grid80.cube")
        
        Returns:
            Tuple of (orbital_index, grid_size) if parsing succeeds, None otherwise
        """
        try:
            # Remove .cube extension and split
            parts = filename.replace('.cube', '').split('_')
            if len(parts) >= 3 and parts[0] == 'orbital' and parts[2].startswith('grid'):
                orbital_index = int(parts[1])
                grid_size = int(parts[2].replace('grid', ''))
                return orbital_index, grid_size
        except (ValueError, IndexError):
            pass
        return None

    def _get_cube_file_pattern(self, orbital_dir: Path, orbital_index: Optional[int] = None) -> str:
        """
        Generate CUBE file pattern for glob matching.
        
        Args:
            orbital_dir: Directory containing CUBE files
            orbital_index: Specific orbital index, or None for all orbitals
        
        Returns:
            Glob pattern string
        """
        if orbital_index is not None:
            return str(orbital_dir / f"orbital_{orbital_index}_grid*.cube")
        else:
            return str(orbital_dir / "orbital_*_grid*.cube")

    def get_cube_files_info(self, calc_dir: str) -> List[Dict[str, Any]]:
        """Get information about CUBE files in a calculation directory."""
        cube_files = []
        calc_path = Path(calc_dir)
        orbital_dir = calc_path / "orbital"
        
        if not orbital_dir.exists():
            return cube_files
        
        # Find all cube files matching our naming pattern
        import glob
        pattern = self._get_cube_file_pattern(orbital_dir)
        
        for file_path in glob.glob(pattern):
            filename = os.path.basename(file_path)
            
            parsed_result = self._parse_cube_filename(filename)
            if parsed_result is not None:
                orbital_index, grid_size = parsed_result
                    
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
        pattern = self._get_cube_file_pattern(orbital_dir, orbital_index)
        
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

    def move_calculations_directory(self, new_path: str) -> Dict[str, Any]:
        """
        Move all calculation data to a new directory.

        Args:
            new_path: Full directory path for calculations (should include PySCF_calculations)

        Returns:
            Dictionary with move operation results

        Raises:
            ValueError: If new_path is invalid
            OSError: If move operation fails
        """
        new_path_obj = Path(new_path).resolve()
        old_path_obj = self.base_dir.resolve()

        # Validate new path
        if new_path_obj == old_path_obj:
            logger.info("New path is the same as current path, no move needed")
            return {
                "success": True,
                "moved_count": 0,
                "message": "Path unchanged"
            }

        # Check if new path is a subdirectory of old path or vice versa
        # Exception: Allow migration from parent to /PySCF_calculations subfolder
        is_migration_to_subfolder = (
            new_path_obj.name == "PySCF_calculations" and
            new_path_obj.parent == old_path_obj
        )

        if not is_migration_to_subfolder:
            try:
                new_path_obj.relative_to(old_path_obj)
                raise ValueError("New path cannot be a subdirectory of the current path")
            except ValueError as e:
                if "subdirectory" in str(e):
                    raise
                # Not a subdirectory, which is what we want
                pass

            try:
                old_path_obj.relative_to(new_path_obj)
                raise ValueError("Current path cannot be a subdirectory of the new path")
            except ValueError as e:
                if "subdirectory" in str(e):
                    raise
                # Not a subdirectory, which is what we want
                pass

        logger.info(f"Moving calculations from {old_path_obj} to {new_path_obj}")

        # Create new directory if it doesn't exist
        new_path_obj.mkdir(parents=True, exist_ok=True)

        # Check write permissions
        if not os.access(new_path_obj, os.W_OK):
            raise OSError(f"No write permission for directory: {new_path_obj}")

        # Check if new directory is empty
        if list(new_path_obj.iterdir()):
            logger.warning(f"New directory is not empty: {new_path_obj}")

        # Get all calculation directories
        calculations = []
        if old_path_obj.exists():
            for item in old_path_obj.iterdir():
                if item.is_dir():
                    calculations.append(item)

        if not calculations:
            logger.info("No calculations to move")
            # Update base_dir even if no calculations
            self.base_dir = new_path_obj
            return {
                "success": True,
                "moved_count": 0,
                "message": "No calculations to move"
            }

        # Perform the move operation
        moved_count = 0
        failed_moves = []

        for calc_dir in calculations:
            try:
                dest_dir = new_path_obj / calc_dir.name

                # If destination exists, skip with warning
                if dest_dir.exists():
                    logger.warning(f"Destination already exists, skipping: {dest_dir}")
                    failed_moves.append({
                        "name": calc_dir.name,
                        "reason": "Destination already exists"
                    })
                    continue

                # Move directory
                shutil.move(str(calc_dir), str(dest_dir))
                moved_count += 1
                logger.info(f"Moved calculation: {calc_dir.name} -> {dest_dir}")

            except Exception as e:
                logger.error(f"Failed to move {calc_dir.name}: {e}")
                failed_moves.append({
                    "name": calc_dir.name,
                    "reason": str(e)
                })

        # Update base_dir
        self.base_dir = new_path_obj

        # Try to remove old directory if empty
        try:
            if old_path_obj.exists() and not list(old_path_obj.iterdir()):
                old_path_obj.rmdir()
                logger.info(f"Removed old empty directory: {old_path_obj}")
        except Exception as e:
            logger.warning(f"Could not remove old directory: {e}")

        result = {
            "success": len(failed_moves) == 0,
            "moved_count": moved_count,
            "failed_count": len(failed_moves),
            "new_path": str(new_path_obj),
            "old_path": str(old_path_obj)
        }

        if failed_moves:
            result["failed_moves"] = failed_moves
            result["message"] = f"Moved {moved_count} calculations, {len(failed_moves)} failed"
        else:
            result["message"] = f"Successfully moved {moved_count} calculations"

        logger.info(f"Move operation completed: {result['message']}")
        return result

    def set_base_directory(self, new_path: str) -> None:
        """
        Update the base directory path.

        Args:
            new_path: Full directory path for calculations (should include PySCF_calculations)
        """
        self.base_dir = Path(new_path)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Base directory updated to: {self.base_dir}")