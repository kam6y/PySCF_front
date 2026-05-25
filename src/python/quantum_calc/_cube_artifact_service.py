"""CUBE artifact management (orbital cube files and cleanup)."""

import glob
import logging
import os
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

logger = logging.getLogger(__name__)


class CubeArtifactService:
    """Manages CUBE files and related artifacts for a calculation directory."""

    def __init__(self, base_dir: Optional[str] = None):
        if base_dir is None:
            # Default to ~/PySCF_calculations if no base_dir is provided
            home = Path.home()
            self.base_dir = home / "PySCF_calculations"
        else:
            # Use the provided path as-is (should already include PySCF_calculations)
            self.base_dir = Path(base_dir)

        self.base_dir.mkdir(parents=True, exist_ok=True)

    def set_base_directory(self, new_path: str) -> None:
        self.base_dir = Path(new_path)
        self.base_dir.mkdir(parents=True, exist_ok=True)

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

    def _delete_files(self, file_paths: List[str]) -> int:
        """Delete files and return the count of successful removals."""
        deleted_count = 0
        for file_path in file_paths:
            try:
                os.unlink(file_path)
                deleted_count += 1
            except Exception:
                continue
        return deleted_count

    def get_cube_files_info(self, calc_dir: str) -> List[Dict[str, Any]]:
        """Get information about CUBE files in a calculation directory."""
        cube_files = []
        calc_path = Path(calc_dir)
        orbital_dir = calc_path / "orbital"

        if not orbital_dir.exists():
            return cube_files

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

    def delete_cube_files(self, calc_dir: str, orbital_index: Optional[int] = None) -> int:
        """
        Delete CUBE files from a calculation directory.

        Args:
            calc_dir: Calculation directory path
            orbital_index: Specific orbital index to delete, or None to delete all

        Returns:
            Number of files deleted
        """
        calc_path = Path(calc_dir)
        orbital_dir = calc_path / "orbital"

        if not orbital_dir.exists():
            return 0

        pattern = self._get_cube_file_pattern(orbital_dir, orbital_index)
        return self._delete_files(glob.glob(pattern))
