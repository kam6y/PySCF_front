"""Calculation directory migration utilities."""

import logging
import os
import shutil
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


class CalculationDirectoryMigration:
    """Moves the calculations base directory and updates the active path."""

    def __init__(self, base_dir: Optional[str] = None):
        if base_dir is None:
            # Default to ~/PySCF_calculations if no base_dir is provided
            home = Path.home()
            self.base_dir = home / "PySCF_calculations"
        else:
            # Use the provided path as-is (should already include PySCF_calculations)
            self.base_dir = Path(base_dir)

        self.base_dir.mkdir(parents=True, exist_ok=True)

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
