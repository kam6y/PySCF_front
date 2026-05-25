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
            for item in sorted(old_path_obj.iterdir(), key=lambda path: path.name):
                if is_migration_to_subfolder and item.resolve() == new_path_obj:
                    continue
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

        conflicting_moves = []
        for calc_dir in calculations:
            dest_dir = new_path_obj / calc_dir.name
            if dest_dir.exists():
                logger.warning(f"Destination already exists, cannot move: {dest_dir}")
                conflicting_moves.append({
                    "name": calc_dir.name,
                    "reason": "Destination already exists"
                })

        if conflicting_moves:
            return {
                "success": False,
                "moved_count": 0,
                "failed_count": len(conflicting_moves),
                "new_path": str(new_path_obj),
                "old_path": str(old_path_obj),
                "failed_moves": conflicting_moves,
                "message": (
                    f"Cannot move calculations: {len(conflicting_moves)} "
                    "destination conflicts detected"
                )
            }

        # Perform the move operation
        moved_directories = []
        failed_moves = []

        for calc_dir in calculations:
            try:
                dest_dir = new_path_obj / calc_dir.name

                # Move directory
                shutil.move(str(calc_dir), str(dest_dir))
                moved_directories.append((dest_dir, calc_dir))
                logger.info(f"Moved calculation: {calc_dir.name} -> {dest_dir}")

            except Exception as e:
                logger.error(f"Failed to move {calc_dir.name}: {e}")
                failed_moves.append({
                    "name": calc_dir.name,
                    "reason": str(e)
                })
                break

        rollback_errors = []
        if failed_moves:
            logger.warning("Move failed; rolling back moved calculation directories")
            for dest_dir, original_dir in reversed(moved_directories):
                try:
                    shutil.move(str(dest_dir), str(original_dir))
                    logger.info(f"Rolled back calculation: {dest_dir.name} -> {original_dir}")
                except Exception as e:
                    logger.error(f"Failed to roll back {dest_dir.name}: {e}")
                    rollback_errors.append({
                        "name": dest_dir.name,
                        "reason": str(e)
                    })

            remaining_moved_count = sum(
                1 for dest_dir, _ in moved_directories if dest_dir.exists()
            )
            result = {
                "success": False,
                "moved_count": remaining_moved_count,
                "failed_count": len(failed_moves),
                "new_path": str(new_path_obj),
                "old_path": str(old_path_obj),
                "failed_moves": failed_moves,
                "message": (
                    f"Move failed after {len(moved_directories)} calculations; "
                    f"rolled back {len(moved_directories) - len(rollback_errors)}"
                )
            }

            if rollback_errors:
                result["rollback_errors"] = rollback_errors

            logger.info(f"Move operation completed: {result['message']}")
            return result

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
            "moved_count": len(moved_directories),
            "failed_count": len(failed_moves),
            "new_path": str(new_path_obj),
            "old_path": str(old_path_obj)
        }

        if failed_moves:
            result["failed_moves"] = failed_moves
            result["message"] = (
                f"Moved {len(moved_directories)} calculations, "
                f"{len(failed_moves)} failed"
            )
        else:
            result["message"] = f"Successfully moved {len(moved_directories)} calculations"

        logger.info(f"Move operation completed: {result['message']}")
        return result
