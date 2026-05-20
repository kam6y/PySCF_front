from __future__ import annotations

import logging
import os
from typing import Any, Dict, Optional

from .exceptions import PauseRequestedException
from .pause_manager import pause_manager

logger = logging.getLogger(__name__)

class CheckpointResumeMixin:
    """Mixin for checkpoint path handling and pause/resume callbacks."""
    def get_checkpoint_path(self) -> str:
        """Get the path to the checkpoint file."""
        return os.path.join(self.working_dir, "calculation.chk")

    def _check_pause_requested(self) -> None:
        """
        Check if pause has been requested for this calculation.

        Raises:
            PauseRequestedException: If pause has been requested
        """
        # Check flag file in working directory
        if pause_manager.check_pause_flag_file(self.working_dir):
            logger.info("Pause requested via flag file, raising PauseRequestedException")
            raise PauseRequestedException("Calculation paused by user request")

    def _scf_callback(self, envs: Dict[str, Any]) -> bool:
        """
        Callback function for SCF iterations.

        This is called after each SCF iteration by PySCF.
        Checks for pause requests and raises PauseRequestedException if needed.

        Args:
            envs: Environment dictionary from PySCF containing iteration info

        Returns:
            False to continue calculation, True to stop (but we use exceptions instead)

        Raises:
            PauseRequestedException: If pause has been requested
        """
        # Check for pause request
        self._check_pause_requested()

        # Log progress every few iterations
        if 'cycle' in envs and envs['cycle'] % 5 == 0:
            logger.debug(f"SCF iteration {envs['cycle']}, checking pause status")

        return False  # Continue calculation

    def _geometry_optimization_callback(self, envs: Dict[str, Any]) -> bool:
        """
        Callback function for geometry optimization steps.

        This is called after each geometry optimization step.
        Saves trajectory and checks for pause requests.

        Args:
            envs: Environment dictionary from PySCF geometric optimizer

        Returns:
            False to continue optimization, True to stop (but we use exceptions instead)

        Raises:
            PauseRequestedException: If pause has been requested
        """
        # Check for pause request
        self._check_pause_requested()

        # Save geometry trajectory step if file_manager is available
        if hasattr(self, 'file_manager'):
            mol = envs.get('mol')
            if mol is None:
                return False

            try:
                step_num = envs.get('cycle', 0)
                # Convert current geometry to XYZ string
                atom_symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
                coords = mol.atom_coords(unit="ANG")

                lines = [str(mol.natm)]
                lines.append(f"Optimization step {step_num}")
                for symbol, coord in zip(atom_symbols, coords):
                    lines.append(f"{symbol:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}")

                geometry_xyz = "\n".join(lines)
                self.file_manager.save_geometry_trajectory_step(self.working_dir, step_num, geometry_xyz)
                logger.debug(f"Saved geometry trajectory step {step_num}")
            except Exception as e:
                logger.warning(f"Failed to save geometry trajectory step: {e}")

        return False  # Continue optimization

    def resume_from_checkpoint(self, pause_state: Optional[Dict[str, Any]] = None) -> None:
        """
        Resume calculation from checkpoint file.

        This method configures the SCF calculation to use the checkpoint file
        as an initial guess, enabling true checkpoint-based resume.

        Args:
            pause_state: Optional pause state information containing checkpoint details
        """
        chk_path = self.get_checkpoint_path()

        if not os.path.exists(chk_path):
            logger.warning(f"Checkpoint file not found: {chk_path}")
            logger.info("Will start calculation from scratch")
            return

        logger.info(f"Resuming calculation from checkpoint: {chk_path}")

        # Configure SCF to use checkpoint file as initial guess
        if hasattr(self, 'mf') and self.mf is not None:
            self.mf.init_guess = 'chkfile'
            self.mf.chkfile = chk_path
            logger.info("Configured SCF to use checkpoint file as initial guess")
