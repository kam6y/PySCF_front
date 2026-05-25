"""Calculation artifact operations such as orbitals and CUBE files."""

import logging
import os
from typing import Any

from quantum_calc import CalculationError, FileManagerError
from quantum_calc.orbital_generator import MolecularOrbitalGenerator

from .calculation_service_context import CalculationServiceContext
from .exceptions import NotFoundError, ServiceError, ValidationError

logger = logging.getLogger(__name__)


class CalculationArtifactService:
    """Handles molecular orbital and persisted CUBE artifact operations."""

    ORBITAL_CUBE_GRID_SIZE_MIN = 40
    ORBITAL_CUBE_GRID_SIZE_MAX = 120
    ORBITAL_CUBE_ISOVALUE_POS_MIN = 0.001
    ORBITAL_CUBE_ISOVALUE_POS_MAX = 0.1
    ORBITAL_CUBE_ISOVALUE_NEG_MIN = -0.1
    ORBITAL_CUBE_ISOVALUE_NEG_MAX = -0.001

    def __init__(self, context: CalculationServiceContext) -> None:
        self.context = context

    def get_molecular_orbitals(self, calculation_id: str) -> dict[str, Any]:
        """
        Get molecular orbital information for a calculation.

        Raises:
            NotFoundError: If calculation or orbital data is unavailable.
            ValidationError: If calculation is not completed.
            ServiceError: For calculation or file access failures.
        """
        try:
            calc_path = self.context.resolve_calculation_path(calculation_id)

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            status = self.context.repository.read_calculation_status(calc_path)
            if status != "completed":
                raise ValidationError(
                    f'Calculation "{calculation_id}" is not completed. '
                    f"Status: {status}"
                )

            orbital_generator = MolecularOrbitalGenerator(calc_path)

            if not orbital_generator.validate_calculation():
                raise NotFoundError(
                    "Orbital data is not available or calculation data is invalid."
                )

            orbital_summary = orbital_generator.get_orbital_summary()

            logger.info(
                f"Retrieved orbital information for calculation {calculation_id}"
            )
            logger.info(
                f"Total orbitals: {orbital_summary['total_orbitals']}, "
                f"HOMO: {orbital_summary['homo_index']}, "
                f"LUMO: {orbital_summary['lumo_index']}"
            )

            return orbital_summary

        except CalculationError as e:
            logger.error(f"Calculation error getting orbitals for {calculation_id}: {e}")
            raise ServiceError(str(e))
        except FileManagerError as e:
            logger.error(f"File manager error getting orbitals for {calculation_id}: {e}")
            raise ServiceError("Unable to access calculation files.")

    def generate_orbital_cube(
        self,
        calculation_id: str,
        orbital_index: int,
        grid_size: int = 80,
        isovalue_pos: float | None = None,
        isovalue_neg: float | None = None,
    ) -> dict[str, Any]:
        """
        Generate CUBE file for a specific molecular orbital.

        Raises:
            NotFoundError: If calculation or orbital data is unavailable.
            ValidationError: If parameters are invalid or calculation is incomplete.
            ServiceError: For generation or file access failures.
        """
        try:
            self.validate_orbital_cube_parameters(
                grid_size=grid_size,
                isovalue_pos=isovalue_pos,
                isovalue_neg=isovalue_neg,
            )

            calc_path = self.context.resolve_calculation_path(calculation_id)

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            status = self.context.repository.read_calculation_status(calc_path)
            if status != "completed":
                raise ValidationError(
                    f'Calculation "{calculation_id}" is not completed. '
                    f"Status: {status}"
                )

            orbital_generator = MolecularOrbitalGenerator(calc_path)

            if not orbital_generator.validate_calculation():
                raise NotFoundError(
                    "Orbital data is not available or calculation data is invalid."
                )

            logger.info(
                f"Generating CUBE file for calculation {calculation_id}, "
                f"orbital {orbital_index}"
            )
            logger.info(
                f"Parameters: grid_size={grid_size}, "
                f"isovalue_pos={isovalue_pos}, isovalue_neg={isovalue_neg}"
            )

            cube_data = orbital_generator.generate_cube_file(
                orbital_index=orbital_index,
                grid_size=grid_size,
                isovalue_pos=isovalue_pos,
                isovalue_neg=isovalue_neg,
                return_content=True,
                save_to_disk=True,
            )

            if cube_data.get("cached", False):
                logger.info(
                    f"Using cached CUBE file for calculation {calculation_id}, "
                    f"orbital {orbital_index}"
                )
            else:
                logger.info(
                    f"Successfully generated CUBE file for calculation {calculation_id}, "
                    f"orbital {orbital_index}"
                )
            logger.info(
                f"File size: {cube_data['generation_params']['file_size_kb']:.1f} KB"
            )

            return cube_data

        except CalculationError as e:
            logger.error(
                f"Calculation error generating CUBE for {calculation_id}, "
                f"orbital {orbital_index}: {e}"
            )
            if "invalid orbital index" in str(e).lower():
                raise ValidationError(str(e))
            raise ServiceError(str(e))
        except FileManagerError as e:
            logger.error(
                f"File manager error generating CUBE for {calculation_id}, "
                f"orbital {orbital_index}: {e}"
            )
            raise ServiceError("Unable to access calculation files.")
        except ValueError as e:
            logger.warning(f"Invalid orbital index for calculation {calculation_id}: {e}")
            raise ValidationError("Invalid orbital index.")

    def list_cube_files(self, calculation_id: str) -> dict[str, Any]:
        """
        List all CUBE files for a calculation.

        Raises:
            NotFoundError: If calculation not found.
            ServiceError: For unexpected file access failures.
        """
        try:
            calc_path = self.context.resolve_calculation_path(calculation_id)

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            cube_files = self.context.cube_service.get_cube_files_info(calc_path)

            logger.info(
                f"Found {len(cube_files)} CUBE files for calculation {calculation_id}"
            )

            return {
                "calculation_id": calculation_id,
                "cube_files": cube_files,
                "total_files": len(cube_files),
                "total_size_kb": sum(f["file_size_kb"] for f in cube_files),
            }
        except (NotFoundError, ValidationError):
            raise
        except Exception as e:
            logger.error(
                f"Error listing CUBE files for {calculation_id}: {e}",
                exc_info=True,
            )
            raise ServiceError("An internal error occurred.")

    def delete_cube_files(
        self,
        calculation_id: str,
        orbital_index: int | None = None,
    ) -> dict[str, Any]:
        """
        Delete CUBE files for a calculation.

        Raises:
            NotFoundError: If calculation not found.
            ServiceError: For unexpected file access failures.
        """
        try:
            calc_path = self.context.resolve_calculation_path(calculation_id)

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            deleted_count = self.context.cube_service.delete_cube_files(
                calc_path,
                orbital_index,
            )

            if deleted_count > 0:
                if orbital_index is not None:
                    logger.info(
                        f"Deleted {deleted_count} CUBE files for orbital "
                        f"{orbital_index} in calculation {calculation_id}"
                    )
                    message = (
                        f"Deleted {deleted_count} CUBE files for orbital "
                        f"{orbital_index}."
                    )
                else:
                    logger.info(
                        f"Deleted {deleted_count} CUBE files for calculation "
                        f"{calculation_id}"
                    )
                    message = f"Deleted {deleted_count} CUBE files."
            else:
                if orbital_index is not None:
                    message = f"No CUBE files found for orbital {orbital_index}."
                else:
                    message = "No CUBE files found."

            return {
                "calculation_id": calculation_id,
                "orbital_index": orbital_index,
                "deleted_files": deleted_count,
                "message": message,
            }
        except (NotFoundError, ValidationError):
            raise
        except Exception as e:
            logger.error(
                f"Error deleting CUBE files for {calculation_id}: {e}",
                exc_info=True,
            )
            raise ServiceError("An internal error occurred.")

    def validate_orbital_cube_parameters(
        self,
        grid_size: int,
        isovalue_pos: float | None,
        isovalue_neg: float | None,
    ) -> None:
        """Validate orbital CUBE generation parameters against the API contract."""
        if not self.ORBITAL_CUBE_GRID_SIZE_MIN <= grid_size <= self.ORBITAL_CUBE_GRID_SIZE_MAX:
            raise ValidationError(
                "grid_size must be between "
                f"{self.ORBITAL_CUBE_GRID_SIZE_MIN} and "
                f"{self.ORBITAL_CUBE_GRID_SIZE_MAX}."
            )

        if (
            isovalue_pos is not None
            and not (
                self.ORBITAL_CUBE_ISOVALUE_POS_MIN
                <= isovalue_pos
                <= self.ORBITAL_CUBE_ISOVALUE_POS_MAX
            )
        ):
            raise ValidationError(
                "isovalue_pos must be between "
                f"{self.ORBITAL_CUBE_ISOVALUE_POS_MIN} and "
                f"{self.ORBITAL_CUBE_ISOVALUE_POS_MAX}."
            )

        if (
            isovalue_neg is not None
            and not (
                self.ORBITAL_CUBE_ISOVALUE_NEG_MIN
                <= isovalue_neg
                <= self.ORBITAL_CUBE_ISOVALUE_NEG_MAX
            )
        ):
            raise ValidationError(
                "isovalue_neg must be between "
                f"{self.ORBITAL_CUBE_ISOVALUE_NEG_MIN} and "
                f"{self.ORBITAL_CUBE_ISOVALUE_NEG_MAX}."
            )
