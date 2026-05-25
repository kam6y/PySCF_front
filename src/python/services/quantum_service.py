"""
Quantum chemistry calculation service facade.

This module preserves the public service API used by routes and agent tools while
delegating each responsibility to a focused calculation service.
"""

import logging
from typing import Any

from quantum_calc import CalculationRepository, CubeArtifactService

from .calculation_analysis_service import CalculationAnalysisService
from .calculation_artifact_service import CalculationArtifactService
from .calculation_command_service import CalculationCommandService
from .calculation_query_service import CalculationQueryService
from .calculation_service_context import CalculationServiceContext

logger = logging.getLogger(__name__)


class QuantumService:
    """Facade for quantum chemistry calculation operations."""

    GPU4PYSCF_SUPPORTED_METHODS = CalculationCommandService.GPU4PYSCF_SUPPORTED_METHODS
    TERMINAL_STATUSES = CalculationServiceContext.TERMINAL_STATUSES
    NON_TERMINAL_STATUSES = CalculationServiceContext.NON_TERMINAL_STATUSES
    RESTART_INTERRUPTED_MESSAGE = (
        CalculationServiceContext.RESTART_INTERRUPTED_MESSAGE
    )
    ORBITAL_CUBE_GRID_SIZE_MIN = CalculationArtifactService.ORBITAL_CUBE_GRID_SIZE_MIN
    ORBITAL_CUBE_GRID_SIZE_MAX = CalculationArtifactService.ORBITAL_CUBE_GRID_SIZE_MAX
    ORBITAL_CUBE_ISOVALUE_POS_MIN = (
        CalculationArtifactService.ORBITAL_CUBE_ISOVALUE_POS_MIN
    )
    ORBITAL_CUBE_ISOVALUE_POS_MAX = (
        CalculationArtifactService.ORBITAL_CUBE_ISOVALUE_POS_MAX
    )
    ORBITAL_CUBE_ISOVALUE_NEG_MIN = (
        CalculationArtifactService.ORBITAL_CUBE_ISOVALUE_NEG_MIN
    )
    ORBITAL_CUBE_ISOVALUE_NEG_MAX = (
        CalculationArtifactService.ORBITAL_CUBE_ISOVALUE_NEG_MAX
    )

    def __init__(self) -> None:
        """Initialize QuantumService and its focused service boundaries."""
        self.context = CalculationServiceContext.from_settings()
        self.query_service = CalculationQueryService(self.context)
        self.command_service = CalculationCommandService(
            self.context,
            self.query_service,
        )
        self.artifact_service = CalculationArtifactService(self.context)
        self.analysis_service = CalculationAnalysisService(self.context)

    @property
    def repository(self) -> CalculationRepository:
        """Expose the shared repository for existing service tests and callers."""
        return self.context.repository

    @repository.setter
    def repository(self, repository: CalculationRepository) -> None:
        self.context.repository = repository

    @property
    def cube_service(self) -> CubeArtifactService:
        """Expose the shared cube service for existing service tests and callers."""
        return self.context.cube_service

    @cube_service.setter
    def cube_service(self, cube_service: CubeArtifactService) -> None:
        self.context.cube_service = cube_service

    def update_calculations_directory(self, new_directory: str) -> None:
        """
        Update the calculations directory for all calculation service boundaries.

        Args:
            new_directory: New directory path for calculations.
        """
        logger.info(f"Updating QuantumService calculations directory to: {new_directory}")
        self.context.set_base_directory(new_directory)

        from quantum_calc.file_watcher import update_watcher_base_directory

        update_watcher_base_directory(new_directory)

    def get_supported_parameters(self) -> dict[str, Any]:
        """Get supported quantum chemistry parameters."""
        return self.query_service.get_supported_parameters()

    def validate_calculation_parameters(
        self,
        params: dict[str, Any],
    ) -> str | None:
        """Validate calculation parameters for compatibility and correctness."""
        return self.command_service.validate_calculation_parameters(params)

    def start_calculation(self, params: dict[str, Any]) -> dict[str, Any]:
        """Start a new quantum chemistry calculation."""
        return self.command_service.start_calculation(params)

    def list_calculations(
        self,
        name_query: str | None = None,
        status: str | None = None,
        calculation_method: str | None = None,
        basis_function: str | None = None,
        date_from: str | None = None,
        date_to: str | None = None,
    ) -> dict[str, Any]:
        """List available calculations with optional filtering."""
        return self.query_service.list_calculations(
            name_query=name_query,
            status=status,
            calculation_method=calculation_method,
            basis_function=basis_function,
            date_from=date_from,
            date_to=date_to,
        )

    def get_calculation_details(self, calculation_id: str) -> dict[str, Any]:
        """Get detailed information about a specific calculation."""
        return self.query_service.get_calculation_details(calculation_id)

    def update_calculation(
        self,
        calculation_id: str,
        new_name: str,
    ) -> dict[str, Any]:
        """Update calculation metadata."""
        return self.command_service.update_calculation(calculation_id, new_name)

    def delete_calculation(self, calculation_id: str) -> dict[str, Any]:
        """Delete a calculation and its files."""
        return self.command_service.delete_calculation(calculation_id)

    def get_calculation_status(self) -> dict[str, Any]:
        """Get status information about the calculation system."""
        return self.query_service.get_calculation_status()

    def get_molecular_orbitals(self, calculation_id: str) -> dict[str, Any]:
        """Get molecular orbital information for a calculation."""
        return self.artifact_service.get_molecular_orbitals(calculation_id)

    def generate_orbital_cube(
        self,
        calculation_id: str,
        orbital_index: int,
        grid_size: int = 80,
        isovalue_pos: float | None = None,
        isovalue_neg: float | None = None,
    ) -> dict[str, Any]:
        """Generate CUBE file for a specific molecular orbital."""
        return self.artifact_service.generate_orbital_cube(
            calculation_id,
            orbital_index,
            grid_size=grid_size,
            isovalue_pos=isovalue_pos,
            isovalue_neg=isovalue_neg,
        )

    def list_cube_files(self, calculation_id: str) -> dict[str, Any]:
        """List all CUBE files for a calculation."""
        return self.artifact_service.list_cube_files(calculation_id)

    def delete_cube_files(
        self,
        calculation_id: str,
        orbital_index: int | None = None,
    ) -> dict[str, Any]:
        """Delete CUBE files for a calculation."""
        return self.artifact_service.delete_cube_files(calculation_id, orbital_index)

    def generate_ir_spectrum(
        self,
        calculation_id: str,
        broadening_fwhm: float = 100.0,
        x_min: float = 400.0,
        x_max: float = 4000.0,
        show_peaks: bool = True,
    ) -> dict[str, Any]:
        """Generate IR spectrum for a calculation."""
        return self.analysis_service.generate_ir_spectrum(
            calculation_id,
            broadening_fwhm=broadening_fwhm,
            x_min=x_min,
            x_max=x_max,
            show_peaks=show_peaks,
        )

    def pause_calculation(self, calculation_id: str) -> dict[str, Any]:
        """Pause a running calculation."""
        return self.command_service.pause_calculation(calculation_id)

    def resume_calculation(self, calculation_id: str) -> dict[str, Any]:
        """Resume a paused calculation."""
        return self.command_service.resume_calculation(calculation_id)

    def _build_calculation_instance(
        self,
        calculation_id: str,
        parameters: dict[str, Any],
        status: str,
        waiting_reason: str | None = None,
    ) -> dict[str, Any]:
        """Build calculation instance dict for responses."""
        return self.context.build_calculation_instance(
            calculation_id,
            parameters,
            status,
            waiting_reason,
        )

    def _validate_orbital_cube_parameters(
        self,
        grid_size: int,
        isovalue_pos: float | None,
        isovalue_neg: float | None,
    ) -> None:
        """Validate orbital CUBE generation parameters against the API contract."""
        self.artifact_service.validate_orbital_cube_parameters(
            grid_size=grid_size,
            isovalue_pos=isovalue_pos,
            isovalue_neg=isovalue_neg,
        )

    def _resolve_calculation_path(self, calculation_id: str) -> str:
        """Resolve a calculation ID from an external request into a safe path."""
        return self.context.resolve_calculation_path(calculation_id)

    def _recover_stale_non_terminal_calculations(
        self,
        process_manager: Any | None = None,
    ) -> None:
        """Mark persisted non-terminal calculations as error when no worker owns them."""
        self.context.recover_stale_non_terminal_calculations(process_manager)
