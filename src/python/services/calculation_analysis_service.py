"""Analysis-derived calculation operations."""

import logging
import os
from datetime import datetime
from typing import Any

from quantum_calc.ir_spectrum import create_ir_spectrum_from_calculation_results

from .calculation_service_context import CalculationServiceContext
from .exceptions import NotFoundError, ServiceError, ValidationError

logger = logging.getLogger(__name__)


class CalculationAnalysisService:
    """Handles derived analysis endpoints such as IR spectra."""

    def __init__(self, context: CalculationServiceContext) -> None:
        self.context = context

    def generate_ir_spectrum(
        self,
        calculation_id: str,
        broadening_fwhm: float = 100.0,
        x_min: float = 400.0,
        x_max: float = 4000.0,
        show_peaks: bool = True,
    ) -> dict[str, Any]:
        """
        Generate IR spectrum for a calculation.

        Raises:
            NotFoundError: If calculation or results are missing.
            ValidationError: If calculation data cannot produce an IR spectrum.
            ServiceError: If spectrum generation fails.
        """
        try:
            calc_path = self.context.resolve_calculation_path(calculation_id)

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            status = self.context.repository.read_calculation_status(calc_path)
            if status != "completed":
                raise ValidationError(
                    f"Calculation is not completed (status: {status}). "
                    "IR spectrum cannot be generated."
                )

            results = self.context.repository.read_calculation_results(calc_path)
            if not results:
                raise NotFoundError("Calculation results not found.")

            if not results.get("vibrational_frequencies"):
                raise ValidationError(
                    "No vibrational frequency data found. Frequency analysis may "
                    "not have been performed or failed."
                )

            if broadening_fwhm <= 0:
                raise ValidationError("Broadening FWHM must be positive.")

            if x_min >= x_max:
                raise ValidationError("x_min must be less than x_max.")

            logger.info(f"Generating IR spectrum for calculation {calculation_id}")
            logger.info(
                f"Parameters: FWHM={broadening_fwhm}, range=({x_min}, {x_max}), "
                f"show_peaks={show_peaks}"
            )

            ir_result = create_ir_spectrum_from_calculation_results(
                results,
                broadening_fwhm=broadening_fwhm,
                x_range=(x_min, x_max),
            )

            if not ir_result.get("success"):
                error_msg = ir_result.get(
                    "error",
                    "Unknown error occurred during IR spectrum generation",
                )
                logger.error(
                    f"IR spectrum generation failed for {calculation_id}: {error_msg}"
                )
                raise ServiceError(f"IR spectrum generation failed: {error_msg}")

            spectrum_data = ir_result.get("spectrum_data", {})
            plot_image = ir_result.get("plot_image_base64")

            logger.info(f"IR spectrum generated successfully for {calculation_id}")
            logger.info(
                f"Spectrum contains {len(spectrum_data.get('peaks', []))} peaks"
            )

            return {
                "calculation_id": calculation_id,
                "spectrum": {
                    "x_axis": spectrum_data.get("x_axis", []),
                    "y_axis": spectrum_data.get("spectrum", []),
                    "peaks": spectrum_data.get("peaks", []),
                    "metadata": spectrum_data.get("metadata", {}),
                },
                "plot_image_base64": plot_image,
                "generation_info": {
                    "broadening_fwhm_cm": broadening_fwhm,
                    "frequency_range_cm": [x_min, x_max],
                    "peaks_marked": show_peaks,
                    "generated_at": datetime.now().isoformat(),
                },
            }
        except ValueError as e:
            logger.error(
                f"Invalid parameters for IR spectrum generation ({calculation_id}): {e}"
            )
            raise ValidationError(f"Invalid parameters: {str(e)}")
