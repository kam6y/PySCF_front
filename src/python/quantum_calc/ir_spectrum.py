"""IR Spectrum calculation functionality for quantum chemistry calculations.

This module provides functionality to generate theoretical IR spectra from vibrational
frequency calculations, including:
- Scale factor corrections for different computational methods
- Lorentzian broadening for realistic peak shapes

Based on harmonic approximation with corrections for experimental comparison.
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from data.scale_factors import get_recommended_scale_factor

logger = logging.getLogger(__name__)


class IRSpectrumCalculator:
    """
    Calculate and generate IR spectra from vibrational frequency data.
    """
    
    def __init__(self, method: str = 'B3LYP', basis_set: str = '6-31G*'):
        """
        Initialize IR spectrum calculator.
        
        Args:
            method: Computational method used for frequency calculation
            basis_set: Basis set used for frequency calculation
        """
        self.method = method
        self.basis_set = basis_set
        self.scale_factor, self.scale_message = get_recommended_scale_factor(method, basis_set)
        
        logger.info(f"IR spectrum calculator initialized for {method}/{basis_set}")
        logger.info(f"Scale factor: {self.scale_factor} ({self.scale_message})")
    
    def calculate_ir_spectrum(self, 
                            frequencies: List[float],
                            intensities: Optional[List[float]] = None,
                            broadening_fwhm: float = 100.0,
                            x_range: Tuple[float, float] = (400, 4000),
                            num_points: int = 3000) -> Dict[str, Any]:
        """
        Calculate IR spectrum from vibrational frequencies and intensities.
        
        Args:
            frequencies: Vibrational frequencies in cm⁻¹ (unscaled)
            intensities: IR intensities (if None, uniform intensities are used)
            broadening_fwhm: Full width at half maximum for Lorentzian broadening (cm⁻¹)
            x_range: Frequency range for spectrum (min_freq, max_freq) in cm⁻¹
            num_points: Number of points in the spectrum
            
        Returns:
            Dictionary containing spectrum data and metadata
        """
        try:
            # Input validation
            if not frequencies:
                raise ValueError("No frequencies provided")
            
            frequencies = np.array(frequencies)
            if len(frequencies) == 0:
                raise ValueError("Empty frequencies array")
            
            # Apply scale factor
            scaled_frequencies = frequencies * self.scale_factor
            
            # Use uniform intensities if not provided
            if intensities is None:
                intensities = np.ones_like(scaled_frequencies)
                logger.info("Using uniform intensities for all peaks")
            else:
                intensities = np.array(intensities)
                if len(intensities) != len(scaled_frequencies):
                    raise ValueError(f"Intensities length ({len(intensities)}) must match frequencies length ({len(scaled_frequencies)})")
            
            # Filter frequencies within the desired range
            freq_mask = (scaled_frequencies >= x_range[0]) & (scaled_frequencies <= x_range[1])
            active_frequencies = scaled_frequencies[freq_mask]
            active_intensities = intensities[freq_mask]
            
            logger.info(f"Processing {len(active_frequencies)} frequencies in range {x_range}")
            
            # Generate x-axis (wavenumber grid)
            x_axis = np.linspace(x_range[0], x_range[1], num_points)
            
            # Calculate spectrum using Lorentzian broadening
            spectrum = self._calculate_lorentzian_spectrum(
                x_axis, active_frequencies, active_intensities, broadening_fwhm
            )
            
            # Prepare peak information
            peaks = []
            for freq, intensity in zip(active_frequencies, active_intensities):
                peaks.append({
                    'frequency_cm': float(freq),
                    'intensity': float(intensity),
                    'original_frequency_cm': float(freq / self.scale_factor)
                })
            
            # Sort peaks by frequency
            peaks.sort(key=lambda p: p['frequency_cm'])
            
            result = {
                'x_axis': x_axis.tolist(),
                'spectrum': spectrum.tolist(),
                'peaks': peaks,
                'metadata': {
                    'method': self.method,
                    'basis_set': self.basis_set,
                    'scale_factor': self.scale_factor,
                    'scale_message': self.scale_message,
                    'broadening_fwhm_cm': broadening_fwhm,
                    'frequency_range_cm': x_range,
                    'num_peaks_total': len(frequencies),
                    'num_peaks_in_range': len(active_frequencies),
                    'num_points': num_points
                }
            }
            
            logger.info(f"IR spectrum calculated successfully with {len(peaks)} peaks")
            return result
            
        except Exception as e:
            logger.error(f"IR spectrum calculation failed: {str(e)}")
            raise
    
    def _calculate_lorentzian_spectrum(self, 
                                     x_axis: np.ndarray,
                                     frequencies: np.ndarray, 
                                     intensities: np.ndarray,
                                     fwhm: float) -> np.ndarray:
        """
        Calculate spectrum using Lorentzian line shape.
        
        Args:
            x_axis: Wavenumber grid
            frequencies: Peak frequencies
            intensities: Peak intensities
            fwhm: Full width at half maximum
            
        Returns:
            Spectrum intensity array
        """
        spectrum = np.zeros_like(x_axis)
        gamma = fwhm / 2  # Half width at half maximum
        
        for freq, intensity in zip(frequencies, intensities):
            # Lorentzian line shape: I(x) = I₀ * (γ²) / ((x - x₀)² + γ²)
            lorentzian = intensity * (gamma**2) / ((x_axis - freq)**2 + gamma**2)
            spectrum += lorentzian
        
        return spectrum
    



def create_ir_spectrum_from_calculation_results(calculation_results: Dict[str, Any],
                                              **kwargs) -> Dict[str, Any]:
    """
    Create IR spectrum from PySCF calculation results.
    
    Args:
        calculation_results: Results dictionary from quantum calculation
        **kwargs: Additional arguments for spectrum generation
        
    Returns:
        Dictionary containing spectrum data only
    """
    try:
        # Extract frequency data
        frequencies = calculation_results.get('vibrational_frequencies')
        if not frequencies:
            raise ValueError("No vibrational frequencies found in calculation results")
        
        # Extract method and basis set
        parameters = calculation_results.get('parameters', {})
        method = parameters.get('calculation_method', 'B3LYP')
        basis_set = parameters.get('basis_function', '6-31G*')
        
        # Handle different method name formats
        if method.upper() == 'DFT':
            # For DFT calculations, use the exchange-correlation functional
            xc_functional = parameters.get('exchange_correlation', 'B3LYP')
            method = xc_functional
        
        logger.info(f"Creating IR spectrum from calculation results: {method}/{basis_set}")
        logger.info(f"Found {len(frequencies)} vibrational frequencies")
        
        # Initialize calculator
        calculator = IRSpectrumCalculator(method=method, basis_set=basis_set)
        
        # Generate spectrum data only
        spectrum_data = calculator.calculate_ir_spectrum(frequencies, **kwargs)
        
        return {
            'spectrum_data': spectrum_data,
            'success': True
        }
        
    except Exception as e:
        logger.error(f"Failed to create IR spectrum from calculation results: {str(e)}")
        return {
            'spectrum_data': None,
            'success': False,
            'error': str(e)
        }