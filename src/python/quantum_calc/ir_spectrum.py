"""IR Spectrum calculation and plotting functionality for quantum chemistry calculations.

This module provides functionality to generate theoretical IR spectra from vibrational
frequency calculations, including:
- Scale factor corrections for different computational methods
- Lorentzian broadening for realistic peak shapes
- Plot generation with matplotlib
- Base64 encoding for web display

Based on harmonic approximation with corrections for experimental comparison.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import List, Dict, Any, Optional, Tuple
import logging
import base64
import io
from data.scale_factors import get_recommended_scale_factor

logger = logging.getLogger(__name__)

# Set matplotlib backend to non-interactive for server environments
plt.switch_backend('Agg')


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
    
    def generate_plot(self, 
                     spectrum_data: Dict[str, Any],
                     title: Optional[str] = None,
                     figure_size: Tuple[float, float] = (12, 6),
                     dpi: int = 100,
                     show_peaks: bool = True,
                     peak_threshold: float = 0.05) -> str:
        """
        Generate IR spectrum plot and return as base64 encoded image.
        
        Args:
            spectrum_data: Output from calculate_ir_spectrum()
            title: Plot title (auto-generated if None)
            figure_size: Figure size in inches (width, height)
            dpi: Figure DPI for resolution
            show_peaks: Whether to mark individual peaks
            peak_threshold: Minimum relative intensity to show peak markers
            
        Returns:
            Base64 encoded PNG image
        """
        try:
            plt.style.use('default')
            fig, ax = plt.subplots(figsize=figure_size, dpi=dpi)
            
            # Extract data
            x_axis = np.array(spectrum_data['x_axis'])
            spectrum = np.array(spectrum_data['spectrum'])
            peaks = spectrum_data['peaks']
            metadata = spectrum_data['metadata']
            
            # Plot spectrum
            ax.plot(x_axis, spectrum, 'b-', linewidth=1.5, label='IR Spectrum')
            
            # Mark significant peaks if requested
            if show_peaks and peaks:
                max_intensity = np.max(spectrum) if len(spectrum) > 0 else 1.0
                threshold = max_intensity * peak_threshold
                
                marked_peaks = 0
                for peak in peaks:
                    freq = peak['frequency_cm']
                    intensity = peak['intensity']
                    
                    # Find the spectrum intensity at this frequency
                    freq_idx = np.argmin(np.abs(x_axis - freq))
                    spectrum_intensity = spectrum[freq_idx]
                    
                    if spectrum_intensity >= threshold:
                        ax.axvline(x=freq, color='red', linestyle='--', alpha=0.6, linewidth=0.8)
                        ax.plot(freq, spectrum_intensity, 'ro', markersize=4)
                        
                        # Add frequency label for the strongest peaks
                        if spectrum_intensity >= max_intensity * 0.3:  # Only label strong peaks
                            ax.annotate(f'{freq:.0f}', 
                                      xy=(freq, spectrum_intensity),
                                      xytext=(5, 5), 
                                      textcoords='offset points',
                                      fontsize=8, 
                                      color='red')
                        marked_peaks += 1
                
                logger.info(f"Marked {marked_peaks} significant peaks")
            
            # Set labels and title
            ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
            ax.set_ylabel('Intensity (arb. units)', fontsize=12)
            
            if title is None:
                title = f"IR Spectrum - {metadata['method']}/{metadata['basis_set']}"
                title += f"\nScale factor: {metadata['scale_factor']:.3f}, "
                title += f"Broadening FWHM: {metadata['broadening_fwhm_cm']:.0f} cm⁻¹"
            
            ax.set_title(title, fontsize=14, pad=20)
            
            # Invert x-axis (standard for IR spectra)
            ax.invert_xaxis()
            
            # Set y-axis to start from 0
            ax.set_ylim(bottom=0)
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
            # Add metadata text box
            info_text = f"Method: {metadata['method']}\n"
            info_text += f"Basis set: {metadata['basis_set']}\n"
            info_text += f"Scale factor: {metadata['scale_factor']:.3f}\n"
            info_text += f"Peaks shown: {metadata['num_peaks_in_range']}/{metadata['num_peaks_total']}"
            
            ax.text(0.98, 0.98, info_text,
                   transform=ax.transAxes,
                   verticalalignment='top',
                   horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                   fontsize=9)
            
            # Adjust layout
            plt.tight_layout()
            
            # Convert to base64
            buffer = io.BytesIO()
            plt.savefig(buffer, format='png', dpi=dpi, bbox_inches='tight')
            buffer.seek(0)
            
            image_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
            
            plt.close(fig)  # Clean up
            
            logger.info("IR spectrum plot generated successfully")
            return image_base64
            
        except Exception as e:
            logger.error(f"Plot generation failed: {str(e)}")
            plt.close('all')  # Clean up any open figures
            raise
    
    def generate_spectrum_with_plot(self,
                                  frequencies: List[float],
                                  intensities: Optional[List[float]] = None,
                                  broadening_fwhm: float = 100.0,
                                  x_range: Tuple[float, float] = (400, 4000),
                                  **plot_kwargs) -> Dict[str, Any]:
        """
        Generate IR spectrum data and plot in one call.
        
        Args:
            frequencies: Vibrational frequencies in cm⁻¹ (unscaled)
            intensities: IR intensities (if None, uniform intensities are used)
            broadening_fwhm: Full width at half maximum for Lorentzian broadening (cm⁻¹)
            x_range: Frequency range for spectrum (min_freq, max_freq) in cm⁻¹
            **plot_kwargs: Additional arguments for plot generation
            
        Returns:
            Dictionary containing spectrum data, plot image, and metadata
        """
        try:
            # Calculate spectrum
            spectrum_data = self.calculate_ir_spectrum(
                frequencies, intensities, broadening_fwhm, x_range
            )
            
            # Generate plot
            plot_image = self.generate_plot(spectrum_data, **plot_kwargs)
            
            # Combine results
            result = {
                'spectrum_data': spectrum_data,
                'plot_image_base64': plot_image,
                'success': True
            }
            
            logger.info("IR spectrum with plot generated successfully")
            return result
            
        except Exception as e:
            logger.error(f"IR spectrum generation failed: {str(e)}")
            return {
                'spectrum_data': None,
                'plot_image_base64': None,
                'success': False,
                'error': str(e)
            }


def create_ir_spectrum_from_calculation_results(calculation_results: Dict[str, Any],
                                              **kwargs) -> Dict[str, Any]:
    """
    Create IR spectrum from PySCF calculation results.
    
    Args:
        calculation_results: Results dictionary from quantum calculation
        **kwargs: Additional arguments for spectrum generation
        
    Returns:
        Dictionary containing spectrum data and plot
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
        
        # Generate spectrum
        result = calculator.generate_spectrum_with_plot(frequencies, **kwargs)
        
        return result
        
    except Exception as e:
        logger.error(f"Failed to create IR spectrum from calculation results: {str(e)}")
        return {
            'spectrum_data': None,
            'plot_image_base64': None,
            'success': False,
            'error': str(e)
        }