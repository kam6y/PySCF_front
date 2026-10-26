from __future__ import annotations

import logging
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

class FrequencyAnalysisMixin:
    """Mixin for vibrational frequency and thermochemistry analysis."""
    def _safe_float_conversion(self, value: Any, param_name: str) -> Optional[float]:
        """
        Safely convert PySCF thermochemistry values to float.
        
        Args:
            value: Value from PySCF thermo function (could be scalar, tuple, array)
            param_name: Parameter name for logging purposes
            
        Returns:
            Float value or None if conversion fails
        """
        try:
            # Handle different data types returned by PySCF
            if isinstance(value, (list, tuple)):
                # Take the first element if it's a sequence
                if len(value) > 0:
                    return float(value[0])
                else:
                    logger.warning(f"Empty sequence for {param_name}")
                    return None
            elif hasattr(value, '__iter__') and not isinstance(value, str):
                # Handle numpy arrays and other iterables
                try:
                    import numpy as np
                    if isinstance(value, np.ndarray):
                        return float(value.flat[0])  # Get first element safely
                    else:
                        # Generic iterable
                        first_val = next(iter(value))
                        return float(first_val)
                except (StopIteration, IndexError):
                    logger.warning(f"Empty iterable for {param_name}")
                    return None
            else:
                # Direct scalar conversion
                return float(value)
        except (ValueError, TypeError) as e:
            logger.warning(f"Failed to convert {param_name} (type: {type(value)}, value: {value}): {e}")
            return None

    def _calculate_frequencies(self) -> Optional[Dict[str, Any]]:
        """Calculate vibrational frequencies and thermochemical properties."""
        if not hasattr(self, 'mf') or self.mf is None:
            return None
        if not hasattr(self, 'mol') or self.mol is None:
            return None
        
        # Initialize partial result structure
        frequency_results = {
            'frequency_analysis_performed': False,
            'vibrational_frequencies': None,
            'imaginary_frequencies_count': None,
            'zero_point_energy': None,
            'thermal_energy_298K': None,
            'entropy_298K': None,
            'gibbs_free_energy_298K': None,
            'heat_capacity_298K': None,
            'normal_modes': None,  # Vibrational mode displacement vectors
            'molecule_structure': None  # Molecular geometry for mode visualization
        }
        
        try:
            logger.info("Starting vibrational frequency analysis...")
            
            # Import thermochemistry module early to catch import errors
            from pyscf.hessian import thermo
            
            # Step 1: Calculate Hessian matrix
            try:
                logger.info("Computing Hessian matrix...")
                hessian_calc = self.mf.Hessian()
                hessian = hessian_calc.kernel()
                logger.info("Hessian calculation completed successfully")
            except Exception as e:
                logger.error(f"Hessian calculation failed: {str(e)}")
                return frequency_results
            
            # Step 2: Perform harmonic analysis
            try:
                logger.info("Performing harmonic analysis...")
                # Use self.mf.mol to ensure harmonic analysis uses the same geometry as Hessian calculation
                freq_info = thermo.harmonic_analysis(self.mf.mol, hessian)
                
                # Extract frequencies in cm^-1
                frequencies_au = freq_info['freq_au']
                frequencies_cm = freq_info['freq_wavenumber']
                normal_modes = freq_info['norm_mode']  # Shape: (3*natom, nmodes)

                # Get molecular structure for mode visualization
                # Use self.mf.mol to ensure we get optimized geometry if optimization was performed
                mol_for_structure = self.mf.mol
                atom_coords = mol_for_structure.atom_coords()  # Shape: (natom, 3)
                atom_symbols = [mol_for_structure.atom_symbol(i) for i in range(mol_for_structure.natm)]

                # Filter out low frequencies (below 80 cm^-1 threshold)
                # and count imaginary frequencies
                imaginary_count = 0
                positive_frequencies = []
                positive_mode_indices = []  # Track which modes correspond to positive frequencies

                for mode_idx, freq_cm in enumerate(frequencies_cm):
                    if freq_cm < 0:
                        # Imaginary frequency (negative eigenvalue)
                        imaginary_count += 1
                    elif freq_cm >= 80.0:  # PySCF default threshold
                        # Real, significant frequency
                        positive_frequencies.append(float(freq_cm))
                        positive_mode_indices.append(mode_idx)

                logger.info(f"Found {len(positive_frequencies)} positive frequencies (≥80 cm⁻¹)")
                logger.info(f"Found {imaginary_count} imaginary frequencies")

                # Extract normal modes for positive frequencies
                # normal_modes shape: (nmodes, 3*natom) - each row is a mode
                # We need to reshape to (natom, 3) for each mode
                natom = mol_for_structure.natm
                formatted_normal_modes = []

                # Debug: Log array shape for troubleshooting
                logger.debug(f"normal_modes shape: {normal_modes.shape}")
                logger.debug(f"natom: {natom}, expected mode vector size: {3*natom}")

                for mode_idx in positive_mode_indices:
                    mode_vector = normal_modes[mode_idx]  # Shape: (3*natom,)
                    # Reshape to (natom, 3) where each row is [dx, dy, dz] for an atom
                    mode_vector_reshaped = mode_vector.reshape(natom, 3)
                    formatted_normal_modes.append(mode_vector_reshaped.tolist())

                # Format molecular structure
                molecule_structure = {
                    'symbols': atom_symbols,
                    'coordinates': atom_coords.tolist()  # Shape: (natom, 3)
                }

                # Update results with frequency information
                frequency_results.update({
                    'frequency_analysis_performed': True,
                    'vibrational_frequencies': positive_frequencies,
                    'imaginary_frequencies_count': imaginary_count,
                    'normal_modes': formatted_normal_modes,  # List of (natom, 3) arrays
                    'molecule_structure': molecule_structure
                })

                # Log optimization quality assessment
                if imaginary_count == 0:
                    logger.info("Geometry optimization successful: no imaginary frequencies detected")
                elif imaginary_count == 1:
                    logger.warning("One imaginary frequency detected: possible transition state")
                else:
                    logger.warning(f"Multiple imaginary frequencies ({imaginary_count}) detected: poor optimization")

            except Exception as e:
                logger.error(f"Harmonic analysis failed: {str(e)}")
                return frequency_results

            # Step 2.5: Calculate IR intensities
            try:
                logger.info("Calculating IR intensities...")

                # Import the appropriate Infrared class based on the method type
                from pyscf.prop import infrared

                # Determine which Infrared class to use based on mean field type
                # Check if UHF/UKS (unrestricted) or RHF/RKS (restricted)
                if hasattr(self.mf, 'mo_occ'):
                    mo_occ = self._as_numpy_array(self.mf.mo_occ)
                    is_unrestricted = isinstance(self.mf.mo_occ, (tuple, list)) or (
                        hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2
                    )
                else:
                    # Fallback: assume restricted
                    is_unrestricted = False

                # Check if DFT (has xc attribute) or HF
                is_dft = hasattr(self.mf, 'xc') and self.mf.xc is not None

                # Select the appropriate Infrared class
                if is_unrestricted:
                    if is_dft:
                        logger.info("Using UKS Infrared calculator")
                        ir_calculator = infrared.uks.Infrared(self.mf)
                    else:
                        logger.info("Using UHF Infrared calculator")
                        ir_calculator = infrared.uhf.Infrared(self.mf)
                else:
                    if is_dft:
                        logger.info("Using RKS Infrared calculator")
                        ir_calculator = infrared.rks.Infrared(self.mf)
                    else:
                        logger.info("Using RHF Infrared calculator")
                        ir_calculator = infrared.rhf.Infrared(self.mf)

                # Run IR intensity calculation
                ir_result = ir_calculator.run()

                # Extract IR intensities (units: km/mol)
                # ir_inten is a 1D array with length equal to number of vibrational modes
                ir_intensities_all = ir_result.ir_inten

                logger.info(f"IR intensities calculated: {len(ir_intensities_all)} modes total")

                # Extract only the IR intensities corresponding to positive frequencies
                # positive_mode_indices contains the indices of modes with freq >= 80 cm^-1
                ir_intensities_filtered = [float(ir_intensities_all[idx]) for idx in positive_mode_indices]

                # Store IR intensities in results
                frequency_results['ir_intensities'] = ir_intensities_filtered

                logger.info(f"IR intensities extracted for {len(ir_intensities_filtered)} positive frequency modes")

            except ImportError as e:
                logger.warning(f"IR intensity calculation skipped: pyscf.prop.infrared module not available ({str(e)})")
                logger.warning("To enable IR intensities, install: pip install git+https://github.com/pyscf/properties")
                frequency_results['ir_intensities'] = None
            except Exception as e:
                logger.warning(f"IR intensity calculation failed: {str(e)}")
                logger.warning("Continuing without IR intensities - uniform intensities will be used for IR spectrum")
                frequency_results['ir_intensities'] = None
            
            # Step 3: Calculate thermochemical properties (optional, non-critical)
            try:
                temperature = 298.15  # K
                pressure = 101325     # Pa
                
                logger.info("Computing thermochemical properties at 298.15 K...")
                thermo_info = thermo.thermo(self.mf, frequencies_au, temperature, pressure)
                
                # Safely extract thermochemical properties
                logger.debug(f"Thermo info keys: {list(thermo_info.keys())}")
                
                zero_point_energy = self._safe_float_conversion(thermo_info.get('ZPE'), 'ZPE')
                thermal_energy = self._safe_float_conversion(thermo_info.get('E_tot'), 'E_tot')
                entropy = self._safe_float_conversion(thermo_info.get('S_tot'), 'S_tot')
                gibbs_free_energy = self._safe_float_conversion(thermo_info.get('G_tot'), 'G_tot')
                heat_capacity = self._safe_float_conversion(thermo_info.get('Cv_tot'), 'Cv_tot')
                
                # Log successful conversions
                successful_conversions = []
                if zero_point_energy is not None:
                    successful_conversions.append('ZPE')
                if thermal_energy is not None:
                    successful_conversions.append('E_tot')
                if entropy is not None:
                    successful_conversions.append('S_tot')
                if gibbs_free_energy is not None:
                    successful_conversions.append('G_tot')
                if heat_capacity is not None:
                    successful_conversions.append('Cv_tot')
                
                logger.info(f"Successfully converted thermochemical properties: {successful_conversions}")
                
                # Update results with thermochemical data
                frequency_results.update({
                    'zero_point_energy': zero_point_energy,
                    'thermal_energy_298K': thermal_energy,
                    'entropy_298K': entropy,
                    'gibbs_free_energy_298K': gibbs_free_energy,
                    'heat_capacity_298K': heat_capacity
                })
                
                logger.info("Thermochemical analysis completed successfully")
                
            except Exception as e:
                logger.warning(f"Thermochemical analysis failed (frequency data preserved): {str(e)}")
                # Continue without thermochemical data - frequency data is still valuable
            
            logger.info("Vibrational frequency analysis completed successfully")
            return frequency_results
            
        except ImportError as e:
            logger.warning(f"Frequency analysis skipped: missing PySCF hessian module ({str(e)})")
            return None
        except Exception as e:
            logger.error(f"Vibrational frequency analysis failed completely: {str(e)}")
            # Return partial results if any were obtained
            if frequency_results['frequency_analysis_performed']:
                logger.info("Returning partial frequency analysis results")
                return frequency_results
            return None

    def _perform_frequency_analysis(self) -> None:
        """Perform vibrational frequency analysis."""
        logger.info("Performing vibrational frequency analysis...")
        frequency_results = self._calculate_frequencies()
        
        if frequency_results is not None:
            self.results.update(frequency_results)
            logger.info("Vibrational frequency analysis completed successfully")
        else:
            logger.warning("Vibrational frequency analysis failed or was skipped")
