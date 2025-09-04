import React, { useState, useEffect } from 'react';
import styles from './IRSpectrumViewer.module.css';
import { getIRSpectrum } from '../apiClient';
import { IRSpectrumData } from '../types/api-types';


interface IRSpectrumViewerProps {
  calculationId: string;
  onError?: (error: string) => void;
}

export const IRSpectrumViewer: React.FC<IRSpectrumViewerProps> = ({
  calculationId,
  onError
}) => {
  const [spectrumData, setSpectrumData] = useState<IRSpectrumData | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [settings, setSettings] = useState({
    broadening_fwhm: 100.0,
    x_min: 400.0,
    x_max: 4000.0,
    show_peaks: true
  });
  const [showSettings, setShowSettings] = useState(false);

  const fetchIRSpectrum = async () => {
    if (!calculationId) return;

    setIsLoading(true);
    setError(null);

    try {
      const result = await getIRSpectrum(calculationId, {
        broadening_fwhm: settings.broadening_fwhm,
        x_min: settings.x_min,
        x_max: settings.x_max,
        show_peaks: settings.show_peaks
      });
      
      setSpectrumData(result);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Unknown error';
      setError(errorMessage);
      if (onError) {
        onError(errorMessage);
      }
    } finally {
      setIsLoading(false);
    }
  };

  useEffect(() => {
    fetchIRSpectrum();
  }, [calculationId]);

  const handleSettingsUpdate = () => {
    fetchIRSpectrum();
    setShowSettings(false);
  };

  const formatFrequency = (freq: number) => {
    return freq.toFixed(1);
  };

  const formatIntensity = (intensity: number) => {
    return intensity.toFixed(2);
  };

  if (isLoading) {
    return (
      <div className={styles.container}>
        <div className={styles.loadingContainer}>
          <div className={styles.loadingText}>
            ⚛️ Generating IR spectrum...
          </div>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className={styles.container}>
        <div className={styles.errorContainer}>
          <div className={styles.errorText}>❌ {error}</div>
          <button onClick={fetchIRSpectrum} className={styles.retryButton}>
            Retry
          </button>
        </div>
      </div>
    );
  }

  if (!spectrumData) {
    return (
      <div className={styles.container}>
        <div className={styles.noDataContainer}>
          No IR spectrum data available
        </div>
      </div>
    );
  }

  const { spectrum, plot_image_base64 } = spectrumData;
  const { metadata, peaks } = spectrum;

  return (
    <div className={styles.container}>
      <div className={styles.header}>
        <h3>IR Spectrum Analysis</h3>
        <button
          onClick={() => setShowSettings(!showSettings)}
          className={styles.settingsButton}
        >
          Settings
        </button>
      </div>

      {showSettings && (
        <div className={styles.settingsPanel}>
          <div className={styles.settingsGrid}>
            <div className={styles.settingItem}>
              <label>Broadening FWHM (cm⁻¹):</label>
              <input
                type="number"
                value={settings.broadening_fwhm}
                onChange={(e) => setSettings(prev => ({
                  ...prev,
                  broadening_fwhm: parseFloat(e.target.value) || 100
                }))}
                min={0.1}
                max={1000}
                step={10}
                className={styles.settingInput}
              />
            </div>
            <div className={styles.settingItem}>
              <label>Wavenumber range (cm⁻¹):</label>
              <div className={styles.rangeInputs}>
                <input
                  type="number"
                  value={settings.x_min}
                  onChange={(e) => setSettings(prev => ({
                    ...prev,
                    x_min: parseFloat(e.target.value) || 400
                  }))}
                  min={0}
                  max={10000}
                  step={100}
                  className={styles.settingInput}
                  placeholder="Min"
                />
                <input
                  type="number"
                  value={settings.x_max}
                  onChange={(e) => setSettings(prev => ({
                    ...prev,
                    x_max: parseFloat(e.target.value) || 4000
                  }))}
                  min={0}
                  max={10000}
                  step={100}
                  className={styles.settingInput}
                  placeholder="Max"
                />
              </div>
            </div>
            <div className={styles.settingItem}>
              <label>
                <input
                  type="checkbox"
                  checked={settings.show_peaks}
                  onChange={(e) => setSettings(prev => ({
                    ...prev,
                    show_peaks: e.target.checked
                  }))}
                />
                Show peak markers
              </label>
            </div>
          </div>
          <div className={styles.settingsActions}>
            <button onClick={handleSettingsUpdate} className={styles.applyButton}>
              Apply
            </button>
            <button onClick={() => setShowSettings(false)} className={styles.cancelButton}>
              Cancel
            </button>
          </div>
        </div>
      )}

      <div className={styles.content}>
        {/* Plot Section */}
        <div className={styles.plotSection}>
          <div className={styles.plotContainer}>
            <img
              src={`data:image/png;base64,${plot_image_base64}`}
              alt="IR Spectrum"
              className={styles.spectrumPlot}
            />
          </div>
        </div>

        {/* Calculation Info */}
        <div className={styles.infoSection}>
          <div className={styles.infoGrid}>
            <div>
              <strong>Method:</strong> {metadata.method}
            </div>
            <div>
              <strong>Basis Set:</strong> {metadata.basis_set}
            </div>
            <div>
              <strong>Scale Factor:</strong> {metadata.scale_factor.toFixed(3)}
            </div>
            <div>
              <strong>Peaks:</strong> {metadata.num_peaks_in_range}/{metadata.num_peaks_total}
            </div>
          </div>
        </div>

        {/* Main Peaks */}
        <div className={styles.peaksSection}>
          <h4>Vibrational Frequencies</h4>
          <div className={styles.peaksTableContainer}>
            <table className={styles.peaksTable}>
              <thead>
                <tr>
                  <th>Frequency (cm⁻¹)</th>
                  <th>Intensity</th>
                </tr>
              </thead>
              <tbody>
                {peaks
                  .sort((a, b) => b.frequency_cm - a.frequency_cm)
                  .slice(0, 10) // Show only top 10 peaks
                  .map((peak, index) => (
                    <tr key={index}>
                      <td className={styles.frequencyCell}>
                        {formatFrequency(peak.frequency_cm)}
                      </td>
                      <td className={styles.intensityCell}>
                        {formatIntensity(peak.intensity)}
                      </td>
                    </tr>
                  ))}
              </tbody>
            </table>
          </div>
          {peaks.length > 10 && (
            <div className={styles.peaksNote}>
              Showing top 10 peaks out of {peaks.length} total
            </div>
          )}
        </div>
      </div>
    </div>
  );
};