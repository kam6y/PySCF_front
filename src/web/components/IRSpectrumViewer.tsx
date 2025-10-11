import React, {
  useState,
  useEffect,
  useMemo,
  useCallback,
  useRef,
} from 'react';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceLine,
} from 'recharts';
import styles from './IRSpectrumViewer.module.css';
import { getIRSpectrum } from '../apiClient';
import type { components } from '../types/generated-api';

type IRSpectrumData = components['schemas']['IRSpectrumData'];
type IRPeak = components['schemas']['IRPeak'];

interface IRSpectrumViewerProps {
  calculationId: string | null;
  onError?: (error: string) => void;
}

interface ChartDataPoint {
  wavenumber: number;
  intensity: number;
}

interface CustomTooltipProps {
  active?: boolean;
  payload?: Array<{
    value: number;
    dataKey: string;
  }>;
  label?: number;
  peaks?: IRPeak[];
}

const CustomTooltip: React.FC<CustomTooltipProps> = ({
  active,
  payload,
  label,
  peaks,
}) => {
  if (active && payload && payload.length && label !== undefined) {
    const nearbyPeak = peaks?.find(
      peak => Math.abs(peak.frequency_cm - label) < 10
    );

    return (
      <div className={styles.customTooltip}>
        <p>{`Wavenumber: ${label.toFixed(1)} cm⁻¹`}</p>
        <p>{`Intensity: ${payload[0].value.toFixed(3)}`}</p>
        {nearbyPeak && (
          <div className={styles.peakInfo}>
            <p>
              <strong>Peak: {nearbyPeak.frequency_cm.toFixed(1)} cm⁻¹</strong>
            </p>
            <p>Original: {nearbyPeak.original_frequency_cm.toFixed(1)} cm⁻¹</p>
          </div>
        )}
      </div>
    );
  }
  return null;
};

export const IRSpectrumViewer: React.FC<IRSpectrumViewerProps> = ({
  calculationId,
  onError,
}) => {
  const [spectrumData, setSpectrumData] = useState<IRSpectrumData | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [settings, setSettings] = useState({
    broadening_fwhm: 20.0,
    x_min: 400.0,
    x_max: 4000.0,
    show_peaks: true,
  });
  const [showSettings, setShowSettings] = useState(false);
  const abortControllerRef = useRef<AbortController | null>(null);
  const isMountedRef = useRef(true);

  const fetchIRSpectrum = useCallback(async () => {
    if (!calculationId || !isMountedRef.current) return;

    // 前のリクエストをキャンセル
    if (abortControllerRef.current) {
      abortControllerRef.current.abort();
    }

    abortControllerRef.current = new AbortController();
    setIsLoading(true);
    setError(null);

    try {
      const result = await getIRSpectrum(calculationId, {
        broadening_fwhm: settings.broadening_fwhm,
        x_min: settings.x_min,
        x_max: settings.x_max,
        show_peaks: settings.show_peaks,
      });

      if (isMountedRef.current) {
        setSpectrumData(result);
      }
    } catch (err) {
      if (err instanceof Error && err.name === 'AbortError') {
        return; // キャンセルされた場合は何もしない
      }

      if (isMountedRef.current) {
        const errorMessage =
          err instanceof Error ? err.message : 'Unknown error';
        setError(errorMessage);
        if (onError) {
          onError(errorMessage);
        }
      }
    } finally {
      if (isMountedRef.current) {
        setIsLoading(false);
      }
    }
  }, [calculationId, settings, onError]);

  useEffect(() => {
    fetchIRSpectrum();
  }, [fetchIRSpectrum]);

  // クリーンアップ用のuseEffect
  useEffect(() => {
    isMountedRef.current = true;

    return () => {
      isMountedRef.current = false;
      if (abortControllerRef.current) {
        abortControllerRef.current.abort();
      }
    };
  }, []);

  const handleSettingsUpdate = useCallback(() => {
    setShowSettings(false);
    // fetchIRSpectrumはuseEffectで自動実行される
  }, []);

  const handleRetry = useCallback(() => {
    fetchIRSpectrum();
  }, [fetchIRSpectrum]);

  // メモリ効率を改善：最初から必要な範囲だけを生成
  const filteredChartData = useMemo(() => {
    if (!spectrumData?.spectrum) return [];

    const { x_axis, y_axis } = spectrumData.spectrum;
    const filteredData: ChartDataPoint[] = [];

    for (let i = 0; i < x_axis.length; i++) {
      const wavenumber = x_axis[i];
      if (wavenumber >= settings.x_min && wavenumber <= settings.x_max) {
        filteredData.push({
          wavenumber,
          intensity: y_axis[i],
        });
      }
    }

    return filteredData;
  }, [spectrumData, settings.x_min, settings.x_max]);

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
          <div className={styles.loadingText}>⚛️ Generating IR spectrum...</div>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className={styles.container}>
        <div className={styles.errorContainer}>
          <div className={styles.errorText}>❌ {error}</div>
          <button onClick={handleRetry} className={styles.retryButton}>
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

  const { spectrum } = spectrumData;
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
                onChange={e =>
                  setSettings(prev => ({
                    ...prev,
                    broadening_fwhm: parseFloat(e.target.value) || 100,
                  }))
                }
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
                  onChange={e =>
                    setSettings(prev => ({
                      ...prev,
                      x_min: parseFloat(e.target.value) || 400,
                    }))
                  }
                  min={0}
                  max={10000}
                  step={100}
                  className={styles.settingInput}
                  placeholder="Min"
                />
                <input
                  type="number"
                  value={settings.x_max}
                  onChange={e =>
                    setSettings(prev => ({
                      ...prev,
                      x_max: parseFloat(e.target.value) || 4000,
                    }))
                  }
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
                  onChange={e =>
                    setSettings(prev => ({
                      ...prev,
                      show_peaks: e.target.checked,
                    }))
                  }
                />
                Show peak markers
              </label>
            </div>
          </div>
          <div className={styles.settingsActions}>
            <button
              onClick={handleSettingsUpdate}
              className={styles.applyButton}
            >
              Apply
            </button>
            <button
              onClick={() => setShowSettings(false)}
              className={styles.cancelButton}
            >
              Cancel
            </button>
          </div>
        </div>
      )}

      <div className={styles.content}>
        <div className={styles.chartSection}>
          <ResponsiveContainer width="100%" height={400}>
            <LineChart
              data={filteredChartData}
              margin={{
                top: 20,
                right: 30,
                left: 40,
                bottom: 40,
              }}
            >
              <CartesianGrid strokeDasharray="3 3" opacity={0.3} />
              <XAxis
                dataKey="wavenumber"
                type="number"
                scale="linear"
                domain={[settings.x_max, settings.x_min]}
                tick={{ fontSize: 11 }}
                tickFormatter={value => Math.round(value).toString()}
                ticks={[
                  settings.x_max,
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.1
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.2
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.3
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.4
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.5
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.6
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.7
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.8
                  ),
                  Math.round(
                    settings.x_max - (settings.x_max - settings.x_min) * 0.9
                  ),
                  settings.x_min,
                ]}
                label={{
                  value: 'Wavenumber (cm⁻¹)',
                  position: 'insideBottom',
                  offset: -25,
                  style: { textAnchor: 'middle' },
                }}
              />
              <YAxis
                tick={{ fontSize: 12 }}
                domain={[0, 'dataMax']}
                tickFormatter={value => {
                  if (value >= 1e6) return (value / 1e6).toFixed(1) + 'M';
                  if (value >= 1e3) return (value / 1e3).toFixed(1) + 'K';
                  return value.toFixed(1);
                }}
                label={{
                  value: 'Intensity (arb. units)',
                  angle: -90,
                  position: 'insideLeft',
                  style: { textAnchor: 'middle' },
                }}
              />
              <Tooltip
                content={<CustomTooltip peaks={peaks} />}
                cursor={{ strokeDasharray: '3 3' }}
              />

              <Line
                type="monotone"
                dataKey="intensity"
                stroke="#2563eb"
                strokeWidth={1.5}
                dot={false}
                name="IR Spectrum"
                isAnimationActive={false}
              />

              {settings.show_peaks &&
                peaks.map((peak: IRPeak, index: number) => {
                  if (
                    peak.frequency_cm >= settings.x_min &&
                    peak.frequency_cm <= settings.x_max
                  ) {
                    return (
                      <ReferenceLine
                        key={index}
                        x={peak.frequency_cm}
                        stroke="#dc2626"
                        strokeDasharray="5 5"
                        strokeOpacity={0.7}
                        label={{
                          value: peak.frequency_cm.toFixed(0),
                          position: 'top',
                          style: { fill: '#dc2626', fontSize: '10px' },
                        }}
                      />
                    );
                  }
                  return null;
                })}
            </LineChart>
          </ResponsiveContainer>
        </div>

        <div className={styles.metadataSection}>
          <h4>Analysis Information</h4>
          <div className={styles.metadataGrid}>
            <div className={styles.metadataItem}>
              <span className={styles.metadataLabel}>Method:</span>
              <span className={styles.metadataValue}>{metadata.method}</span>
            </div>
            <div className={styles.metadataItem}>
              <span className={styles.metadataLabel}>Basis Set:</span>
              <span className={styles.metadataValue}>{metadata.basis_set}</span>
            </div>
            <div className={styles.metadataItem}>
              <span className={styles.metadataLabel}>Scale Factor:</span>
              <span className={styles.metadataValue}>
                {metadata.scale_factor.toFixed(3)}
              </span>
            </div>
            <div className={styles.metadataItem}>
              <span className={styles.metadataLabel}>Broadening FWHM:</span>
              <span className={styles.metadataValue}>
                {metadata.broadening_fwhm_cm.toFixed(0)} cm⁻¹
              </span>
            </div>
            <div className={styles.metadataItem}>
              <span className={styles.metadataLabel}>Peaks Shown:</span>
              <span className={styles.metadataValue}>
                {metadata.num_peaks_in_range}/{metadata.num_peaks_total}
              </span>
            </div>
          </div>
        </div>

        {settings.show_peaks && peaks.length > 0 && (
          <div className={styles.peaksSection}>
            <h4>Peak Information</h4>
            <div className={styles.peaksTable}>
              <div className={styles.peaksHeader}>
                <span>Frequency (cm⁻¹)</span>
                <span>Intensity</span>
                <span>Original Freq.</span>
              </div>
              {peaks
                .filter(
                  (peak: IRPeak) =>
                    peak.frequency_cm >= settings.x_min &&
                    peak.frequency_cm <= settings.x_max
                )
                .sort((a: IRPeak, b: IRPeak) => b.intensity - a.intensity)
                .slice(0, 10)
                .map((peak: IRPeak, index: number) => (
                  <div key={index} className={styles.peaksRow}>
                    <span>{formatFrequency(peak.frequency_cm)}</span>
                    <span>{formatIntensity(peak.intensity)}</span>
                    <span>{formatFrequency(peak.original_frequency_cm)}</span>
                  </div>
                ))}
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
