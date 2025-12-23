import React, {
  useState,
  useEffect,
  useMemo,
  useCallback,
  useRef,
  useDeferredValue,
} from 'react';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
} from 'recharts';
import styles from './IRSpectrumChart.module.css';
import { getIRSpectrum } from '../apiClient';
import type { components } from '../types/generated-api';
import {
  IR_SPECTRUM_DEFAULTS,
  IR_SPECTRUM_CONSTRAINTS,
  IR_SPECTRUM_API_RANGE,
  type IRSettings,
} from '../utils/irSpectrumConstants';

type IRSpectrumData = components['schemas']['IRSpectrumData'];
type IRPeak = components['schemas']['IRPeak'];
type IRSettingsUpdate = Partial<IRSettings>;

interface IRSpectrumChartProps {
  calculationId: string | null;
  onError?: (error: string) => void;
  onSpectrumDataLoaded?: (data: IRSpectrumData) => void;
  selectedPeakIndex?: number | null;
  settings: IRSettings;
  onSettingsChange: React.Dispatch<React.SetStateAction<IRSettings>>;
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

/**
 * Generate evenly spaced tick values for a range
 */
const generateTicks = (
  min: number,
  max: number,
  count: number = 11
): number[] => {
  const step = (max - min) / (count - 1);
  return Array.from({ length: count }, (_, i) => Math.round(max - step * i));
};

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

export const IRSpectrumChart: React.FC<IRSpectrumChartProps> = React.memo(
  ({
    calculationId,
    onError,
    onSpectrumDataLoaded,
    selectedPeakIndex,
    settings,
    onSettingsChange,
  }) => {
    const [spectrumData, setSpectrumData] = useState<IRSpectrumData | null>(
      null
    );
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const abortControllerRef = useRef<AbortController | null>(null);
    const isMountedRef = useRef(true);
    const onSpectrumDataLoadedRef = useRef(onSpectrumDataLoaded);
    const sliderMin = IR_SPECTRUM_CONSTRAINTS.x_min.min;
    const sliderMax = IR_SPECTRUM_CONSTRAINTS.x_max.max;
    const sliderStep = IR_SPECTRUM_CONSTRAINTS.x_min.step;
    const MIN_RANGE_GAP = Math.max(sliderStep, 10);

    // スライダー操作中のパフォーマンス向上のため、x_min/x_maxの更新を遅延
    const deferredSettings = useDeferredValue(settings);

    // Keep ref up to date
    useEffect(() => {
      onSpectrumDataLoadedRef.current = onSpectrumDataLoaded;
    }, [onSpectrumDataLoaded]);

    const updateSettings = useCallback(
      (updates: IRSettingsUpdate) => {
        onSettingsChange(prev => ({
          ...prev,
          ...updates,
        }));
      },
      [onSettingsChange]
    );

    const handleShowPeaksToggle = useCallback(
      (checked: boolean) => {
        updateSettings({ show_peaks: checked });
      },
      [updateSettings]
    );

    const handleRangeUpdate = useCallback(
      (type: 'min' | 'max', rawValue: number) => {
        if (!Number.isFinite(rawValue)) {
          return;
        }

        if (type === 'min') {
          const upperBound = Math.max(
            sliderMin,
            settings.x_max - MIN_RANGE_GAP
          );
          const nextValue = Math.min(Math.max(rawValue, sliderMin), upperBound);
          updateSettings({ x_min: nextValue });
        } else {
          const lowerBound = Math.min(
            sliderMax,
            settings.x_min + MIN_RANGE_GAP
          );
          const nextValue = Math.max(Math.min(rawValue, sliderMax), lowerBound);
          updateSettings({ x_max: nextValue });
        }
      },
      [
        MIN_RANGE_GAP,
        settings.x_max,
        settings.x_min,
        sliderMax,
        sliderMin,
        updateSettings,
      ]
    );

    const handleRangeSliderChange = useCallback(
      (type: 'min' | 'max') => (event: React.ChangeEvent<HTMLInputElement>) => {
        handleRangeUpdate(type, event.target.valueAsNumber);
      },
      [handleRangeUpdate]
    );

    const handleBroadeningInputChange = useCallback(
      (event: React.ChangeEvent<HTMLInputElement>) => {
        const parsed = parseFloat(event.target.value);
        if (Number.isFinite(parsed)) {
          const clamped = Math.min(
            Math.max(parsed, IR_SPECTRUM_CONSTRAINTS.broadening_fwhm.min),
            IR_SPECTRUM_CONSTRAINTS.broadening_fwhm.max
          );
          updateSettings({ broadening_fwhm: clamped });
        } else {
          updateSettings({
            broadening_fwhm: IR_SPECTRUM_DEFAULTS.broadening_fwhm,
          });
        }
      },
      [updateSettings]
    );

    const fetchIRSpectrum = useCallback(async () => {
      if (!calculationId || !isMountedRef.current) return;

      if (abortControllerRef.current) {
        abortControllerRef.current.abort();
      }

      abortControllerRef.current = new AbortController();
      setIsLoading(true);
      setError(null);

      try {
        // 常に固定範囲（0-4500 cm⁻¹）でサーバーにリクエスト
        // ユーザーの設定範囲（settings.x_min/x_max）はクライアント側フィルタリングで適用
        const result = await getIRSpectrum(calculationId, {
          broadening_fwhm: settings.broadening_fwhm,
          x_min: IR_SPECTRUM_API_RANGE.x_min,
          x_max: IR_SPECTRUM_API_RANGE.x_max,
          show_peaks: true,
        });

        if (isMountedRef.current) {
          setSpectrumData(result);
          if (onSpectrumDataLoadedRef.current) {
            onSpectrumDataLoadedRef.current(result);
          }
        }
      } catch (err) {
        if (err instanceof Error && err.name === 'AbortError') {
          return;
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
    }, [settings.broadening_fwhm, calculationId, onError]);

    useEffect(() => {
      fetchIRSpectrum();
    }, [fetchIRSpectrum]);

    useEffect(() => {
      isMountedRef.current = true;

      return () => {
        isMountedRef.current = false;
        if (abortControllerRef.current) {
          abortControllerRef.current.abort();
        }
      };
    }, []);

    const handleRetry = useCallback(() => {
      fetchIRSpectrum();
    }, [fetchIRSpectrum]);

    const filteredChartData = useMemo(() => {
      if (!spectrumData?.spectrum) return [];

      const { x_axis, y_axis } = spectrumData.spectrum;
      const filteredData: ChartDataPoint[] = [];

      for (let i = 0; i < x_axis.length; i++) {
        const wavenumber = x_axis[i];
        if (
          wavenumber >= deferredSettings.x_min &&
          wavenumber <= deferredSettings.x_max
        ) {
          filteredData.push({
            wavenumber,
            intensity: y_axis[i],
          });
        }
      }

      return filteredData;
    }, [spectrumData, deferredSettings.x_min, deferredSettings.x_max]);

    const visiblePeaksCount = useMemo(() => {
      if (!spectrumData?.spectrum?.peaks) return 0;
      return spectrumData.spectrum.peaks.filter(
        peak =>
          peak.frequency_cm >= deferredSettings.x_min &&
          peak.frequency_cm <= deferredSettings.x_max
      ).length;
    }, [spectrumData, deferredSettings.x_min, deferredSettings.x_max]);

    const sliderRange = Math.max(sliderMax - sliderMin, 1);
    const minPercent = ((settings.x_min - sliderMin) / sliderRange) * 100;
    const maxPercent = ((settings.x_max - sliderMin) / sliderRange) * 100;

    if (isLoading) {
      return (
        <div className={styles.loadingContainer}>
          <div className={styles.loadingText}>⚛️ Generating IR spectrum...</div>
        </div>
      );
    }

    if (error) {
      return (
        <div className={styles.errorContainer}>
          <div className={styles.errorText}>❌ {error}</div>
          <button onClick={handleRetry} className={styles.retryButton}>
            Retry
          </button>
        </div>
      );
    }

    if (!spectrumData) {
      return (
        <div className={styles.noDataContainer}>
          No IR spectrum data available
        </div>
      );
    }

    const { spectrum } = spectrumData;
    const { metadata, peaks } = spectrum;

    return (
      <section>
        <div className={styles.settingsPanel}>
          <div className={styles.settingsGrid}>
            <div className={styles.settingItem}>
              <label>Broadening FWHM (cm⁻¹):</label>
              <input
                type="number"
                value={settings.broadening_fwhm}
                onChange={handleBroadeningInputChange}
                min={IR_SPECTRUM_CONSTRAINTS.broadening_fwhm.min}
                max={IR_SPECTRUM_CONSTRAINTS.broadening_fwhm.max}
                step={IR_SPECTRUM_CONSTRAINTS.broadening_fwhm.step}
                className={styles.settingInput}
              />
            </div>
            <div className={styles.settingItem}>
              <label>Wavenumber range (cm⁻¹):</label>
              <div className={styles.rangeSlider}>
                <div className={styles.rangeTrack} />
                <div
                  className={styles.rangeProgress}
                  style={{
                    right: `${Math.max(0, Math.min(100, minPercent))}%`,
                    left: `${Math.max(0, Math.min(100, 100 - maxPercent))}%`,
                  }}
                />
                <input
                  type="range"
                  min={sliderMin}
                  max={sliderMax}
                  step={sliderStep}
                  value={settings.x_min}
                  onChange={handleRangeSliderChange('min')}
                  className={styles.rangeInput}
                />
                <input
                  type="range"
                  min={sliderMin}
                  max={sliderMax}
                  step={sliderStep}
                  value={settings.x_max}
                  onChange={handleRangeSliderChange('max')}
                  className={styles.rangeInput}
                />
              </div>
              <div className={styles.rangeValues}>
                <span>{Math.round(settings.x_max)} cm⁻¹</span>
                <span>{Math.round(settings.x_min)} cm⁻¹</span>
              </div>
            </div>
            <div className={styles.settingItem}>
              <label>
                <input
                  type="checkbox"
                  checked={settings.show_peaks}
                  onChange={e => handleShowPeaksToggle(e.target.checked)}
                />
                Show peak markers
              </label>
            </div>
          </div>
        </div>
        <div className={styles.chartContainer}>
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
                domain={[deferredSettings.x_min, deferredSettings.x_max]}
                reversed={true}
                tick={{ fontSize: 11 }}
                tickFormatter={value => Math.round(value).toString()}
                ticks={generateTicks(
                  deferredSettings.x_min,
                  deferredSettings.x_max,
                  11
                )}
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
                peaks.map((peak, index) => {
                  if (
                    peak.frequency_cm >= deferredSettings.x_min &&
                    peak.frequency_cm <= deferredSettings.x_max
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
                {visiblePeaksCount}/{metadata.num_peaks_total ?? '--'}
              </span>
            </div>
          </div>
        </div>
      </section>
    );
  }
);
