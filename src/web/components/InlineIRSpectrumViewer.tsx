// src/web/components/InlineIRSpectrumViewer.tsx

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
  ResponsiveContainer,
  ReferenceLine,
} from 'recharts';
import { getIRSpectrum } from '../apiClient';
import type { components } from '../types/generated-api';
import styles from './InlineIRSpectrumViewer.module.css';

type IRSpectrumData = components['schemas']['IRSpectrumData'];
type IRPeak = components['schemas']['IRPeak'];

interface InlineIRSpectrumViewerProps {
  calculation_id: string;
  broadening_fwhm?: number;
  x_min?: number;
  x_max?: number;
  show_peaks?: boolean;
  height?: number;
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

export const InlineIRSpectrumViewer: React.FC<InlineIRSpectrumViewerProps> =
  React.memo(
    ({
      calculation_id,
      broadening_fwhm = 20.0,
      x_min = 400.0,
      x_max = 4000.0,
      show_peaks = true,
      height = 350,
      onError,
    }) => {
      const [spectrumData, setSpectrumData] = useState<IRSpectrumData | null>(
        null
      );
      const [isLoading, setIsLoading] = useState(false);
      const [error, setError] = useState<string | null>(null);
      const abortControllerRef = useRef<AbortController | null>(null);
      const isMountedRef = useRef(true);

      const fetchIRSpectrum = useCallback(async () => {
        if (!calculation_id || !isMountedRef.current) return;

        // Cancel previous request
        if (abortControllerRef.current) {
          abortControllerRef.current.abort();
        }

        abortControllerRef.current = new AbortController();
        setIsLoading(true);
        setError(null);

        try {
          const result = await getIRSpectrum(calculation_id, {
            broadening_fwhm,
            x_min,
            x_max,
            show_peaks,
          });

          if (isMountedRef.current) {
            setSpectrumData(result);
          }
        } catch (err) {
          if (err instanceof Error && err.name === 'AbortError') {
            return; // Ignore cancellation
          }

          if (isMountedRef.current) {
            const errorMessage =
              err instanceof Error ? err.message : 'Unknown error';
            setError(errorMessage);
            onError?.(errorMessage);
          }
        } finally {
          if (isMountedRef.current) {
            setIsLoading(false);
          }
        }
      }, [calculation_id, broadening_fwhm, x_min, x_max, show_peaks, onError]);

      useEffect(() => {
        fetchIRSpectrum();
      }, [fetchIRSpectrum]);

      // Cleanup on unmount
      useEffect(() => {
        isMountedRef.current = true;

        return () => {
          isMountedRef.current = false;
          if (abortControllerRef.current) {
            abortControllerRef.current.abort();
          }
        };
      }, []);

      // Filter chart data for the specified range
      const filteredChartData = useMemo(() => {
        if (!spectrumData?.spectrum) return [];

        const { x_axis, y_axis } = spectrumData.spectrum;
        const filteredData: ChartDataPoint[] = [];

        for (let i = 0; i < x_axis.length; i++) {
          const wavenumber = x_axis[i];
          if (wavenumber >= x_min && wavenumber <= x_max) {
            filteredData.push({
              wavenumber,
              intensity: y_axis[i],
            });
          }
        }

        return filteredData;
      }, [spectrumData, x_min, x_max]);

      if (isLoading) {
        return (
          <div className={styles.container}>
            <div className={styles.loadingContainer}>
              <div className={styles.spinner}></div>
              <div>Generating IR spectrum...</div>
            </div>
          </div>
        );
      }

      if (error) {
        return (
          <div className={styles.container}>
            <div className={styles.errorContainer}>
              <div>❌ Failed to load IR spectrum</div>
              <div className={styles.errorDetail}>{error}</div>
            </div>
          </div>
        );
      }

      if (!spectrumData) {
        return (
          <div className={styles.container}>
            <div className={styles.errorContainer}>
              <div>❌ No IR spectrum data available</div>
            </div>
          </div>
        );
      }

      const { spectrum } = spectrumData;
      const { metadata, peaks } = spectrum;

      return (
        <div className={styles.container}>
          {/* Header with metadata */}
          <div className={styles.header}>
            <div className={styles.title}>IR Spectrum</div>
            <div className={styles.metadata}>
              <span>
                {metadata.method} / {metadata.basis_set}
              </span>
              <span className={styles.separator}>|</span>
              <span>Scale: {metadata.scale_factor.toFixed(3)}</span>
              <span className={styles.separator}>|</span>
              <span>FWHM: {metadata.broadening_fwhm_cm.toFixed(0)} cm⁻¹</span>
              <span className={styles.separator}>|</span>
              <span>
                Peaks: {metadata.num_peaks_in_range}/{metadata.num_peaks_total}
              </span>
            </div>
          </div>

          {/* Chart */}
          <div className={styles.chartContainer}>
            <ResponsiveContainer width="100%" height={height}>
              <LineChart
                data={filteredChartData}
                margin={{
                  top: 10,
                  right: 20,
                  left: 30,
                  bottom: 30,
                }}
              >
                <CartesianGrid strokeDasharray="3 3" opacity={0.3} />
                <XAxis
                  dataKey="wavenumber"
                  type="number"
                  scale="linear"
                  domain={[x_max, x_min]}
                  tick={{ fontSize: 10 }}
                  tickFormatter={value => Math.round(value).toString()}
                  label={{
                    value: 'Wavenumber (cm⁻¹)',
                    position: 'insideBottom',
                    offset: -20,
                    style: { fontSize: 11 },
                  }}
                />
                <YAxis
                  tick={{ fontSize: 10 }}
                  domain={[0, 'dataMax']}
                  tickFormatter={value => {
                    if (value >= 1e6) return (value / 1e6).toFixed(1) + 'M';
                    if (value >= 1e3) return (value / 1e3).toFixed(1) + 'K';
                    return value.toFixed(0);
                  }}
                  label={{
                    value: 'Intensity',
                    angle: -90,
                    position: 'insideLeft',
                    style: { fontSize: 11 },
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

                {show_peaks &&
                  peaks.map((peak: IRPeak, index: number) => {
                    if (
                      peak.frequency_cm >= x_min &&
                      peak.frequency_cm <= x_max
                    ) {
                      return (
                        <ReferenceLine
                          key={index}
                          x={peak.frequency_cm}
                          stroke="#dc2626"
                          strokeDasharray="5 5"
                          strokeOpacity={0.6}
                          label={{
                            value: peak.frequency_cm.toFixed(0),
                            position: 'top',
                            style: { fill: '#dc2626', fontSize: '9px' },
                          }}
                        />
                      );
                    }
                    return null;
                  })}
              </LineChart>
            </ResponsiveContainer>
          </div>

          {/* Peak information table (top 5 peaks) */}
          {show_peaks && peaks.length > 0 && (
            <div className={styles.peaksSection}>
              <div className={styles.peaksTitle}>Major Peaks (Top 5)</div>
              <div className={styles.peaksTable}>
                <div className={styles.peaksHeader}>
                  <span>Frequency (cm⁻¹)</span>
                  <span>Intensity</span>
                  <span>Original</span>
                </div>
                {peaks
                  .filter(
                    (peak: IRPeak) =>
                      peak.frequency_cm >= x_min && peak.frequency_cm <= x_max
                  )
                  .sort((a: IRPeak, b: IRPeak) => b.intensity - a.intensity)
                  .slice(0, 5)
                  .map((peak: IRPeak, index: number) => (
                    <div key={index} className={styles.peaksRow}>
                      <span>{peak.frequency_cm.toFixed(1)}</span>
                      <span>{peak.intensity.toFixed(2)}</span>
                      <span>{peak.original_frequency_cm.toFixed(1)}</span>
                    </div>
                  ))}
              </div>
            </div>
          )}
        </div>
      );
    }
  );
