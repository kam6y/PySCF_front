import React, { useState, useEffect, useCallback } from 'react';
import styles from './VibrationModeViewer.module.css';
import type { components } from '../types/generated-api';
import { MoleculeViewer } from './MoleculeViewer';

type IRSpectrumData = components['schemas']['IRSpectrumData'];
type IRPeak = components['schemas']['IRPeak'];
type AtomDisplacement = components['schemas']['AtomDisplacement'];

interface VibrationModeViewerProps {
  spectrumData: IRSpectrumData | null;
  optimizedGeometry?: string;
  selectedPeakIndex: number | null;
  selectedVibrationMode: AtomDisplacement[] | null;
  onPeakSelect: (peak: IRPeak, peakIndex: number) => void;
  onClearSelection: () => void;
  settings: { x_min: number; x_max: number; show_peaks: boolean };
}

export const VibrationModeViewer: React.FC<VibrationModeViewerProps> =
  React.memo(
    ({
      spectrumData,
      optimizedGeometry,
      selectedPeakIndex,
      selectedVibrationMode,
      onPeakSelect,
      onClearSelection,
      settings,
    }) => {
      const [moleculeXYZ, setMoleculeXYZ] = useState<string | null>(null);

      /**
       * 分子構造の座標データを生成
       *
       * 優先度の高い順に以下のデータソースを使用：
       *
       * 【優先度1】optimizedGeometry（ジオメトリ最適化後の構造）
       *   - 最適化計算が実行された場合に利用可能
       *   - より正確な分子構造を反映
       *   - 振動解析の基準構造として使用される
       *
       * 【優先度2】mode_displacements（振動モードの変位データ）
       *   - 各ピークに付属する原子座標情報
       *   - 最適化構造が利用できない場合のフォールバック
       *   - 入力構造をベースにした振動モード情報
       *
       * この優先順位により、常に最も信頼性の高い構造データを表示します。
       */
      useEffect(() => {
        // 優先度1: 最適化済み構造を使用（利用可能な場合）
        if (optimizedGeometry) {
          setMoleculeXYZ(optimizedGeometry);
          return;
        }

        // 優先度2: スペクトルデータの mode_displacements からフォールバック
        if (!spectrumData) {
          setMoleculeXYZ(null);
          return;
        }

        const { peaks } = spectrumData.spectrum;
        // 振動モードデータを持つ最初のピークから分子構造を読み込み
        const firstPeakWithMode = peaks.find(
          p => p.mode_displacements && p.mode_displacements.length > 0
        );
        if (firstPeakWithMode && firstPeakWithMode.mode_displacements) {
          // mode_displacements から XYZ 形式の文字列を生成
          const atoms = firstPeakWithMode.mode_displacements;
          const numAtoms = atoms.length;
          const xyzLines = [`${numAtoms}`, 'Molecule Structure'];
          atoms.forEach(atom => {
            xyzLines.push(
              `${atom.element} ${atom.x.toFixed(6)} ${atom.y.toFixed(6)} ${atom.z.toFixed(6)}`
            );
          });
          const xyzString = xyzLines.join('\n');
          setMoleculeXYZ(xyzString);
        } else {
          setMoleculeXYZ(null);
        }
      }, [spectrumData, optimizedGeometry]);

      const handlePeakClick = useCallback(
        (peak: IRPeak, peakIndex: number) => {
          if (peak.mode_displacements && peak.mode_displacements.length > 0) {
            onPeakSelect(peak, peakIndex);
          }
        },
        [onPeakSelect]
      );

      const formatFrequency = (freq: number) => {
        return freq.toFixed(1);
      };

      const formatIntensity = (intensity: number) => {
        return intensity.toFixed(2);
      };

      if (!spectrumData) {
        return (
          <div className={styles.container}>
            <div className={styles.noDataMessage}>
              No spectrum data available. Please load an IR spectrum first.
            </div>
          </div>
        );
      }

      const { spectrum } = spectrumData;
      const { peaks } = spectrum;

      return (
        <div className={styles.container}>
          {/* Flex wrapper for side-by-side layout */}
          <div className={styles.vibrationFlexWrapper}>
            {/* Left Column: Peak Information */}
            <div className={styles.peakListColumn}>
              <h3 className={styles.peakTitle}>Peak Information</h3>
              <div className={styles.sectionDescription}>
                Select a peak to visualize the corresponding vibrational mode.
              </div>
              {settings.show_peaks && peaks.length > 0 ? (
                <div className={styles.peaksTableWrapper}>
                  <table className={styles.peaksTable}>
                    <thead>
                      <tr>
                        <th>Frequency (cm⁻¹)</th>
                        <th>Intensity</th>
                        <th>Original Freq.</th>
                      </tr>
                    </thead>
                    <tbody>
                      {peaks
                        .filter(
                          peak =>
                            peak.frequency_cm >= settings.x_min &&
                            peak.frequency_cm <= settings.x_max
                        )
                        .sort((a, b) => b.intensity - a.intensity)
                        .map((peak, index) => {
                          const hasVibrationData =
                            peak.mode_displacements &&
                            peak.mode_displacements.length > 0;
                          const isSelected = selectedPeakIndex === index;
                          return (
                            <tr
                              key={index}
                              className={`${styles.peaksRow} ${
                                hasVibrationData ? styles.clickableRow : ''
                              } ${isSelected ? styles.selectedRow : ''}`}
                              onClick={() =>
                                hasVibrationData && handlePeakClick(peak, index)
                              }
                              style={{
                                cursor: hasVibrationData
                                  ? 'pointer'
                                  : 'default',
                              }}
                            >
                              <td>{formatFrequency(peak.frequency_cm)}</td>
                              <td>{formatIntensity(peak.intensity)}</td>
                              <td>
                                {formatFrequency(peak.original_frequency_cm)}
                              </td>
                            </tr>
                          );
                        })}
                    </tbody>
                  </table>
                </div>
              ) : (
                <div className={styles.noDataMessage}>
                  No peaks available. Please adjust the settings.
                </div>
              )}
            </div>

            {/* Right Column: Vibration Mode Visualization */}
            <div className={styles.vibrationViewerColumn}>
              <div className={styles.viewerHeader}>
                <h3>Vibration Mode Visualization</h3>
                {selectedVibrationMode && (
                  <button
                    onClick={onClearSelection}
                    className={styles.clearButton}
                    title="Clear selection"
                  >
                    ✕
                  </button>
                )}
              </div>
              <div className={styles.sectionDescription}>
                Interactive 3D visualization of the selected vibrational mode.
              </div>
              <div className={styles.moleculeViewerContainer}>
                <MoleculeViewer
                  width="100%"
                  height="500px"
                  backgroundColor="#f8f9fa"
                  xyzData={moleculeXYZ}
                  vibrationMode={selectedVibrationMode}
                  animationAmplitude={0.3}
                  className={styles.irSpectrumMoleculeViewer}
                />
                {!selectedVibrationMode && (
                  <div className={styles.placeholder}>
                    <p>Select a peak to view its vibration mode</p>
                  </div>
                )}
              </div>
            </div>
          </div>
        </div>
      );
    }
  );
