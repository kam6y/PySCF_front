import React, { type Dispatch, type SetStateAction } from 'react';
import styles from '../../pages/CalculationResultsPage.module.css';
import { IRSpectrumChart } from '../IRSpectrumChart';
import { VibrationModeViewer } from '../VibrationModeViewer';
import { LazyViewer } from '../LazyViewer';
import { CalculationResults } from '../../types/api-types';
import { IRSettings } from '../../utils/irSpectrumConstants';
import type { components } from '../../types/generated-api';

type IRSpectrumData = components['schemas']['IRSpectrumData'];
type IRPeak = components['schemas']['IRPeak'];
type AtomDisplacement = components['schemas']['AtomDisplacement'];

interface VibrationalAnalysisSectionProps {
  calculationId: string;
  results: CalculationResults;
  irSpectrumData: IRSpectrumData | null;
  selectedIRPeakIndex: number | null;
  selectedVibrationMode: AtomDisplacement[] | null;
  irSettings: IRSettings;
  onSetIRSettings: Dispatch<SetStateAction<IRSettings>>;
  onSpectrumDataLoaded: (data: IRSpectrumData) => void;
  onPeakSelect: (peak: IRPeak, peakIndex: number) => void;
  onClearSelection: () => void;
  onError: (error: string) => void;
}

export const VibrationalAnalysisSection =
  React.memo<VibrationalAnalysisSectionProps>(
    ({
      calculationId,
      results,
      irSpectrumData,
      selectedIRPeakIndex,
      selectedVibrationMode,
      irSettings,
      onSetIRSettings,
      onSpectrumDataLoaded,
      onPeakSelect,
      onClearSelection,
      onError,
    }) => {
      return (
        <section
          className={`${styles.calculationSection} ${styles.vibrationalSection}`}
        >
          <h2 className={styles.primaryHeader}>Vibrational Analysis</h2>

          {/* IR Spectrum - Split into two sections */}
          {results.vibrational_frequencies &&
            results.vibrational_frequencies.length > 0 && (
              <>
                <div className={styles.irSpectrumSubsection}>
                  <h3>Infrared (IR) Spectrum</h3>
                  <div className={styles.sectionDescription}>
                    Theoretical infrared spectrum generated from vibrational
                    frequency calculations with scale factor corrections and
                    Lorentzian broadening for realistic peak shapes.
                  </div>
                  <LazyViewer>
                    <IRSpectrumChart
                      calculationId={calculationId}
                      onError={onError}
                      onSpectrumDataLoaded={onSpectrumDataLoaded}
                      settings={irSettings}
                      onSettingsChange={onSetIRSettings}
                    />
                  </LazyViewer>
                </div>

                <LazyViewer>
                  <VibrationModeViewer
                    spectrumData={irSpectrumData}
                    optimizedGeometry={results.optimized_geometry}
                    selectedPeakIndex={selectedIRPeakIndex}
                    selectedVibrationMode={selectedVibrationMode}
                    onPeakSelect={onPeakSelect}
                    onClearSelection={onClearSelection}
                    settings={irSettings}
                  />
                </LazyViewer>
              </>
            )}
        </section>
      );
    }
  );
