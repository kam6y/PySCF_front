import React, { useState, useEffect, useCallback, useMemo } from 'react';
import styles from './CalculationResultsPage.module.css';
import { CalculationInstance } from '../types/api-types';
import { useProcessedCalculationResults } from '../hooks/useProcessedCalculationResults';
import { StyleSpec } from '../../types/3dmol';
import type { components } from '../types/generated-api';
import {
  IR_SPECTRUM_DEFAULTS,
  type IRSettings,
} from '../utils/irSpectrumConstants';
import { CalculationSettingsSummary } from '../components/calculation-results/CalculationSettingsSummary';
import { OptimizedStructureSection } from '../components/calculation-results/OptimizedStructureSection';
import { ElectronicPropertiesSection } from '../components/calculation-results/ElectronicPropertiesSection';
import { CASResultsSection } from '../components/calculation-results/CASResultsSection';
import { EnergeticsSection } from '../components/calculation-results/EnergeticsSection';
import { VibrationalAnalysisSection } from '../components/calculation-results/VibrationalAnalysisSection';
import { MolecularOrbitalsSection } from '../components/calculation-results/MolecularOrbitalsSection';
import { TDDFTResultsSection } from '../components/calculation-results/TDDFTResultsSection';
import { TechnicalDetailsSection } from '../components/calculation-results/TechnicalDetailsSection';

type IRSpectrumData = components['schemas']['IRSpectrumData'];
type IRPeak = components['schemas']['IRPeak'];
type AtomDisplacement = components['schemas']['AtomDisplacement'];

interface CalculationResultsPageProps {
  activeCalculation?: CalculationInstance;
  isLoadingDetails?: boolean;
  detailsError?: string | null;
  onCalculationUpdate: (updatedCalculation: CalculationInstance) => void;
}

export const CalculationResultsPage = ({
  activeCalculation,
  isLoadingDetails = false,
  detailsError = null,
  onCalculationUpdate,
}: CalculationResultsPageProps) => {
  const [error, setError] = useState<string | null>(null);
  const [selectedOrbitalIndex, setSelectedOrbitalIndex] = useState<
    number | null
  >(null);

  // Molecule viewer state for optimized structure section
  const [currentStyle, setCurrentStyle] = useState<StyleSpec | null>({
    stick: {},
  });
  const [showAxes, setShowAxes] = useState(false);
  const [showCoordinates, setShowCoordinates] = useState(false);
  const [useAtomicRadii, setUseAtomicRadii] = useState(false);

  // IR Spectrum shared state
  const [irSpectrumData, setIRSpectrumData] = useState<IRSpectrumData | null>(
    null
  );
  const [selectedIRPeakIndex, setSelectedIRPeakIndex] = useState<number | null>(
    null
  );
  const [selectedVibrationMode, setSelectedVibrationMode] = useState<
    AtomDisplacement[] | null
  >(null);
  const [irSettings, setIRSettings] = useState<IRSettings>({
    broadening_fwhm: IR_SPECTRUM_DEFAULTS.broadening_fwhm,
    x_min: IR_SPECTRUM_DEFAULTS.x_min,
    x_max: IR_SPECTRUM_DEFAULTS.x_max,
    show_peaks: IR_SPECTRUM_DEFAULTS.show_peaks,
  });

  useEffect(() => {
    setError(detailsError);
  }, [detailsError]);

  const handleOrbitalSelect = useCallback((orbitalIndex: number) => {
    setSelectedOrbitalIndex(orbitalIndex);
  }, []);

  const handleIRPeakSelect = useCallback((peak: IRPeak, peakIndex: number) => {
    setSelectedIRPeakIndex(peakIndex);
    setSelectedVibrationMode(peak.mode_displacements || null);
  }, []);

  const handleClearVibrationSelection = useCallback(() => {
    setSelectedIRPeakIndex(null);
    setSelectedVibrationMode(null);
  }, []);

  const handleSpectrumDataLoaded = useCallback((data: IRSpectrumData) => {
    setIRSpectrumData(data);
    setIRSettings({
      broadening_fwhm: data.spectrum.metadata.broadening_fwhm_cm,
      x_min: IR_SPECTRUM_DEFAULTS.x_min,
      x_max: IR_SPECTRUM_DEFAULTS.x_max,
      show_peaks: IR_SPECTRUM_DEFAULTS.show_peaks,
    });
  }, []);

  // Process and memoize calculation results data
  const processedData = useProcessedCalculationResults(activeCalculation);

  // Determine if optimized structure section should be shown
  const shouldShowOptimizedStructure = useMemo(() => {
    if (!processedData) return false;

    const { parameters, results } = processedData;

    return (
      ['HF', 'DFT', 'MP2'].includes(parameters.calculation_method) &&
      (parameters as any).optimize_geometry !== false &&
      !!results.optimized_geometry
    );
  }, [processedData]);

  // Keep prop for interface compatibility (currently unused in this container)
  void onCalculationUpdate;

  // Show loading state
  if (isLoadingDetails) {
    return (
      <div className={styles.pageContainer}>
        <div className={styles.pageContent}>
          <div className={styles.loadingContainer}>
            <div className={styles.loadingText}>
              Loading calculation details...
            </div>
          </div>
        </div>
      </div>
    );
  }

  // Show error state
  if (error) {
    return (
      <div className={styles.pageContainer}>
        <div className={styles.pageContent}>
          <div className={styles.errorContainer}>❌ {error}</div>
        </div>
      </div>
    );
  }

  // Show message when no calculation is selected
  if (!activeCalculation) {
    return (
      <div className={styles.pageContainer}>
        <div className={styles.pageContent}>
          <div className={styles.noCalculationContainer}>
            No calculation selected. Please select a calculation from the
            sidebar to view its results.
          </div>
        </div>
      </div>
    );
  }

  // Show message for incomplete calculations
  if (activeCalculation.status !== 'completed' || !activeCalculation.results) {
    const statusMessages = {
      pending: 'This calculation is pending. Please run the calculation first.',
      running:
        'This calculation is currently running. Please wait for completion.',
      pausing: 'This calculation is pausing. Please wait...',
      paused:
        'This calculation is paused. You can resume it from where it was paused.',
      error:
        'This calculation failed. Please check the settings and try again.',
    };

    return (
      <div className={styles.pageContainer}>
        <div className={styles.pageContent}>
          <div className={styles.statusMessageContainer}>
            {statusMessages[
              activeCalculation.status as keyof typeof statusMessages
            ] || 'Calculation results are not available.'}
          </div>
          <div className={styles.statusMessageMeta}>
            <strong>Calculation:</strong> {activeCalculation.name}
            <br />
            <strong>Status:</strong>{' '}
            <span
              className={`${styles.statusBadge} ${styles[activeCalculation.status]}`}
            >
              {activeCalculation.status}
            </span>
          </div>
        </div>
      </div>
    );
  }

  // Early return if processedData is null (already handled by loading/error states above)
  if (!processedData) {
    return null;
  }

  const { results, parameters } = processedData;

  return (
    <div className={styles.pageContainer}>
      <div className={styles.pageContent}>
        <CalculationSettingsSummary
          activeCalculation={activeCalculation}
          processedData={{
            results,
            parameters,
            shouldShowTDDFTSection: !!processedData.shouldShowTDDFTSection,
            shouldShowCASSection: processedData.shouldShowCASSection,
          }}
        />

        {shouldShowOptimizedStructure && (
          <OptimizedStructureSection
            results={results}
            parameters={parameters}
            currentStyle={currentStyle}
            onStyleChange={setCurrentStyle}
            showAxes={showAxes}
            onShowAxesChange={setShowAxes}
            showCoordinates={showCoordinates}
            onShowCoordinatesChange={setShowCoordinates}
            useAtomicRadii={useAtomicRadii}
            onUseAtomicRadiiChange={setUseAtomicRadii}
          />
        )}

        {processedData.shouldShowElectronicProperties && (
          <ElectronicPropertiesSection
            results={results}
            parameters={parameters}
            activeCalculation={activeCalculation}
          />
        )}

        {processedData.shouldShowCASSection && (
          <CASResultsSection results={results} parameters={parameters} />
        )}

        {processedData.shouldShowEnergeticsSection && (
          <EnergeticsSection
            results={results}
            parameters={parameters}
            processedData={{
              shouldShowCCSDSection: processedData.shouldShowCCSDSection,
            }}
          />
        )}

        {processedData.shouldShowVibrationalSection && (
          <VibrationalAnalysisSection
            calculationId={activeCalculation.id}
            results={results}
            irSpectrumData={irSpectrumData}
            selectedIRPeakIndex={selectedIRPeakIndex}
            selectedVibrationMode={selectedVibrationMode}
            irSettings={irSettings}
            onSetIRSettings={setIRSettings}
            onSpectrumDataLoaded={handleSpectrumDataLoaded}
            onPeakSelect={handleIRPeakSelect}
            onClearSelection={handleClearVibrationSelection}
            onError={setError}
          />
        )}

        <MolecularOrbitalsSection
          calculationId={activeCalculation.id}
          selectedOrbitalIndex={selectedOrbitalIndex}
          onOrbitalSelect={handleOrbitalSelect}
          onError={setError}
        />

        {processedData.shouldShowTDDFTSection &&
          results.excitation_energies && (
            <TDDFTResultsSection results={results} parameters={parameters} />
          )}

        <TechnicalDetailsSection results={results} parameters={parameters} />
      </div>
    </div>
  );
};
