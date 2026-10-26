import React, { useMemo } from 'react';
import styles from '../../pages/CalculationResultsPage.module.css';
import { MoleculeViewerSection } from '../MoleculeViewerSection';
import {
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';
import { StyleSpec } from '../../../types/3dmol';

interface OptimizedStructureSectionProps {
  results: CalculationResults;
  parameters: CalculationParameters;
  currentStyle: StyleSpec | null;
  onStyleChange: (style: StyleSpec) => void;
  showAxes: boolean;
  onShowAxesChange: (show: boolean) => void;
  showCoordinates: boolean;
  onShowCoordinatesChange: (show: boolean) => void;
  showAtomNumbers: boolean;
  onShowAtomNumbersChange: (show: boolean) => void;
  useAtomicRadii: boolean;
  onUseAtomicRadiiChange: (use: boolean) => void;
  selectedAtomIndices: number[];
  onAtomClick: (atomIndex: number) => void;
}

export const OptimizedStructureSection =
  React.memo<OptimizedStructureSectionProps>(
    ({
      results,
      parameters,
      currentStyle,
      onStyleChange,
      showAxes,
      onShowAxesChange,
      showCoordinates,
      onShowCoordinatesChange,
      showAtomNumbers,
      onShowAtomNumbersChange,
      useAtomicRadii,
      onUseAtomicRadiiChange,
      selectedAtomIndices,
      onAtomClick,
    }) => {
      const structureTitle = useMemo(() => {
        switch (parameters.calculation_method) {
          case 'HF':
            return 'HF-Optimized Molecular Structure';
          case 'MP2':
            return 'MP2-Optimized Molecular Structure';
          case 'DFT':
          default:
            return 'DFT-Optimized Molecular Structure';
        }
      }, [parameters.calculation_method]);

      return (
        <section
          className={`${styles.calculationSection} ${styles.structureSection}`}
        >
          <h2 className={styles.primaryHeader}>{structureTitle}</h2>

          {/* 2-Column Layout: Left (Info + Coordinates) and Right (3D Viewer) */}
          <div className={styles.structureContentWrapper}>
            {/* Left Column: Description, Molecular Info, and XYZ Coordinates */}
            <div className={styles.structureLeftColumn}>
              {/* Frequency Quality Indicators */}
              {results.frequency_analysis_performed &&
                results.imaginary_frequencies_count != null &&
                results.imaginary_frequencies_count >= 0 && (
                  <>
                    {/* Frequency Data */}
                    <div className={styles.frequencyStatus}>
                      <div>
                        <strong>Imaginary Frequencies:</strong>{' '}
                        <code>{results.imaginary_frequencies_count}</code>
                      </div>
                    </div>

                    {/* Imaginary Frequencies Warning */}
                    {results.imaginary_frequencies_count > 0 && (
                      <div className={styles.imaginaryFrequencyWarning}>
                        <strong>⚠️ Optimization Quality Warning:</strong>
                        <div className={styles.warningContent}>
                          This structure has{' '}
                          {results.imaginary_frequencies_count} imaginary
                          {results.imaginary_frequencies_count === 1
                            ? ' frequency'
                            : ' frequencies'}
                          , which may indicate:
                          <ul>
                            <li>
                              The structure is at a transition state or saddle
                              point
                            </li>
                            <li>
                              The optimization did not fully converge to a
                              minimum
                            </li>
                            <li>Further optimization may be needed</li>
                          </ul>
                        </div>
                      </div>
                    )}
                  </>
                )}

              {/* XYZ Coordinates Display */}
              <div className={styles.xyzCoordinatesContainer}>
                <strong>XYZ Coordinates:</strong>
                <pre className={styles.xyzCoordinates}>
                  {results.optimized_geometry}
                </pre>
              </div>
            </div>

            {/* Right Column: 3D Molecular Viewer */}
            <div className={styles.structureRightColumn}>
              <div className={styles.viewer3DContainer}>
                <h3>3D Molecular Visualization</h3>
                <MoleculeViewerSection
                  hasValidMolecule={!!results.optimized_geometry}
                  xyzData={results.optimized_geometry}
                  currentStyle={currentStyle}
                  onStyleChange={onStyleChange}
                  showAxes={showAxes}
                  onShowAxesChange={onShowAxesChange}
                  showCoordinates={showCoordinates}
                  onShowCoordinatesChange={onShowCoordinatesChange}
                  showAtomNumbers={showAtomNumbers}
                  onShowAtomNumbersChange={onShowAtomNumbersChange}
                  useAtomicRadii={useAtomicRadii}
                  onUseAtomicRadiiChange={onUseAtomicRadiiChange}
                  selectedAtomIndices={selectedAtomIndices}
                  onAtomClick={onAtomClick}
                />
              </div>
            </div>
          </div>
        </section>
      );
    }
  );
