import React from 'react';
import { components } from '../types/generated-api';

type EnhancedCIAnalysis = components['schemas']['EnhancedCIAnalysis'];
import styles from './CIAnalysisViewer.module.css';

interface CIAnalysisViewerProps {
  enhancedCiAnalysis?: EnhancedCIAnalysis | null;
  ciCoefficientsAvailable?: boolean;
  kernelReturnInfo?: Record<string, any> | null;
}

export const CIAnalysisViewer: React.FC<CIAnalysisViewerProps> = ({
  enhancedCiAnalysis,
  ciCoefficientsAvailable,
  kernelReturnInfo,
}) => {
  if (!enhancedCiAnalysis || !enhancedCiAnalysis.available) {
    return (
      <div className={styles.unavailable}>
        <h4>CI Coefficient Analysis</h4>
        <p className={styles.unavailableMessage}>
          {enhancedCiAnalysis?.reason ||
            (ciCoefficientsAvailable === false
              ? 'No CI coefficients available from calculation'
              : 'CI analysis not available')}
        </p>
        {kernelReturnInfo && (
          <div className={styles.kernelInfo}>
            <h5>Kernel Return Information</h5>
            <pre className={styles.kernelInfoPre}>
              {JSON.stringify(kernelReturnInfo, null, 2)}
            </pre>
          </div>
        )}
      </div>
    );
  }

  if (enhancedCiAnalysis.error) {
    return (
      <div className={styles.error}>
        <h4>CI Coefficient Analysis</h4>
        <p className={styles.errorMessage}>Error: {enhancedCiAnalysis.error}</p>
      </div>
    );
  }

  const getCharacterDescription = (character: string) => {
    switch (character) {
      case 'single_configuration':
        return 'Single configuration dominates the wavefunction';
      case 'few_configuration':
        return 'Few configurations contribute significantly';
      case 'multiconfigurational':
        return 'Many configurations contribute to the wavefunction';
      default:
        return 'Unknown wavefunction character';
    }
  };

  const getCharacterClass = (character: string) => {
    switch (character) {
      case 'single_configuration':
        return styles.singleConfig;
      case 'few_configuration':
        return styles.fewConfig;
      case 'multiconfigurational':
        return styles.multiConfig;
      default:
        return styles.unknownConfig;
    }
  };

  return (
    <div className={styles.container}>
      <h4>CI Coefficient Analysis</h4>

      {/* Overview Section */}
      <div className={styles.overview}>
        <div className={styles.overviewItem}>
          <span className={styles.label}>Total Coefficients:</span>
          <span className={styles.value}>
            {enhancedCiAnalysis.total_coefficients?.toLocaleString()}
          </span>
        </div>
        <div className={styles.overviewItem}>
          <span className={styles.label}>CI Vector Shape:</span>
          <span className={styles.value}>
            {enhancedCiAnalysis.ci_vector_shape}
          </span>
        </div>
        <div className={styles.overviewItem}>
          <span className={styles.label}>Effective Configurations:</span>
          <span className={styles.value}>
            {enhancedCiAnalysis.effective_configurations}
          </span>
        </div>
      </div>

      {/* Wavefunction Character */}
      {enhancedCiAnalysis.wavefunction_character && (
        <div className={styles.wavefunctionCharacter}>
          <h5>Wavefunction Character</h5>
          <div
            className={`${styles.characterBadge} ${getCharacterClass(enhancedCiAnalysis.wavefunction_character)}`}
          >
            {enhancedCiAnalysis.wavefunction_character
              .replace('_', ' ')
              .toUpperCase()}
          </div>
          <p className={styles.characterDescription}>
            {getCharacterDescription(enhancedCiAnalysis.wavefunction_character)}
          </p>
        </div>
      )}

      {/* Key Statistics */}
      <div className={styles.statistics}>
        <h5>Key Statistics</h5>
        <div className={styles.statGrid}>
          <div className={styles.statItem}>
            <span className={styles.statLabel}>Leading Contribution:</span>
            <span className={styles.statValue}>
              {enhancedCiAnalysis.leading_contribution_percent?.toFixed(1)}%
            </span>
          </div>
          <div className={styles.statItem}>
            <span className={styles.statLabel}>
              Multiconfigurational Character:
            </span>
            <span className={styles.statValue}>
              {enhancedCiAnalysis.multiconfigurational_character?.toFixed(1)}%
            </span>
          </div>
          <div className={styles.statItem}>
            <span className={styles.statLabel}>Wavefunction Entropy:</span>
            <span className={styles.statValue}>
              {enhancedCiAnalysis.wavefunction_entropy?.toFixed(4)}
            </span>
          </div>
          <div className={styles.statItem}>
            <span className={styles.statLabel}>Normalization:</span>
            <span className={styles.statValue}>
              {enhancedCiAnalysis.normalization?.toFixed(6)}
            </span>
          </div>
        </div>
      </div>

      {/* Major Configurations */}
      {enhancedCiAnalysis.major_configurations &&
        enhancedCiAnalysis.major_configurations.length > 0 && (
          <div className={styles.majorConfigurations}>
            <h5>Major Configurations</h5>
            <div className={styles.configTable}>
              <div className={styles.configTableHeader}>
                <span>Index</span>
                <span>Coefficient</span>
                <span>Contribution (%)</span>
                <span>Cumulative (%)</span>
              </div>
              {enhancedCiAnalysis.major_configurations
                .slice(0, 20)
                .map((config: any, index: number) => (
                  <div
                    key={config.configuration_index}
                    className={styles.configTableRow}
                  >
                    <span className={styles.configIndex}>
                      {config.configuration_index}
                    </span>
                    <span className={styles.configCoeff}>
                      {config.coefficient?.toFixed(6)}
                    </span>
                    <span className={styles.configContrib}>
                      {config.contribution_percent?.toFixed(2)}
                    </span>
                    <span className={styles.configCumulative}>
                      {config.cumulative_percent?.toFixed(2)}
                    </span>
                  </div>
                ))}
              {enhancedCiAnalysis.major_configurations.length > 20 && (
                <div className={styles.configTableNote}>
                  ... and {enhancedCiAnalysis.major_configurations.length - 20}{' '}
                  more configurations
                </div>
              )}
            </div>
          </div>
        )}

      {/* Technical Details */}
      {kernelReturnInfo && (
        <div className={styles.technicalDetails}>
          <h5>Technical Details</h5>
          <div className={styles.detailsGrid}>
            <div className={styles.detailItem}>
              <span className={styles.detailLabel}>Source:</span>
              <span className={styles.detailValue}>
                {enhancedCiAnalysis.source}
              </span>
            </div>
            {kernelReturnInfo.tuple_length && (
              <div className={styles.detailItem}>
                <span className={styles.detailLabel}>
                  Kernel Return Elements:
                </span>
                <span className={styles.detailValue}>
                  {kernelReturnInfo.tuple_length}
                </span>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
};
