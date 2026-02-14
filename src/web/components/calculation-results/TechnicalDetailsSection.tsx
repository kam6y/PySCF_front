import React from 'react';
import styles from '../../pages/CalculationResultsPage.module.css';
import {
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';

interface TechnicalDetailsSectionProps {
  results: CalculationResults;
  parameters: CalculationParameters;
}

export const TechnicalDetailsSection = React.memo<TechnicalDetailsSectionProps>(({
  results,
  parameters,
}) => {
  return (
<section
  className={`${styles.calculationSection} ${styles.technicalDetailsSection}`}
>
  <h2 className={styles.primaryHeader}>Technical Details</h2>

  {/* SCF Convergence Information */}
  {(results.scf_iterations != null ||
    results.final_energy_change != null ||
    results.final_density_change != null) && (
    <div className={styles.technicalSubsection}>
      <h3>SCF Convergence Information</h3>
      <div className={styles.convergenceGrid}>
        {results.scf_iterations != null && (
          <div>
            <strong>SCF Iterations:</strong>{' '}
            <code>{results.scf_iterations}</code>
          </div>
        )}
        {results.final_energy_change != null && (
          <div>
            <strong>Final Energy Change:</strong>{' '}
            <code>{results.final_energy_change.toExponential(4)}</code>
          </div>
        )}
        {results.final_density_change != null && (
          <div>
            <strong>Final Density Change:</strong>{' '}
            <code>{results.final_density_change.toExponential(4)}</code>
          </div>
        )}
        {results.max_cycle != null && (
          <div>
            <strong>Max SCF Cycles (Setting):</strong>{' '}
            <code>{results.max_cycle}</code>
          </div>
        )}
      </div>
    </div>
  )}

  {/* Calculation Parameters */}
  {!results.gpu_enabled && (
    <div className={styles.technicalSubsection}>
      <h3>Calculation Parameters</h3>
      <div className={styles.parametersGrid}>
        <div>
          <strong>CPU Cores:</strong>{' '}
          {parameters.cpu_cores || 'Default'}
        </div>
        <div>
          <strong>Memory:</strong>{' '}
          {parameters.memory_mb
            ? `${parameters.memory_mb} MB`
            : 'Default'}
        </div>
      </div>
    </div>
  )}

  {/* Checkpoint File Information */}
  <div className={styles.technicalSubsection}>
    <h3>Checkpoint File Information</h3>
    <div className={styles.checkpointInfo}>
      <strong>Working Directory:</strong>{' '}
      <code>{results.working_directory}</code>
    </div>
    <div className={styles.checkpointInfo}>
      <strong>Checkpoint File:</strong>{' '}
      <code>{results.checkpoint_file}</code>
    </div>
    <div className={styles.checkpointStatus}>
      <strong>File Status:</strong>{' '}
      {results.checkpoint_exists
        ? '✅ File exists'
        : '❌ File not found'}
    </div>
    {results.checkpoint_exists && (
      <div className={styles.infoBox}>
        <h4>📁 File Access Methods</h4>
        <div>
          <p>
            <strong>Finder:</strong> Open {results.working_directory} in
            Finder
          </p>
          <p>
            <strong>Terminal:</strong>{' '}
            <code>cd {results.working_directory}</code>
          </p>
          <p>
            <strong>Checkpoint File:</strong>{' '}
            <code>calculation.chk</code>
          </p>
          <p className={styles.sectionDescription}>
            Note: This directory contains molecular orbital data and
            wave function information
          </p>
        </div>
      </div>
    )}
  </div>
</section>
  );
});
