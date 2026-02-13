import styles from '../../pages/CalculationResultsPage.module.css';
import {
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';
import { CIAnalysisViewer } from '../CIAnalysisViewer';

interface CASResultsSectionProps {
  results: CalculationResults;
  parameters: CalculationParameters;
}

export const CASResultsSection = ({
  results,
  parameters,
}: CASResultsSectionProps) => {
  return (
<section
  className={`${styles.calculationSection} ${styles.casSection}`}
>
  <h2 className={styles.secondaryHeader}>
    {parameters.calculation_method === 'CASCI' ? 'CASCI' : 'CASSCF'}{' '}
    Advanced Results
  </h2>

  {/* Energy Information */}
  <div className={styles.casEnergyGrid}>
    <div>
      <strong>SCF Reference Energy:</strong>{' '}
      <code>{results.scf_energy?.toFixed(8) || 'N/A'} hartree</code>
    </div>
    {results.casci_energy && (
      <div>
        <strong>CASCI Energy:</strong>{' '}
        <code>{results.casci_energy.toFixed(8)} hartree</code>
      </div>
    )}
    {results.casscf_energy && (
      <div>
        <strong>CASSCF Energy:</strong>{' '}
        <code>{results.casscf_energy.toFixed(8)} hartree</code>
      </div>
    )}
    {results.correlation_energy && (
      <div>
        <strong>Correlation Energy:</strong>{' '}
        <code>{results.correlation_energy.toFixed(8)} hartree</code>
        <div className={styles.metaText}>
          ({parameters.calculation_method} - SCF)
        </div>
      </div>
    )}
  </div>

  {/* Active Space Information */}
  <div className={styles.activeSpaceInfo}>
    <h3>Active Space Configuration</h3>
    <div className={styles.activeSpaceGrid}>
      <div>
        <strong>Active Orbitals:</strong>{' '}
        <code>{(results as any).ncas || 'N/A'}</code>
      </div>
      <div>
        <strong>Active Electrons:</strong>{' '}
        <code>{(results as any).nelecas || 'N/A'}</code>
      </div>
      {(results as any).method && (
        <div>
          <strong>Reference Method:</strong>{' '}
          <code>{(results as any).method}</code>
        </div>
      )}
    </div>
  </div>

  {/* Convergence Information (CASSCF only) */}
  {parameters.calculation_method === 'CASSCF' && (
    <div className={styles.convergenceInfo}>
      <h3>Convergence Information</h3>
      <div className={styles.convergenceGrid}>
        <div>
          <strong>Converged:</strong>{' '}
          <span
            className={
              results.converged
                ? styles.convergedTrue
                : styles.convergedFalse
            }
          >
            {results.converged ? '✅ Yes' : '❌ No'}
          </span>
        </div>
        {results.macro_iterations !== undefined && (
          <div>
            <strong>Macro Iterations:</strong>{' '}
            <code>{results.macro_iterations}</code>
          </div>
        )}
      </div>
    </div>
  )}

  {/* Natural Orbital Analysis */}
  {(results as any).natural_orbital_analysis &&
    (results as any).natural_orbital_analysis.enabled && (
      <div className={styles.analysisSection}>
        <h3>Natural Orbital Analysis</h3>
        <div className={styles.analysisDescription}>
          Natural orbitals provide insight into the
          multi-configurational character of the wavefunction through
          orbital occupation patterns.
        </div>

        {(results as any).natural_orbital_analysis
          .occupation_numbers && (
          <div className={styles.orbitalOccupationSection}>
            <h4>Orbital Occupation Numbers</h4>
            <div className={styles.occupationGrid}>
              {(
                results as any
              ).natural_orbital_analysis.occupation_numbers.map(
                (occ: number, index: number) => (
                  <div key={index} className={styles.occupationItem}>
                    <div className={styles.orbitalIndex}>
                      {index + 1}
                    </div>
                    <div className={styles.occupationBar}>
                      <div
                        className={styles.occupationFill}
                        style={{
                          width: `${Math.min((occ / 2.0) * 100, 100)}%`,
                          backgroundColor:
                            occ > 1.5
                              ? '#4caf50'
                              : occ > 0.1
                                ? '#ff9800'
                                : '#f44336',
                        }}
                      />
                    </div>
                    <div className={styles.occupationValue}>
                      {occ.toFixed(3)}
                    </div>
                  </div>
                )
              )}
            </div>

            <div className={styles.occupationSummary}>
              <div className={styles.occupationStats}>
                <div>
                  <strong>Strongly Occupied:</strong>{' '}
                  <code>
                    {(results as any).natural_orbital_analysis
                      .strongly_occupied_count || 0}
                  </code>{' '}
                  (occupation &gt; 1.5)
                </div>
                <div>
                  <strong>Weakly Occupied:</strong>{' '}
                  <code>
                    {(results as any).natural_orbital_analysis
                      .weakly_occupied_count || 0}
                  </code>{' '}
                  (0.1 &lt; occupation &le; 1.5)
                </div>
                <div>
                  <strong>Virtual:</strong>{' '}
                  <code>
                    {(results as any).natural_orbital_analysis
                      .virtual_count || 0}
                  </code>{' '}
                  (occupation &le; 0.1)
                </div>
              </div>

              <div className={styles.electronAnalysis}>
                <div>
                  <strong>Total Active Electrons:</strong>{' '}
                  <code>
                    {(
                      results as any
                    ).natural_orbital_analysis.total_active_electrons?.toFixed(
                      2
                    ) || 'N/A'}
                  </code>
                </div>
                <div>
                  <strong>Effective Electron Pairs:</strong>{' '}
                  <code>
                    {(
                      results as any
                    ).natural_orbital_analysis.effective_electron_pairs?.toFixed(
                      2
                    ) || 'N/A'}
                  </code>
                </div>
                <div>
                  <strong>Effective Unpaired Electrons:</strong>{' '}
                  <code>
                    {(
                      results as any
                    ).natural_orbital_analysis.effective_unpaired_electrons?.toFixed(
                      2
                    ) || 'N/A'}
                  </code>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    )}

  {/* CI Coefficient Analysis */}
  {(results as any).ci_coefficient_analysis &&
    (results as any).ci_coefficient_analysis.available && (
      <div className={styles.analysisSection}>
        <h3>CI Coefficient Analysis</h3>
        <div className={styles.analysisDescription}>
          Analysis of configuration interaction coefficients reveals
          the relative importance of different electronic
          configurations.
        </div>

        <div className={styles.ciSummary}>
          <div className={styles.ciStats}>
            <div>
              <strong>Leading Configuration:</strong>{' '}
              <code>
                {(
                  results as any
                ).ci_coefficient_analysis.leading_contribution_percent?.toFixed(
                  1
                ) || 'N/A'}
                %
              </code>
            </div>
            <div>
              <strong>Multi-configurational Character:</strong>{' '}
              <code>
                {(
                  results as any
                ).ci_coefficient_analysis.multiconfigurational_character?.toFixed(
                  1
                ) || 'N/A'}
                %
              </code>
            </div>
            <div>
              <strong>Total Configurations:</strong>{' '}
              <code>
                {(results as any).ci_coefficient_analysis
                  .total_configurations || 'N/A'}
              </code>
            </div>
          </div>
        </div>

        {(results as any).ci_coefficient_analysis
          .major_configurations &&
          (results as any).ci_coefficient_analysis
            .major_configurations.length > 0 && (
            <div className={styles.majorConfigsSection}>
              <h4>Major Configurations (&gt; 1% contribution)</h4>
              <div className={styles.tableContainer}>
                <table className={styles.dataTable}>
                  <thead>
                    <tr>
                      <th>Config #</th>
                      <th className={styles.rightAlign}>
                        Coefficient
                      </th>
                      <th className={styles.rightAlign}>
                        Contribution (%)
                      </th>
                      <th className={styles.rightAlign}>
                        Cumulative (%)
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    {(
                      results as any
                    ).ci_coefficient_analysis.major_configurations
                      .slice(0, 10)
                      .map((config: any, index: number) => (
                        <tr key={index}>
                          <td>{config.configuration_index + 1}</td>
                          <td className={styles.rightAlign}>
                            <code>
                              {config.coefficient.toFixed(4)}
                            </code>
                          </td>
                          <td className={styles.rightAlign}>
                            <code>
                              {config.contribution_percent.toFixed(2)}
                            </code>
                          </td>
                          <td className={styles.rightAlign}>
                            <code>
                              {config.cumulative_percent.toFixed(2)}
                            </code>
                          </td>
                        </tr>
                      ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
      </div>
    )}

  {/* Enhanced CI Analysis from Kernel Return */}
  {(results as any).enhanced_ci_analysis && (
    <div className={styles.analysisSection}>
      <h3>Enhanced CI Coefficient Analysis</h3>
      <div className={styles.analysisDescription}>
        Detailed analysis of CI coefficients extracted directly from
        the kernel return value, providing comprehensive insight into
        wavefunction composition and multiconfigurational character.
      </div>
      <CIAnalysisViewer
        enhancedCiAnalysis={(results as any).enhanced_ci_analysis}
        ciCoefficientsAvailable={
          (results as any).ci_coefficients_available
        }
        kernelReturnInfo={(results as any).kernel_return_info}
      />
    </div>
  )}

  {/* Orbital Analysis */}
  {(results as any).orbital_overlap_analysis &&
    (results as any).orbital_overlap_analysis.available && (
      <div className={styles.analysisSection}>
        <h3>Active Space Orbital Analysis</h3>
        <div className={styles.analysisDescription}>
          Analysis of orbital overlaps and transformations between SCF
          reference and {parameters.calculation_method} optimized
          orbitals.
        </div>

        <div className={styles.orbitalOverlapSummary}>
          <div className={styles.overlapStats}>
            <div>
              <strong>Active Space Size:</strong>{' '}
              <code>
                {(results as any).orbital_overlap_analysis
                  .active_space_orbitals || 'N/A'}
              </code>
            </div>
            <div>
              <strong>Average Max Overlap:</strong>{' '}
              <code>
                {(
                  results as any
                ).orbital_overlap_analysis.average_max_overlap?.toFixed(
                  3
                ) || 'N/A'}
              </code>
            </div>
            <div>
              <strong>Orbital Transformation:</strong>{' '}
              <span
                className={
                  (results as any).orbital_overlap_analysis
                    .orbital_transformation_character === 'minimal'
                    ? styles.minimalTransformation
                    : styles.significantTransformation
                }
              >
                {(results as any).orbital_overlap_analysis
                  .orbital_transformation_character || 'N/A'}
              </span>
            </div>
          </div>
        </div>

        {(results as any).orbital_overlap_analysis
          .active_orbital_analysis && (
          <div className={styles.activeOrbitalSection}>
            <h4>Active Orbital Character</h4>
            <div className={styles.tableContainer}>
              <table className={styles.dataTable}>
                <thead>
                  <tr>
                    <th>Active Orbital</th>
                    <th>Dominant SCF Orbital</th>
                    <th className={styles.rightAlign}>Max Overlap</th>
                    <th>SCF Orbital Type</th>
                  </tr>
                </thead>
                <tbody>
                  {(
                    results as any
                  ).orbital_overlap_analysis.active_orbital_analysis.map(
                    (orbital: any, index: number) => (
                      <tr key={index}>
                        <td>{orbital.active_orbital_index + 1}</td>
                        <td>{orbital.dominant_scf_orbital + 1}</td>
                        <td className={styles.rightAlign}>
                          <code>
                            {orbital.max_overlap.toFixed(3)}
                          </code>
                        </td>
                        <td>
                          <span
                            className={`${styles.orbitalType} ${styles[orbital.scf_orbital_type]}`}
                          >
                            {orbital.scf_orbital_type}
                          </span>
                        </td>
                      </tr>
                    )
                  )}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
    )}

  {/* CASSCF Orbital Rotation Analysis */}
  {parameters.calculation_method === 'CASSCF' &&
    (results as any).orbital_rotation_analysis &&
    (results as any).orbital_rotation_analysis.available && (
      <div className={styles.analysisSection}>
        <h3>Orbital Rotation Analysis</h3>
        <div className={styles.analysisDescription}>
          Analysis of orbital rotations that occurred during CASSCF
          optimization, showing how much the orbitals changed from the
          initial SCF guess.
        </div>

        <div className={styles.rotationSummary}>
          <div className={styles.rotationStats}>
            <div>
              <strong>Overall Rotation Extent:</strong>{' '}
              <span
                className={`${styles.rotationExtent} ${styles[(results as any).orbital_rotation_analysis.rotation_extent]}`}
              >
                {(results as any).orbital_rotation_analysis
                  .rotation_extent || 'N/A'}
              </span>
            </div>
            <div>
              <strong>Max Rotation Magnitude:</strong>{' '}
              <code>
                {(
                  results as any
                ).orbital_rotation_analysis.overall_rotation_magnitude?.toFixed(
                  3
                ) || 'N/A'}
              </code>
            </div>
          </div>
        </div>

        <div className={styles.orbitalSpaceRotations}>
          <h4>Rotation by Orbital Space</h4>
          <div className={styles.rotationGrid}>
            <div>
              <strong>Core Orbitals:</strong>{' '}
              <code>
                {(
                  results as any
                ).orbital_rotation_analysis.core_orbital_rotation_magnitude?.toFixed(
                  3
                ) || 'N/A'}
              </code>
            </div>
            <div>
              <strong>Active Orbitals:</strong>{' '}
              <code>
                {(
                  results as any
                ).orbital_rotation_analysis.active_orbital_rotation_magnitude?.toFixed(
                  3
                ) || 'N/A'}
              </code>
            </div>
            <div>
              <strong>Virtual Orbitals:</strong>{' '}
              <code>
                {(
                  results as any
                ).orbital_rotation_analysis.virtual_orbital_rotation_magnitude?.toFixed(
                  3
                ) || 'N/A'}
              </code>
            </div>
          </div>
        </div>
      </div>
    )}
</section>
  );
};
