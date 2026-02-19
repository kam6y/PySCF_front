import React from 'react';
import styles from '../../pages/CalculationResultsPage.module.css';
import {
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';

interface EnergeticsSectionProps {
  results: CalculationResults;
  parameters: CalculationParameters;
  processedData: {
    shouldShowCCSDSection: boolean;
  };
}

export const EnergeticsSection = React.memo<EnergeticsSectionProps>(
  ({ results, parameters, processedData }) => {
    return (
      <section
        className={`${styles.calculationSection} ${styles.energeticsSection}`}
      >
        <h2 className={styles.primaryHeader}>Energetics</h2>

        {/* Energy Components */}
        {(results.nuclear_repulsion_energy != null ||
          results.electronic_energy != null) && (
          <div className={styles.propertySubsection}>
            <h3>Energy Components</h3>
            <div className={styles.energyComponentsGrid}>
              {results.scf_energy != null && (
                <div>
                  <strong>Total SCF Energy:</strong>
                  <code>{results.scf_energy.toFixed(8)} hartree</code>
                </div>
              )}
              {results.nuclear_repulsion_energy != null && (
                <div>
                  <strong>Nuclear Repulsion Energy:</strong>
                  <code>
                    {results.nuclear_repulsion_energy.toFixed(8)} hartree
                  </code>
                </div>
              )}
              {results.electronic_energy != null && (
                <div>
                  <strong>Electronic Energy:</strong>
                  <code>{results.electronic_energy.toFixed(8)} hartree</code>
                </div>
              )}
            </div>
          </div>
        )}

        {/* Thermochemistry Subsection */}
        {(results.zero_point_energy != null ||
          results.thermal_energy_298K != null) && (
          <div className={styles.propertySubsection}>
            <h3>Thermochemistry</h3>
            <div className={styles.thermochemicalGrid}>
              {results.zero_point_energy !== undefined &&
                results.zero_point_energy !== null && (
                  <div>
                    <strong>Zero-Point Energy:</strong>
                    <br />
                    <code>{results.zero_point_energy.toFixed(8)} hartree</code>
                  </div>
                )}
              {results.thermal_energy_298K !== undefined &&
                results.thermal_energy_298K !== null && (
                  <div>
                    <strong>Thermal Energy (298.15 K):</strong>
                    <br />
                    <code>
                      {results.thermal_energy_298K.toFixed(8)} hartree
                    </code>
                  </div>
                )}
              {results.entropy_298K !== undefined &&
                results.entropy_298K !== null && (
                  <div>
                    <strong>Entropy (298.15 K):</strong>
                    <br />
                    <code>{results.entropy_298K.toFixed(8)} hartree/K</code>
                  </div>
                )}
              {results.gibbs_free_energy_298K !== undefined &&
                results.gibbs_free_energy_298K !== null && (
                  <div>
                    <strong>Gibbs Free Energy (298.15 K):</strong>
                    <br />
                    <code>
                      {results.gibbs_free_energy_298K.toFixed(8)} hartree
                    </code>
                  </div>
                )}
              {results.heat_capacity_298K !== undefined &&
                results.heat_capacity_298K !== null && (
                  <div>
                    <strong>Heat Capacity (298.15 K):</strong>
                    <br />
                    <code>
                      {results.heat_capacity_298K.toFixed(8)} hartree/K
                    </code>
                  </div>
                )}
            </div>
          </div>
        )}

        {/* MP2 Energetics - Conditional subsection */}
        {parameters.calculation_method === 'MP2' && (
          <div className={styles.propertySubsection}>
            <h3>MP2 Energetics</h3>
            <div className={styles.thermochemicalGrid}>
              <div>
                <strong>HF Energy:</strong>{' '}
                <code>
                  {((results as any).hf_energy || results.scf_energy)?.toFixed(
                    6
                  )}{' '}
                  Hartree
                </code>
              </div>
              <div>
                <strong>MP2 Correlation Energy:</strong>{' '}
                <code>
                  {(results as any).mp2_correlation_energy?.toFixed(6)} Hartree
                </code>
              </div>
              <div>
                <strong>MP2 Total Energy:</strong>{' '}
                <code>
                  {(results as any).mp2_total_energy?.toFixed(6)} Hartree
                </code>
              </div>
            </div>

            {/* Correlation Components */}
            {(results.mp2_same_spin_correlation != null ||
              results.mp2_opposite_spin_correlation != null) && (
              <div className={styles.correlationComponents}>
                <h4 className={styles.subsectionHeader}>
                  Correlation Energy Components
                </h4>
                <div className={styles.thermochemicalGrid}>
                  {results.mp2_same_spin_correlation != null && (
                    <div>
                      <strong>Same-Spin Correlation:</strong>{' '}
                      <code>
                        {results.mp2_same_spin_correlation.toFixed(6)} Hartree
                      </code>
                    </div>
                  )}
                  {results.mp2_opposite_spin_correlation != null && (
                    <div>
                      <strong>Opposite-Spin Correlation:</strong>{' '}
                      <code>
                        {results.mp2_opposite_spin_correlation.toFixed(6)}{' '}
                        Hartree
                      </code>
                    </div>
                  )}
                </div>
                <div className={styles.sectionDescription}>
                  ℹ️ These components provide insight into the nature of
                  electron correlation
                </div>
              </div>
            )}
          </div>
        )}

        {/* CCSD Energetics - Conditional subsection */}
        {processedData.shouldShowCCSDSection && (
          <div className={styles.propertySubsection}>
            <h3>CCSD Energetics</h3>
            <div className={styles.thermochemicalGrid}>
              <div>
                <strong>HF Energy:</strong>{' '}
                <code>
                  {((results as any).hf_energy || results.scf_energy)?.toFixed(
                    6
                  )}{' '}
                  Hartree
                </code>
              </div>
              <div>
                <strong>CCSD Correlation Energy:</strong>{' '}
                <code>
                  {(results as any).ccsd_correlation_energy?.toFixed(6)} Hartree
                </code>
              </div>
              <div>
                <strong>CCSD Total Energy:</strong>{' '}
                <code>
                  {(results as any).ccsd_total_energy?.toFixed(6)} Hartree
                </code>
              </div>
              {parameters.calculation_method === 'CCSD_T' &&
                (results as any).ccsd_t_correction && (
                  <>
                    <div>
                      <strong>CCSD(T) Triples Correction:</strong>{' '}
                      <code>
                        {(results as any).ccsd_t_correction?.toFixed(6)} Hartree
                      </code>
                    </div>
                    <div>
                      <strong>CCSD(T) Total Energy:</strong>{' '}
                      <code>
                        {(results as any).ccsd_t_total_energy?.toFixed(6)}{' '}
                        Hartree
                      </code>
                    </div>
                  </>
                )}
            </div>

            {/* CCSD Diagnostic Indicators */}
            {(results.ccsd_t1_diagnostic != null ||
              results.ccsd_d1_diagnostic != null ||
              results.ccsd_d2_diagnostic != null) && (
              <div className={styles.diagnosticsSection}>
                <h4 className={styles.subsectionHeader}>
                  Diagnostic Indicators
                </h4>
                <div className={styles.diagnosticsGrid}>
                  {results.ccsd_t1_diagnostic != null && (
                    <div className={styles.diagnosticBox}>
                      <strong>T1 Diagnostic:</strong>{' '}
                      <code className={styles.diagnosticValue}>
                        {results.ccsd_t1_diagnostic.toFixed(6)}
                      </code>
                    </div>
                  )}
                  {results.ccsd_d1_diagnostic != null && (
                    <div className={styles.diagnosticBox}>
                      <strong>D1 Diagnostic:</strong>{' '}
                      <code className={styles.diagnosticValue}>
                        {results.ccsd_d1_diagnostic.toFixed(6)}
                      </code>
                    </div>
                  )}
                  {results.ccsd_d2_diagnostic != null && (
                    <div className={styles.diagnosticBox}>
                      <strong>D2 Diagnostic:</strong>{' '}
                      <code className={styles.diagnosticValue}>
                        {results.ccsd_d2_diagnostic.toFixed(6)}
                      </code>
                    </div>
                  )}
                </div>
                <div className={styles.referenceInfo}>
                  <h4>Diagnostic Reference Values</h4>
                  <ul>
                    <li>
                      T1: Values &gt; 0.02 may indicate multi-reference
                      character
                    </li>
                    <li>
                      D1: Values &gt; 0.05 may indicate open-shell character
                    </li>
                    <li>
                      D2: Values &gt; 0.15 may indicate strong correlation
                      effects
                    </li>
                  </ul>
                </div>
              </div>
            )}

            {(results as any).frozen_core && (
              <div className={styles.sectionDescription}>
                ℹ️ Frozen core approximation was used in this calculation
              </div>
            )}
          </div>
        )}

        {/* Frontier Orbital Analysis */}
        {(results.homo_energy_ev != null ||
          results.lumo_energy_ev != null ||
          results.homo_lumo_gap_ev != null) && (
          <div className={styles.propertySubsection}>
            <h3>Frontier Orbital Analysis</h3>
            <div className={styles.frontierOrbitalsGrid}>
              {results.homo_index != null && (
                <div className={styles.orbitalEnergyBox}>
                  <strong>HOMO Index:</strong>
                  <code>{results.homo_index}</code>
                </div>
              )}
              {results.lumo_index != null && (
                <div className={styles.orbitalEnergyBox}>
                  <strong>LUMO Index:</strong>
                  <code>{results.lumo_index}</code>
                </div>
              )}
              {results.num_occupied_orbitals != null && (
                <div className={styles.orbitalEnergyBox}>
                  <strong>Occupied Orbitals:</strong>
                  <code>{results.num_occupied_orbitals}</code>
                </div>
              )}
              {results.num_virtual_orbitals != null && (
                <div className={styles.orbitalEnergyBox}>
                  <strong>Virtual Orbitals:</strong>
                  <code>{results.num_virtual_orbitals}</code>
                </div>
              )}
              {results.homo_energy_ev != null && (
                <div className={styles.orbitalEnergyBox}>
                  <strong>HOMO Energy:</strong>
                  <div className={styles.energyValue}>
                    <code className={styles.primaryValue}>
                      {results.homo_energy_ev.toFixed(4)} eV
                    </code>
                    <code className={styles.secondaryUnit}>
                      ({results.homo_energy_hartree?.toFixed(6)} hartree)
                    </code>
                  </div>
                </div>
              )}
              {results.lumo_energy_ev != null && (
                <div className={styles.orbitalEnergyBox}>
                  <strong>LUMO Energy:</strong>
                  <div className={styles.energyValue}>
                    <code className={styles.primaryValue}>
                      {results.lumo_energy_ev.toFixed(4)} eV
                    </code>
                    <code className={styles.secondaryUnit}>
                      ({results.lumo_energy_hartree?.toFixed(6)} hartree)
                    </code>
                  </div>
                </div>
              )}
              {results.homo_lumo_gap_ev != null && (
                <div className={styles.gapBox}>
                  <strong>HOMO-LUMO Gap:</strong>
                  <div className={styles.gapValue}>
                    <code className={styles.primaryGap}>
                      {results.homo_lumo_gap_ev.toFixed(4)} eV
                    </code>
                    <code className={styles.secondaryUnit}>
                      ({results.homo_lumo_gap_hartree?.toFixed(6)} hartree)
                    </code>
                  </div>
                </div>
              )}
            </div>
          </div>
        )}
      </section>
    );
  }
);
