import React from 'react';
import styles from '../../pages/CalculationResultsPage.module.css';
import {
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';

interface TDDFTResultsSectionProps {
  results: CalculationResults;
  parameters: CalculationParameters;
}

export const TDDFTResultsSection = React.memo<TDDFTResultsSectionProps>(
  ({ results, parameters }) => {
    if (!results.excitation_energies) {
      return null;
    }

    const excitationEnergies = results.excitation_energies;

    return (
      <>
        {/* Excited States Summary */}
        <section
          className={`${styles.calculationSection} ${styles.excitedStatesSection}`}
        >
          <h2 className={styles.secondaryHeader}>
            TDDFT Advanced Results - Excited States
          </h2>
          <div className={styles.excitedStatesGrid}>
            <div>
              <strong>Number of States:</strong>{' '}
              <code>{excitationEnergies.length}</code>
            </div>
            <div>
              <strong>TDDFT Method:</strong>{' '}
              <code>{(parameters as any).tddft_method || 'TDDFT'}</code>
            </div>
            <div>
              <strong>Lowest Excitation:</strong>{' '}
              <code>{excitationEnergies[0]?.toFixed(4) || 'N/A'} eV</code>
            </div>
            <div>
              <strong>UV-Vis Range:</strong>{' '}
              <code>{results.excitation_wavelengths?.[0]?.toFixed(0)} nm</code>
            </div>
          </div>
        </section>

        {/* Excitation Energies Table */}
        <section
          className={`${styles.calculationSection} ${styles.excitationTableSection}`}
        >
          <h2>Excitation Energies and Transitions</h2>
          <div className={styles.tableContainer}>
            <table className={styles.dataTable}>
              <thead>
                <tr>
                  <th>State</th>
                  <th className={styles.rightAlign}>Energy (eV)</th>
                  <th className={styles.rightAlign}>Wavelength (nm)</th>
                  <th className={styles.rightAlign}>Osc. Strength</th>
                </tr>
              </thead>
              <tbody>
                {excitationEnergies.map((energy: number, index: number) => {
                  const wavelength = results.excitation_wavelengths?.[index];
                  const oscStrength = results.oscillator_strengths?.[index];

                  return (
                    <tr key={index}>
                      <td>S{index + 1}</td>
                      <td className={`${styles.rightAlign} ${styles.monoFont}`}>
                        {energy.toFixed(4)}
                      </td>
                      <td className={`${styles.rightAlign} ${styles.monoFont}`}>
                        {wavelength ? wavelength.toFixed(1) : 'N/A'}
                      </td>
                      <td className={`${styles.rightAlign} ${styles.monoFont}`}>
                        {oscStrength !== undefined
                          ? oscStrength.toFixed(6)
                          : 'N/A'}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </section>

        {/* UV-Vis Spectrum Visualization */}
        <section
          className={`${styles.calculationSection} ${styles.uvVisSection}`}
        >
          <h2>UV-Vis Spectrum (Simulated)</h2>
          <div className={styles.uvVisChart}>
            <svg width="100%" height="100%" viewBox="0 0 800 300">
              {/* Background Grid */}
              <defs>
                <pattern
                  id="grid"
                  width="40"
                  height="30"
                  patternUnits="userSpaceOnUse"
                >
                  <path
                    d="M 40 0 L 0 0 0 30"
                    fill="none"
                    stroke="#f0f0f0"
                    strokeWidth="1"
                  />
                </pattern>
              </defs>
              <rect width="800" height="300" fill="url(#grid)" />

              {/* Spectrum Bars */}
              {results.excitation_wavelengths?.map(
                (wavelength: number, index: number) => {
                  if (!wavelength || wavelength < 200 || wavelength > 800)
                    return null;

                  const x = ((wavelength - 200) / 600) * 760 + 20;
                  const intensity = results.oscillator_strengths?.[index] || 0;
                  const height = Math.min(intensity * 500, 250);

                  const getColor = (wl: number) => {
                    if (wl < 380) return '#8a2be2';
                    if (wl < 450) return '#4b0082';
                    if (wl < 495) return '#0000ff';
                    if (wl < 570) return '#00ff00';
                    if (wl < 590) return '#ffff00';
                    if (wl < 620) return '#ffa500';
                    if (wl < 750) return '#ff0000';
                    return '#8b4513';
                  };

                  return (
                    <rect
                      key={index}
                      x={x - 1}
                      y={270 - height}
                      width="2"
                      height={height}
                      fill={getColor(wavelength)}
                      opacity="0.7"
                    />
                  );
                }
              )}

              {/* Axis Labels */}
              <text x="20" y="295" fontSize="12" fill="#666">
                200nm
              </text>
              <text x="400" y="295" fontSize="12" fill="#666">
                500nm
              </text>
              <text x="780" y="295" fontSize="12" fill="#666" textAnchor="end">
                800nm
              </text>
              <text
                x="10"
                y="15"
                fontSize="12"
                fill="#666"
                transform="rotate(-90, 10, 15)"
                textAnchor="end"
              >
                Intensity
              </text>
            </svg>
          </div>
          <div className={styles.uvVisDescription}>
            UV-Vis absorption spectrum showing calculated transitions. Colors
            represent approximate wavelength regions.
          </div>
        </section>

        {/* Transition Dipole Moments */}
        {results.transition_dipoles &&
          results.transition_dipoles.length > 0 && (
            <section
              className={`${styles.calculationSection} ${styles.transitionDipoleSection}`}
            >
              <h2>Transition Dipole Moments</h2>
              <div className={styles.tableContainer}>
                <table className={styles.dataTable}>
                  <thead>
                    <tr>
                      <th>State</th>
                      <th className={styles.rightAlign}>μx (a.u.)</th>
                      <th className={styles.rightAlign}>μy (a.u.)</th>
                      <th className={styles.rightAlign}>μz (a.u.)</th>
                      <th className={styles.rightAlign}>|μ| (a.u.)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {results.transition_dipoles.map(
                      (dipole: any, index: number) => {
                        const magnitude = Math.sqrt(
                          dipole.x * dipole.x +
                            dipole.y * dipole.y +
                            dipole.z * dipole.z
                        );

                        return (
                          <tr key={index}>
                            <td>S{index + 1}</td>
                            <td
                              className={`${styles.rightAlign} ${styles.monoFont}`}
                            >
                              {dipole.x.toFixed(6)}
                            </td>
                            <td
                              className={`${styles.rightAlign} ${styles.monoFont}`}
                            >
                              {dipole.y.toFixed(6)}
                            </td>
                            <td
                              className={`${styles.rightAlign} ${styles.monoFont}`}
                            >
                              {dipole.z.toFixed(6)}
                            </td>
                            <td
                              className={`${styles.rightAlign} ${styles.monoFont}`}
                            >
                              <strong>{magnitude.toFixed(6)}</strong>
                            </td>
                          </tr>
                        );
                      }
                    )}
                  </tbody>
                </table>
              </div>
            </section>
          )}

        {/* Natural Transition Orbital Analysis */}
        {results.nto_analysis && results.nto_analysis.length > 0 && (
          <section
            className={`${styles.calculationSection} ${styles.ntoSection}`}
          >
            <h2>Natural Transition Orbital (NTO) Analysis</h2>
            <div className={styles.sectionDescription}>
              NTO analysis provides a more intuitive description of electronic
              excitations by decomposing the transition density matrix into
              dominant hole-particle orbital pairs.
            </div>
            {results.nto_analysis.map((stateData: any, stateIndex: number) => (
              <div key={stateIndex} className={styles.ntoState}>
                <h3 className={styles.ntoStateTitle}>
                  Excited State S{stateData.state} (
                  {stateData.energy?.toFixed(4)} eV)
                </h3>
                {stateData.nto_pairs && stateData.nto_pairs.length > 0 ? (
                  <div className={styles.tableContainer}>
                    <table className={`${styles.dataTable} ${styles.ntoTable}`}>
                      <thead>
                        <tr>
                          <th>NTO Pair</th>
                          <th>Transition</th>
                          <th className={styles.rightAlign}>Weight</th>
                          <th className={styles.rightAlign}>
                            Contribution (%)
                          </th>
                          <th className={styles.centerAlign}>
                            Orbital Indices
                          </th>
                        </tr>
                      </thead>
                      <tbody>
                        {stateData.nto_pairs.map(
                          (pair: any, pairIndex: number) => (
                            <tr key={pairIndex}>
                              <td>
                                <strong>#{pairIndex + 1}</strong>
                              </td>
                              <td className={styles.ntoTransition}>
                                <span className={styles.holeOrbital}>
                                  {pair.hole_orbital}
                                </span>
                                <span className={styles.transitionArrow}>
                                  →
                                </span>
                                <span className={styles.particleOrbital}>
                                  {pair.particle_orbital}
                                </span>
                              </td>
                              <td
                                className={`${styles.rightAlign} ${styles.monoFont}`}
                              >
                                {pair.weight?.toFixed(6) || 'N/A'}
                              </td>
                              <td className={styles.rightAlign}>
                                <div className={styles.ntoContribution}>
                                  <div
                                    className={`${styles.contributionBar} ${
                                      pair.contribution >= 50
                                        ? styles.contributionBarHigh
                                        : pair.contribution >= 25
                                          ? styles.contributionBarMedium
                                          : styles.contributionBarLow
                                    }`}
                                    style={{
                                      width: `${Math.min(pair.contribution || 0, 100)}%`,
                                    }}
                                  />
                                  <span className={styles.contributionValue}>
                                    {pair.contribution?.toFixed(1) || 'N/A'}%
                                  </span>
                                </div>
                              </td>
                              <td
                                className={`${styles.centerAlign} ${styles.ntoIndices}`}
                              >
                                {pair.hole_orbital_index} →{' '}
                                {pair.particle_orbital_index}
                              </td>
                            </tr>
                          )
                        )}
                      </tbody>
                    </table>
                  </div>
                ) : (
                  <div className={styles.ntoNoData}>
                    No significant NTO pairs found for this excited state.
                  </div>
                )}
                <div className={styles.ntoCount}>
                  Total NTO pairs analyzed: {stateData.total_nto_pairs || 0}
                </div>
              </div>
            ))}
            <div className={styles.ntoHelpBox}>
              <h4 className={styles.ntoHelpTitle}>
                💡 How to Read NTO Analysis
              </h4>
              <ul className={styles.ntoHelpList}>
                <li>
                  <strong>Hole Orbitals (red)</strong>: Orbitals from which
                  electrons are excited (mainly HOMO-type)
                </li>
                <li>
                  <strong>Particle Orbitals (blue)</strong>: Orbitals to which
                  electrons are excited (mainly LUMO-type)
                </li>
                <li>
                  <strong>Weight</strong>: Contribution weight of this orbital
                  pair (singular value)
                </li>
                <li>
                  <strong>Contribution</strong>: Percentage contribution of this
                  pair to the total transition
                </li>
                <li>
                  Higher contribution pairs represent the main electronic
                  transitions of the excited state
                </li>
              </ul>
            </div>
          </section>
        )}
      </>
    );
  }
);
