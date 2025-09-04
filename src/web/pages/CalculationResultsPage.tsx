import { useEffect, useState } from 'react';
import styles from './CalculationResultsPage.module.css';
import { CalculationInstance } from '../types/api-types';
import { MolecularOrbitalViewer } from '../components/MolecularOrbitalViewer';
import { MolecularOrbitalEnergyDiagram } from '../components/MolecularOrbitalEnergyDiagram';
import { IRSpectrumViewer } from '../components/IRSpectrumViewer';

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

  useEffect(() => {
    setError(detailsError);
  }, [detailsError]);

  // このuseEffectは不要になったので削除します
  // useEffect(() => {
  //   if (activeCalculation && activeCalculation.status === 'running') {
  //     startPolling();
  //   } else {
  //     stopPolling();
  //   }
  // }, [activeCalculation, startPolling, stopPolling]);

  // Show loading state
  if (isLoadingDetails) {
    return (
      <div className={styles.pageContainer}>
        <div className={styles.pageContent}>
          <h1>Calculation Results</h1>
          <div className={styles.loadingContainer}>
            <div className={styles.loadingText}>
              ⚛️ Loading calculation details...
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
          <h1>Calculation Results</h1>
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
          <h1>Calculation Results</h1>
          <div className={styles.noCalculationContainer}>
            📊 No calculation selected. Please select a calculation from the
            sidebar to view its results.
          </div>
        </div>
      </div>
    );
  }

  // Show message for incomplete calculations
  if (activeCalculation.status !== 'completed' || !activeCalculation.results) {
    const statusMessages = {
      pending:
        '⏳ This calculation is pending. Please run the calculation first.',
      running:
        '⚛️ This calculation is currently running. Please wait for completion.',
      error:
        '❌ This calculation failed. Please check the settings and try again.',
    };

    return (
      <div className={styles.pageContainer}>
        <div className={styles.pageContent}>
          <h1>Calculation Results</h1>
          <div className={styles.statusMessageContainer}>
            {statusMessages[
              activeCalculation.status as keyof typeof statusMessages
            ] || '❓ Calculation results are not available.'}
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

  const results = activeCalculation.results;
  const parameters = activeCalculation.parameters;
  const completedAt = activeCalculation.updatedAt;

  return (
    <div className="page-container">
      <div className="page-content">
        <h1>Quantum Chemistry Calculation Results</h1>

        {/* Calculation Summary */}
        <section
          className={`${styles.calculationSection} ${styles.summarySection}`}
        >
          <h2>Calculation Summary</h2>
          <div className={styles.summaryGrid}>
            <div>
              <strong>Molecule:</strong> {(parameters as any).name || 'Unknown'}
            </div>
            <div>
              <strong>Method:</strong> {parameters.calculation_method}
            </div>
            <div>
              <strong>Basis Set:</strong> {results.basis}
            </div>
            <div>
              <strong>XC Functional:</strong> {results.xc_functional}
            </div>
            <div>
              <strong>Charge:</strong> {results.charge}
            </div>
            <div>
              <strong>Spin Multiplicity:</strong> {results.spin_multiplicity}
            </div>
            <div>
              <strong>Completed:</strong>{' '}
              {new Date(completedAt).toLocaleString()}
            </div>
            <div>
              <strong>Convergence:</strong>{' '}
              {results.converged ? '✅ Converged' : '❌ Not Converged'}
            </div>
          </div>
        </section>

        {/* Energy Results */}
        <section
          className={`${styles.calculationSection} ${styles.energySection}`}
        >
          <h2>Energy Results</h2>
          <div className={styles.energyResult}>
            <strong>SCF Energy:</strong>{' '}
            <code>{results.scf_energy?.toFixed(8) || 'N/A'} hartree</code>
          </div>
        </section>

        {/* Vibrational Frequency Analysis */}
        {results.frequency_analysis_performed && (
          <section
            className={`${styles.calculationSection} ${styles.frequencySection}`}
          >
            <h2>Vibrational Frequency Analysis</h2>

            {/* Optimization Quality Assessment */}
            <div className={styles.frequencyStatus}>
              <strong>Geometry Optimization Status:</strong>{' '}
              {results.imaginary_frequencies_count === 0 ? (
                <span className={styles.successStatus}>
                  ✅ Successful (no imaginary frequencies)
                </span>
              ) : results.imaginary_frequencies_count === 1 ? (
                <span className={styles.warningStatus}>
                  ⚠️ Possible transition state (1 imaginary frequency)
                </span>
              ) : (
                <span className={styles.errorStatus}>
                  ❌ Poor optimization ({results.imaginary_frequencies_count}{' '}
                  imaginary frequencies)
                </span>
              )}
            </div>

            {/* Vibrational Frequencies */}
            {results.vibrational_frequencies &&
              results.vibrational_frequencies.length > 0 && (
                <div className={styles.frequencyStatus}>
                  <strong>Vibrational Frequencies (cm⁻¹):</strong>
                  <div className={styles.frequencyList}>
                    {results.vibrational_frequencies.map((freq, index) => (
                      <span key={index} className={styles.frequencyItem}>
                        {freq.toFixed(1)}
                      </span>
                    ))}
                  </div>
                  <div className={styles.frequencyCount}>
                    Total: {results.vibrational_frequencies.length} normal modes
                    (≥80 cm⁻¹)
                  </div>
                </div>
              )}

            {/* Thermochemical Properties */}
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
          </section>
        )}

        {/* IR Spectrum Analysis */}
        {results.frequency_analysis_performed && results.vibrational_frequencies && results.vibrational_frequencies.length > 0 && (
          <section className={`${styles.calculationSection} ${styles.irSpectrumSection}`}>
            <h2>IR Spectrum Analysis</h2>
            <div className={styles.sectionDescription}>
              Theoretical infrared spectrum generated from vibrational frequency calculations 
              with scale factor corrections and Lorentzian broadening for realistic peak shapes.
            </div>
            <IRSpectrumViewer
              calculationId={activeCalculation.id}
              onError={(error) => setError(error)}
            />
          </section>
        )}

        {/* Orbital Information */}
        <section
          className={`${styles.calculationSection} ${styles.orbitalSection}`}
        >
          <h2>Molecular Orbitals</h2>
          <div className={styles.orbitalInfoGrid}>
            <div>
              <strong>HOMO Index:</strong> <code>{results.homo_index}</code>
            </div>
            <div>
              <strong>LUMO Index:</strong> <code>{results.lumo_index}</code>
            </div>
            <div>
              <strong>Occupied Orbitals:</strong>{' '}
              <code>{results.num_occupied_orbitals}</code>
            </div>
            <div>
              <strong>Virtual Orbitals:</strong>{' '}
              <code>{results.num_virtual_orbitals}</code>
            </div>
          </div>
        </section>

        {/* Molecular Orbital Energy Diagram */}
        <section
          className={`${styles.calculationSection} ${styles.energyDiagramSection}`}
        >
          <h2>Molecular Orbital Energy Level Diagram</h2>
          <div className={styles.sectionDescription}>
            Energy levels of molecular orbitals are illustrated. Click on
            orbitals to view details in 3D visualization.
          </div>
          <MolecularOrbitalEnergyDiagram
            key={`energy-${activeCalculation.id}`}
            calculationId={activeCalculation.id}
            selectedOrbitalIndex={selectedOrbitalIndex}
            onOrbitalSelect={setSelectedOrbitalIndex}
            onError={error => setError(error)}
          />
        </section>

        {/* Molecular Orbital Visualization */}
        <section
          className={`${styles.calculationSection} ${styles.orbitalViewerSection}`}
        >
          <h2>Molecular Orbital Visualization</h2>
          <div className={styles.sectionDescription}>
            3D visualization of molecular orbitals obtained from quantum
            chemistry calculations. Select orbitals to view their shapes and
            distributions.
          </div>
          <MolecularOrbitalViewer
            key={activeCalculation.id}
            calculationId={activeCalculation.id}
            onError={error => setError(error)}
          />
        </section>

        {/* TDDFT Results Section */}
        {parameters.calculation_method === 'TDDFT' &&
          results.excitation_energies && (
            <>
              {/* Excited States Summary */}
              <section
                className={`${styles.calculationSection} ${styles.excitedStatesSection}`}
              >
                <h2>Excited States Summary</h2>
                <div className={styles.excitedStatesGrid}>
                  <div>
                    <strong>Number of States:</strong>{' '}
                    <code>{results.excitation_energies.length}</code>
                  </div>
                  <div>
                    <strong>TDDFT Method:</strong>{' '}
                    <code>{(parameters as any).tddft_method || 'TDDFT'}</code>
                  </div>
                  <div>
                    <strong>Lowest Excitation:</strong>{' '}
                    <code>
                      {results.excitation_energies[0]?.toFixed(4) || 'N/A'} eV
                    </code>
                  </div>
                  <div>
                    <strong>UV-Vis Range:</strong>{' '}
                    <code>
                      {results.excitation_wavelengths?.[0]?.toFixed(0)} nm
                    </code>
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
                        <th>Transition Type (estimation)</th>
                      </tr>
                    </thead>
                    <tbody>
                      {results.excitation_energies.map(
                        (energy: number, index: number) => {
                          const wavelength =
                            results.excitation_wavelengths?.[index];
                          const oscStrength =
                            results.oscillator_strengths?.[index];
                          const transition = results.major_transitions?.[index];

                          return (
                            <tr key={index}>
                              <td>S{index + 1}</td>
                              <td
                                className={`${styles.rightAlign} ${styles.monoFont}`}
                              >
                                {energy.toFixed(4)}
                              </td>
                              <td
                                className={`${styles.rightAlign} ${styles.monoFont}`}
                              >
                                {wavelength ? wavelength.toFixed(1) : 'N/A'}
                              </td>
                              <td
                                className={`${styles.rightAlign} ${styles.monoFont}`}
                              >
                                {oscStrength !== undefined
                                  ? oscStrength.toFixed(6)
                                  : 'N/A'}
                              </td>
                              <td>
                                {transition?.dominant_transition || 'Unknown'}
                              </td>
                            </tr>
                          );
                        }
                      )}
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
                        const intensity =
                          results.oscillator_strengths?.[index] || 0;
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
                    <text
                      x="780"
                      y="295"
                      fontSize="12"
                      fill="#666"
                      textAnchor="end"
                    >
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
                  UV-Vis absorption spectrum showing calculated transitions.
                  Colors represent approximate wavelength regions.
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
                    NTO analysis provides a more intuitive description of
                    electronic excitations by decomposing the transition density
                    matrix into dominant hole-particle orbital pairs.
                  </div>
                  {results.nto_analysis.map(
                    (stateData: any, stateIndex: number) => (
                      <div key={stateIndex} className={styles.ntoState}>
                        <h3 className={styles.ntoStateTitle}>
                          Excited State S{stateData.state} (
                          {stateData.energy?.toFixed(4)} eV)
                        </h3>
                        {stateData.nto_pairs &&
                        stateData.nto_pairs.length > 0 ? (
                          <div className={styles.tableContainer}>
                            <table
                              className={`${styles.dataTable} ${styles.ntoTable}`}
                            >
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
                                        <span
                                          className={styles.transitionArrow}
                                        >
                                          →
                                        </span>
                                        <span
                                          className={styles.particleOrbital}
                                        >
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
                                          <span
                                            className={styles.contributionValue}
                                          >
                                            {pair.contribution?.toFixed(1) ||
                                              'N/A'}
                                            %
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
                          <div
                            style={{
                              padding: '20px',
                              textAlign: 'center',
                              color: '#666',
                              fontStyle: 'italic',
                            }}
                          >
                            No significant NTO pairs found for this excited
                            state.
                          </div>
                        )}
                        <div
                          style={{
                            marginTop: '10px',
                            fontSize: '12px',
                            color: '#666',
                          }}
                        >
                          Total NTO pairs analyzed:{' '}
                          {stateData.total_nto_pairs || 0}
                        </div>
                      </div>
                    )
                  )}
                  <div
                    style={{
                      marginTop: '20px',
                      padding: '15px',
                      backgroundColor: '#e8f5e8',
                      borderRadius: '4px',
                      border: '1px solid #4caf50',
                    }}
                  >
                    <h4 style={{ margin: '0 0 10px 0', color: '#2e7d32' }}>
                      💡 How to Read NTO Analysis
                    </h4>
                    <ul
                      style={{
                        margin: '0',
                        paddingLeft: '20px',
                        fontSize: '14px',
                        lineHeight: '1.5',
                      }}
                    >
                      <li>
                        <strong>Hole軌道（赤色）</strong>: Orbitals from which
                        electrons are excited (mainly HOMO-type)
                      </li>
                      <li>
                        <strong>Particle軌道（青色）</strong>: Orbitals to which
                        electrons are excited (mainly LUMO-type)
                      </li>
                      <li>
                        <strong>Weight</strong>:
                        その軌道ペアの寄与を表す重み（特異値）
                      </li>
                      <li>
                        <strong>Contribution</strong>:
                        全遷移に対するそのペアの寄与率（%）
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
          )}

        {/* CCSD Results Section */}
        {(parameters.calculation_method === 'CCSD' ||
          parameters.calculation_method === 'CCSD_T') && (
          <section
            className={`${styles.calculationSection} ${styles.ccsdSection}`}
          >
            <h2>CCSD Results</h2>
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
            {(results as any).frozen_core && (
              <div
                className={styles.sectionDescription}
                style={{ marginTop: '15px' }}
              >
                ℹ️ Frozen core approximation was used in this calculation
              </div>
            )}
          </section>
        )}

        {/* Checkpoint File Information */}
        <section
          className={`${styles.calculationSection} ${styles.checkpointSection}`}
        >
          <h2>Checkpoint File Information</h2>
          <div>
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
              <div
                style={{
                  padding: '15px',
                  backgroundColor: '#e8f5e8',
                  borderRadius: '4px',
                  border: '1px solid #4caf50',
                }}
              >
                <h4 style={{ margin: '0 0 10px 0', color: '#2e7d32' }}>
                  📁 File Access Methods
                </h4>
                <div style={{ fontSize: '14px', lineHeight: '1.4' }}>
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
                  <p
                    style={{
                      marginTop: '10px',
                      fontStyle: 'italic',
                      color: '#666',
                    }}
                  >
                    ※ This directory contains molecular orbital data and wave
                    function information
                  </p>
                </div>
              </div>
            )}
          </div>
        </section>

        {/* Mulliken Charges */}
        {results.mulliken_charges && results.mulliken_charges.length > 0 && (
          <section
            style={{
              marginBottom: '30px',
              padding: '20px',
              backgroundColor: '#f0fff8',
              borderRadius: '8px',
            }}
          >
            <h2>Mulliken Charge Analysis</h2>
            <div
              style={{
                marginBottom: '15px',
                fontSize: '14px',
                color: '#666',
              }}
            >
              Partial charges of each atom by Mulliken population analysis.
              Positive values indicate electron deficiency (positive charge),
              negative values indicate electron excess (negative charge).
            </div>
            <div style={{ overflowX: 'auto' }}>
              <table
                style={{
                  width: '100%',
                  borderCollapse: 'collapse',
                  backgroundColor: 'white',
                  borderRadius: '4px',
                  overflow: 'hidden',
                }}
              >
                <thead>
                  <tr style={{ backgroundColor: '#e8f5e8' }}>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      Atom Number
                    </th>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      Element
                    </th>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      Mulliken Charge (e)
                    </th>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      Charge Character
                    </th>
                  </tr>
                </thead>
                <tbody>
                  {results.mulliken_charges.map(
                    (chargeData: any, index: number) => {
                      const isPositive = chargeData.charge > 0;
                      const chargeColor = isPositive ? '#d32f2f' : '#1976d2';
                      const chargeBgColor = isPositive ? '#ffebee' : '#e3f2fd';

                      return (
                        <tr
                          key={index}
                          style={{
                            backgroundColor:
                              index % 2 === 0 ? '#fafafa' : 'white',
                          }}
                        >
                          <td
                            style={{
                              padding: '10px',
                              border: '1px solid #ddd',
                              textAlign: 'center',
                            }}
                          >
                            {chargeData.atom_index + 1}
                          </td>
                          <td
                            style={{
                              padding: '10px',
                              border: '1px solid #ddd',
                              textAlign: 'center',
                              fontWeight: 'bold',
                            }}
                          >
                            {chargeData.element}
                          </td>
                          <td
                            style={{
                              padding: '10px',
                              border: '1px solid #ddd',
                              textAlign: 'center',
                              fontFamily: 'monospace',
                              backgroundColor: chargeBgColor,
                              color: chargeColor,
                              fontWeight: 'bold',
                            }}
                          >
                            {chargeData.charge > 0 ? '+' : ''}
                            {chargeData.charge.toFixed(4)}
                          </td>
                          <td
                            style={{
                              padding: '10px',
                              border: '1px solid #ddd',
                              textAlign: 'center',
                              color: chargeColor,
                              fontWeight: 'bold',
                            }}
                          >
                            {isPositive ? 'Positive (δ+)' : 'Negative (δ−)'}
                          </td>
                        </tr>
                      );
                    }
                  )}
                </tbody>
              </table>
            </div>
            <div
              style={{
                marginTop: '15px',
                padding: '10px',
                backgroundColor: '#e8f5e8',
                borderRadius: '4px',
                fontSize: '14px',
              }}
            >
              <strong>Total Charge:</strong>{' '}
              <code>
                {results.mulliken_charges
                  .reduce((sum: number, charge: any) => sum + charge.charge, 0)
                  .toFixed(4)}{' '}
                e
              </code>{' '}
              (Molecular Charge: <code>{results.charge || 0}</code> e)
            </div>
          </section>
        )}

        {/* Molecular Structure */}
        <section
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#fff5f5',
            borderRadius: '8px',
          }}
        >
          {parameters.calculation_method === 'HF' ||
          parameters.calculation_method === 'MP2' ||
          parameters.calculation_method === 'CCSD' ||
          parameters.calculation_method === 'CCSD_T' ? (
            <>
              <h2>HF-Optimized Geometry</h2>
              <div
                style={{
                  fontSize: '14px',
                  color: '#666',
                  marginBottom: '15px',
                }}
              >
                ℹ️ Geometry optimized using Hartree-Fock method
              </div>
            </>
          ) : (
            <>
              <h2>DFT-Optimized Geometry</h2>
              <div
                style={{
                  fontSize: '14px',
                  color: '#666',
                  marginBottom: '15px',
                }}
              >
                ℹ️ Geometry optimized using DFT method
              </div>
            </>
          )}
          <div>
            <strong>Number of Atoms:</strong> {results.atom_count}
          </div>
          <div style={{ marginTop: '15px' }}>
            <strong>XYZ Coordinates:</strong>
            <pre
              style={{
                backgroundColor: '#f8f8f8',
                padding: '15px',
                borderRadius: '4px',
                overflow: 'auto',
                fontSize: '14px',
                fontFamily: 'monospace',
                border: '1px solid #ddd',
                marginTop: '10px',
              }}
            >
              {results.optimized_geometry}
            </pre>
          </div>
        </section>

        {/* Calculation Parameters */}
        <section
          className={`${styles.calculationSection} ${styles.parametersSection}`}
        >
          <h2>Calculation Parameters</h2>
          <div className={styles.parametersGrid}>
            <div>
              <strong>Max SCF Cycles:</strong> {results.max_cycle}
            </div>
            <div>
              <strong>CPU Cores:</strong> {parameters.cpu_cores || 'Default'}
            </div>
            <div>
              <strong>Memory:</strong>{' '}
              {parameters.memory_mb ? `${parameters.memory_mb} MB` : 'Default'}
            </div>
            <div>
              <strong>Solvent Method:</strong> {parameters.solvent_method}
            </div>
            {parameters.solvent !== '-' && (
              <div>
                <strong>Solvent:</strong> {parameters.solvent}
              </div>
            )}
          </div>
        </section>
      </div>
    </div>
  );
};
