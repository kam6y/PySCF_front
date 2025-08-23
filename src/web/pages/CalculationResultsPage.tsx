import { useEffect, useState } from 'react';
import { CalculationInstance } from '../types/api-types';
import { MolecularOrbitalViewer } from '../components/MolecularOrbitalViewer';
import { MolecularOrbitalEnergyDiagram } from '../components/MolecularOrbitalEnergyDiagram';

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
  const [selectedOrbitalIndex, setSelectedOrbitalIndex] = useState<number | null>(null);

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
      <div className="page-container">
        <div className="page-content">
          <h1>Calculation Results</h1>
          <div style={{ textAlign: 'center', padding: '20px' }}>
            ⚛️ Loading calculation details...
          </div>
        </div>
      </div>
    );
  }

  // Show error state
  if (error) {
    return (
      <div className="page-container">
        <div className="page-content">
          <h1>Calculation Results</h1>
          <div
            style={{ color: '#e74c3c', padding: '20px', textAlign: 'center' }}
          >
            ❌ {error}
          </div>
        </div>
      </div>
    );
  }

  // Show message when no calculation is selected
  if (!activeCalculation) {
    return (
      <div className="page-container">
        <div className="page-content">
          <h1>Calculation Results</h1>
          <div style={{ textAlign: 'center', padding: '20px', color: '#666' }}>
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
      <div className="page-container">
        <div className="page-content">
          <h1>Calculation Results</h1>
          <div style={{ textAlign: 'center', padding: '20px', color: '#666' }}>
            {statusMessages[
              activeCalculation.status as keyof typeof statusMessages
            ] || '❓ Calculation results are not available.'}
          </div>
          <div style={{ textAlign: 'center', marginTop: '20px' }}>
            <strong>Calculation:</strong> {activeCalculation.name}
            <br />
            <strong>Status:</strong>{' '}
            <span className={`status-badge ${activeCalculation.status}`}>
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
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#f8f9fa',
            borderRadius: '8px',
          }}
        >
          <h2>Calculation Summary</h2>
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
              gap: '15px',
            }}
          >
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
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#e8f6f3',
            borderRadius: '8px',
          }}
        >
          <h2>Energy Results</h2>
          <div style={{ fontSize: '18px' }}>
            <strong>SCF Energy:</strong>{' '}
            <code>{results.scf_energy?.toFixed(8) || 'N/A'} hartree</code>
          </div>
        </section>

        {/* Vibrational Frequency Analysis */}
        {results.frequency_analysis_performed && (
          <section
            style={{
              marginBottom: '30px',
              padding: '20px',
              backgroundColor: '#f0f8ff',
              borderRadius: '8px',
            }}
          >
            <h2>Vibrational Frequency Analysis</h2>
            
            {/* Optimization Quality Assessment */}
            <div style={{ marginBottom: '15px' }}>
              <strong>Geometry Optimization Status:</strong>{' '}
              {results.imaginary_frequencies_count === 0 ? (
                <span style={{ color: '#059669' }}>✅ Successful (no imaginary frequencies)</span>
              ) : results.imaginary_frequencies_count === 1 ? (
                <span style={{ color: '#d97706' }}>⚠️ Possible transition state (1 imaginary frequency)</span>
              ) : (
                <span style={{ color: '#dc2626' }}>❌ Poor optimization ({results.imaginary_frequencies_count} imaginary frequencies)</span>
              )}
            </div>
            
            {/* Vibrational Frequencies */}
            {results.vibrational_frequencies && results.vibrational_frequencies.length > 0 && (
              <div style={{ marginBottom: '15px' }}>
                <strong>Vibrational Frequencies (cm⁻¹):</strong>
                <div 
                  style={{ 
                    marginTop: '8px',
                    padding: '10px',
                    backgroundColor: '#f8fafc',
                    borderRadius: '4px',
                    maxHeight: '150px',
                    overflowY: 'auto',
                    fontFamily: 'monospace',
                    fontSize: '14px'
                  }}
                >
                  {results.vibrational_frequencies.map((freq, index) => (
                    <span key={index} style={{ marginRight: '12px' }}>
                      {freq.toFixed(1)}
                    </span>
                  ))}
                </div>
                <div style={{ fontSize: '12px', color: '#6b7280', marginTop: '5px' }}>
                  Total: {results.vibrational_frequencies.length} normal modes (≥80 cm⁻¹)
                </div>
              </div>
            )}
            
            {/* Thermochemical Properties */}
            <div style={{ 
              display: 'grid', 
              gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))', 
              gap: '15px',
              fontSize: '14px'
            }}>
              {results.zero_point_energy !== undefined && results.zero_point_energy !== null && (
                <div>
                  <strong>Zero-Point Energy:</strong><br />
                  <code>{results.zero_point_energy.toFixed(8)} hartree</code>
                </div>
              )}
              {results.thermal_energy_298K !== undefined && results.thermal_energy_298K !== null && (
                <div>
                  <strong>Thermal Energy (298.15 K):</strong><br />
                  <code>{results.thermal_energy_298K.toFixed(8)} hartree</code>
                </div>
              )}
              {results.entropy_298K !== undefined && results.entropy_298K !== null && (
                <div>
                  <strong>Entropy (298.15 K):</strong><br />
                  <code>{results.entropy_298K.toFixed(8)} hartree/K</code>
                </div>
              )}
              {results.gibbs_free_energy_298K !== undefined && results.gibbs_free_energy_298K !== null && (
                <div>
                  <strong>Gibbs Free Energy (298.15 K):</strong><br />
                  <code>{results.gibbs_free_energy_298K.toFixed(8)} hartree</code>
                </div>
              )}
              {results.heat_capacity_298K !== undefined && results.heat_capacity_298K !== null && (
                <div>
                  <strong>Heat Capacity (298.15 K):</strong><br />
                  <code>{results.heat_capacity_298K.toFixed(8)} hartree/K</code>
                </div>
              )}
            </div>
          </section>
        )}

        {/* Orbital Information */}
        <section
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#fef9e7',
            borderRadius: '8px',
          }}
        >
          <h2>Molecular Orbitals</h2>
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))',
              gap: '15px',
            }}
          >
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
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#fff9e6',
            borderRadius: '8px',
          }}
        >
          <h2>分子軌道エネルギー準位図</h2>
          <div
            style={{ marginBottom: '15px', fontSize: '14px', color: '#666' }}
          >
            分子軌道のエネルギー準位を図示します。軌道をクリックすると3D可視化で詳細を確認できます。
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
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#f0f8ff',
            borderRadius: '8px',
          }}
        >
          <h2>分子軌道可視化</h2>
          <div
            style={{ marginBottom: '15px', fontSize: '14px', color: '#666' }}
          >
            量子化学計算で得られた分子軌道を3D可視化します。軌道を選択して形状や分布を確認できます。
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
                style={{
                  marginBottom: '30px',
                  padding: '20px',
                  backgroundColor: '#e8f0ff',
                  borderRadius: '8px',
                }}
              >
                <h2>Excited States Summary</h2>
                <div
                  style={{
                    display: 'grid',
                    gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))',
                    gap: '15px',
                  }}
                >
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
                style={{
                  marginBottom: '30px',
                  padding: '20px',
                  backgroundColor: '#fff8dc',
                  borderRadius: '8px',
                }}
              >
                <h2>Excitation Energies and Transitions</h2>
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
                      <tr style={{ backgroundColor: '#f0f8ff' }}>
                        <th
                          style={{
                            padding: '12px',
                            border: '1px solid #ddd',
                            textAlign: 'left',
                          }}
                        >
                          State
                        </th>
                        <th
                          style={{
                            padding: '12px',
                            border: '1px solid #ddd',
                            textAlign: 'right',
                          }}
                        >
                          Energy (eV)
                        </th>
                        <th
                          style={{
                            padding: '12px',
                            border: '1px solid #ddd',
                            textAlign: 'right',
                          }}
                        >
                          Wavelength (nm)
                        </th>
                        <th
                          style={{
                            padding: '12px',
                            border: '1px solid #ddd',
                            textAlign: 'right',
                          }}
                        >
                          Osc. Strength
                        </th>
                        <th
                          style={{
                            padding: '12px',
                            border: '1px solid #ddd',
                            textAlign: 'left',
                          }}
                        >
                          Transition Type (estimation)
                        </th>
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
                                }}
                              >
                                S{index + 1}
                              </td>
                              <td
                                style={{
                                  padding: '10px',
                                  border: '1px solid #ddd',
                                  textAlign: 'right',
                                  fontFamily: 'monospace',
                                }}
                              >
                                {energy.toFixed(4)}
                              </td>
                              <td
                                style={{
                                  padding: '10px',
                                  border: '1px solid #ddd',
                                  textAlign: 'right',
                                  fontFamily: 'monospace',
                                }}
                              >
                                {wavelength ? wavelength.toFixed(1) : 'N/A'}
                              </td>
                              <td
                                style={{
                                  padding: '10px',
                                  border: '1px solid #ddd',
                                  textAlign: 'right',
                                  fontFamily: 'monospace',
                                }}
                              >
                                {oscStrength !== undefined
                                  ? oscStrength.toFixed(6)
                                  : 'N/A'}
                              </td>
                              <td
                                style={{
                                  padding: '10px',
                                  border: '1px solid #ddd',
                                }}
                              >
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
                style={{
                  marginBottom: '30px',
                  padding: '20px',
                  backgroundColor: '#f0fff0',
                  borderRadius: '8px',
                }}
              >
                <h2>UV-Vis Spectrum (Simulated)</h2>
                <div
                  style={{
                    height: '300px',
                    border: '1px solid #ddd',
                    borderRadius: '4px',
                    backgroundColor: 'white',
                    position: 'relative',
                    overflow: 'hidden',
                  }}
                >
                  <svg
                    width="100%"
                    height="100%"
                    viewBox="0 0 800 300"
                    style={{ display: 'block' }}
                  >
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
                <div
                  style={{
                    marginTop: '10px',
                    fontSize: '14px',
                    color: '#666',
                    textAlign: 'center',
                  }}
                >
                  UV-Vis absorption spectrum showing calculated transitions.
                  Colors represent approximate wavelength regions.
                </div>
              </section>

              {/* Transition Dipole Moments */}
              {results.transition_dipoles &&
                results.transition_dipoles.length > 0 && (
                  <section
                    style={{
                      marginBottom: '30px',
                      padding: '20px',
                      backgroundColor: '#fff0f5',
                      borderRadius: '8px',
                    }}
                  >
                    <h2>Transition Dipole Moments</h2>
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
                          <tr style={{ backgroundColor: '#f0f8ff' }}>
                            <th
                              style={{
                                padding: '12px',
                                border: '1px solid #ddd',
                                textAlign: 'left',
                              }}
                            >
                              State
                            </th>
                            <th
                              style={{
                                padding: '12px',
                                border: '1px solid #ddd',
                                textAlign: 'right',
                              }}
                            >
                              μx (a.u.)
                            </th>
                            <th
                              style={{
                                padding: '12px',
                                border: '1px solid #ddd',
                                textAlign: 'right',
                              }}
                            >
                              μy (a.u.)
                            </th>
                            <th
                              style={{
                                padding: '12px',
                                border: '1px solid #ddd',
                                textAlign: 'right',
                              }}
                            >
                              μz (a.u.)
                            </th>
                            <th
                              style={{
                                padding: '12px',
                                border: '1px solid #ddd',
                                textAlign: 'right',
                              }}
                            >
                              |μ| (a.u.)
                            </th>
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
                                    }}
                                  >
                                    S{index + 1}
                                  </td>
                                  <td
                                    style={{
                                      padding: '10px',
                                      border: '1px solid #ddd',
                                      textAlign: 'right',
                                      fontFamily: 'monospace',
                                    }}
                                  >
                                    {dipole.x.toFixed(6)}
                                  </td>
                                  <td
                                    style={{
                                      padding: '10px',
                                      border: '1px solid #ddd',
                                      textAlign: 'right',
                                      fontFamily: 'monospace',
                                    }}
                                  >
                                    {dipole.y.toFixed(6)}
                                  </td>
                                  <td
                                    style={{
                                      padding: '10px',
                                      border: '1px solid #ddd',
                                      textAlign: 'right',
                                      fontFamily: 'monospace',
                                    }}
                                  >
                                    {dipole.z.toFixed(6)}
                                  </td>
                                  <td
                                    style={{
                                      padding: '10px',
                                      border: '1px solid #ddd',
                                      textAlign: 'right',
                                      fontFamily: 'monospace',
                                    }}
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
                  style={{
                    marginBottom: '30px',
                    padding: '20px',
                    backgroundColor: '#f0f8ff',
                    borderRadius: '8px',
                  }}
                >
                  <h2>Natural Transition Orbital (NTO) Analysis</h2>
                  <div
                    style={{
                      marginBottom: '15px',
                      fontSize: '14px',
                      color: '#666',
                    }}
                  >
                    NTO analysis provides a more intuitive description of
                    electronic excitations by decomposing the transition density
                    matrix into dominant hole-particle orbital pairs.
                  </div>
                  {results.nto_analysis.map(
                    (stateData: any, stateIndex: number) => (
                      <div key={stateIndex} style={{ marginBottom: '30px' }}>
                        <h3 style={{ color: '#2c5aa0', marginBottom: '15px' }}>
                          Excited State S{stateData.state} (
                          {stateData.energy?.toFixed(4)} eV)
                        </h3>
                        {stateData.nto_pairs &&
                        stateData.nto_pairs.length > 0 ? (
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
                                <tr style={{ backgroundColor: '#e6f3ff' }}>
                                  <th
                                    style={{
                                      padding: '12px',
                                      border: '1px solid #ddd',
                                      textAlign: 'left',
                                    }}
                                  >
                                    NTO Pair
                                  </th>
                                  <th
                                    style={{
                                      padding: '12px',
                                      border: '1px solid #ddd',
                                      textAlign: 'left',
                                    }}
                                  >
                                    Transition
                                  </th>
                                  <th
                                    style={{
                                      padding: '12px',
                                      border: '1px solid #ddd',
                                      textAlign: 'right',
                                    }}
                                  >
                                    Weight
                                  </th>
                                  <th
                                    style={{
                                      padding: '12px',
                                      border: '1px solid #ddd',
                                      textAlign: 'right',
                                    }}
                                  >
                                    Contribution (%)
                                  </th>
                                  <th
                                    style={{
                                      padding: '12px',
                                      border: '1px solid #ddd',
                                      textAlign: 'center',
                                    }}
                                  >
                                    Orbital Indices
                                  </th>
                                </tr>
                              </thead>
                              <tbody>
                                {stateData.nto_pairs.map(
                                  (pair: any, pairIndex: number) => (
                                    <tr
                                      key={pairIndex}
                                      style={{
                                        backgroundColor:
                                          pairIndex % 2 === 0
                                            ? '#fafafa'
                                            : 'white',
                                      }}
                                    >
                                      <td
                                        style={{
                                          padding: '10px',
                                          border: '1px solid #ddd',
                                        }}
                                      >
                                        <strong>#{pairIndex + 1}</strong>
                                      </td>
                                      <td
                                        style={{
                                          padding: '10px',
                                          border: '1px solid #ddd',
                                        }}
                                      >
                                        <span
                                          style={{
                                            color: '#d32f2f',
                                            fontWeight: 'bold',
                                          }}
                                        >
                                          {pair.hole_orbital}
                                        </span>
                                        <span
                                          style={{
                                            margin: '0 8px',
                                            color: '#666',
                                          }}
                                        >
                                          →
                                        </span>
                                        <span
                                          style={{
                                            color: '#1976d2',
                                            fontWeight: 'bold',
                                          }}
                                        >
                                          {pair.particle_orbital}
                                        </span>
                                      </td>
                                      <td
                                        style={{
                                          padding: '10px',
                                          border: '1px solid #ddd',
                                          textAlign: 'right',
                                          fontFamily: 'monospace',
                                        }}
                                      >
                                        {pair.weight?.toFixed(6) || 'N/A'}
                                      </td>
                                      <td
                                        style={{
                                          padding: '10px',
                                          border: '1px solid #ddd',
                                          textAlign: 'right',
                                        }}
                                      >
                                        <div
                                          style={{
                                            display: 'flex',
                                            alignItems: 'center',
                                            justifyContent: 'flex-end',
                                            gap: '8px',
                                          }}
                                        >
                                          <div
                                            style={{
                                              width: `${Math.min(pair.contribution || 0, 100)}%`,
                                              height: '12px',
                                              backgroundColor:
                                                pair.contribution >= 50
                                                  ? '#4caf50'
                                                  : pair.contribution >= 25
                                                    ? '#ff9800'
                                                    : '#f44336',
                                              borderRadius: '6px',
                                              minWidth: '2px',
                                            }}
                                          />
                                          <span
                                            style={{
                                              fontWeight: 'bold',
                                              minWidth: '50px',
                                            }}
                                          >
                                            {pair.contribution?.toFixed(1) ||
                                              'N/A'}
                                            %
                                          </span>
                                        </div>
                                      </td>
                                      <td
                                        style={{
                                          padding: '10px',
                                          border: '1px solid #ddd',
                                          textAlign: 'center',
                                          fontSize: '12px',
                                          color: '#666',
                                        }}
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
                      💡 NTO解析の読み方
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
                        <strong>Hole軌道（赤色）</strong>:
                        励起により電子が抜ける軌道（主にHOMO系）
                      </li>
                      <li>
                        <strong>Particle軌道（青色）</strong>:
                        励起により電子が移る軌道（主にLUMO系）
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
                        寄与率が高いペアほどその励起状態の主要な電子遷移を表しています
                      </li>
                    </ul>
                  </div>
                </section>
              )}
            </>
          )}

        {/* CCSD Results Section */}
        {(parameters.calculation_method === 'CCSD' || parameters.calculation_method === 'CCSD_T') && (
          <section
            style={{
              marginBottom: '30px',
              padding: '20px',
              backgroundColor: '#fff0e8',
              borderRadius: '8px',
            }}
          >
            <h2>CCSD Results</h2>
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))',
                gap: '15px',
              }}
            >
              <div>
                <strong>HF Energy:</strong>{' '}
                <code>{((results as any).hf_energy || results.scf_energy)?.toFixed(6)} Hartree</code>
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
              {parameters.calculation_method === 'CCSD_T' && (results as any).ccsd_t_correction && (
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
                      {(results as any).ccsd_t_total_energy?.toFixed(6)} Hartree
                    </code>
                  </div>
                </>
              )}
            </div>
            {(results as any).frozen_core && (
              <div style={{ marginTop: '15px', fontSize: '14px', color: '#666' }}>
                ℹ️ Frozen core approximation was used in this calculation
              </div>
            )}
          </section>
        )}

        {/* Checkpoint File Information */}
        <section
          style={{
            marginBottom: '30px',
            padding: '20px',
            backgroundColor: '#f0f3ff',
            borderRadius: '8px',
          }}
        >
          <h2>Checkpoint File Information</h2>
          <div>
            <div style={{ marginBottom: '10px' }}>
              <strong>Working Directory:</strong>{' '}
              <code>{results.working_directory}</code>
            </div>
            <div style={{ marginBottom: '10px' }}>
              <strong>Checkpoint File:</strong>{' '}
              <code>{results.checkpoint_file}</code>
            </div>
            <div style={{ marginBottom: '15px' }}>
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
                  📁 ファイルへのアクセス方法
                </h4>
                <div style={{ fontSize: '14px', lineHeight: '1.4' }}>
                  <p>
                    <strong>Finder:</strong> {results.working_directory}{' '}
                    をFinderで開く
                  </p>
                  <p>
                    <strong>Terminal:</strong>{' '}
                    <code>cd {results.working_directory}</code>
                  </p>
                  <p>
                    <strong>チェックポイントファイル:</strong>{' '}
                    <code>calculation.chk</code>
                  </p>
                  <p
                    style={{
                      marginTop: '10px',
                      fontStyle: 'italic',
                      color: '#666',
                    }}
                  >
                    ※
                    このディレクトリには分子軌道データや波動関数情報が保存されています
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
            <h2>Mulliken電荷解析</h2>
            <div
              style={{
                marginBottom: '15px',
                fontSize: '14px',
                color: '#666',
              }}
            >
              Mulliken人口解析による各原子の部分電荷。正の値は電子不足（正電荷）、負の値は電子過剰（負電荷）を示します。
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
                      原子番号
                    </th>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      元素
                    </th>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      Mulliken電荷 (e)
                    </th>
                    <th
                      style={{
                        padding: '12px',
                        border: '1px solid #ddd',
                        textAlign: 'center',
                      }}
                    >
                      電荷の性質
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
                            {chargeData.charge > 0 ? '+' : ''}{chargeData.charge.toFixed(4)}
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
                            {isPositive ? '陽性（δ+）' : '陰性（δ−）'}
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
              <strong>合計電荷:</strong>{' '}
              <code>
                {results.mulliken_charges
                  .reduce((sum: number, charge: any) => sum + charge.charge, 0)
                  .toFixed(4)} e
              </code>
              {' '}(分子電荷: <code>{results.charge || 0}</code> e)
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
          {(parameters.calculation_method === 'HF' || 
            parameters.calculation_method === 'MP2' || 
            parameters.calculation_method === 'CCSD' || 
            parameters.calculation_method === 'CCSD_T') ? (
            <>
              <h2>HF-Optimized Geometry</h2>
              <div style={{ fontSize: '14px', color: '#666', marginBottom: '15px' }}>
                ℹ️ Geometry optimized using Hartree-Fock method
              </div>
            </>
          ) : (
            <>
              <h2>DFT-Optimized Geometry</h2>
              <div style={{ fontSize: '14px', color: '#666', marginBottom: '15px' }}>
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
          style={{
            padding: '20px',
            backgroundColor: '#f5f5f5',
            borderRadius: '8px',
          }}
        >
          <h2>Calculation Parameters</h2>
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
              gap: '15px',
            }}
          >
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
