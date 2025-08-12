import { useEffect, useState } from "react";
import { CalculationInstance } from "../types/api-types";
// useCalculationPolling を useCalculationSubscription に変更
import { useCalculationSubscription } from "../hooks/useCalculationSubscription";

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
  onCalculationUpdate
}: CalculationResultsPageProps) => {
  const [error, setError] = useState<string | null>(null);

  // --- ここから変更 ---
  // useCalculationPolling を削除し、useCalculationSubscription を呼び出す
  useCalculationSubscription({
    calculationId: activeCalculation?.id || null,
    status: activeCalculation?.status,
    onUpdate: onCalculationUpdate,
    onError: (err) => setError(err)
  });
  // --- ここまで変更 ---
  
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
          <div style={{ color: '#e74c3c', padding: '20px', textAlign: 'center' }}>
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
            📊 No calculation selected. Please select a calculation from the sidebar to view its results.
          </div>
        </div>
      </div>
    );
  }

  // Show message for incomplete calculations
  if (activeCalculation.status !== 'completed' || !activeCalculation.results) {
    const statusMessages = {
      pending: '⏳ This calculation is pending. Please run the calculation first.',
      running: '⚛️ This calculation is currently running. Please wait for completion.',
      error: '❌ This calculation failed. Please check the settings and try again.'
    };
    
    return (
      <div className="page-container">
        <div className="page-content">
          <h1>Calculation Results</h1>
          <div style={{ textAlign: 'center', padding: '20px', color: '#666' }}>
            {statusMessages[activeCalculation.status as keyof typeof statusMessages] || 
             '❓ Calculation results are not available.'}
          </div>
          <div style={{ textAlign: 'center', marginTop: '20px' }}>
            <strong>Calculation:</strong> {activeCalculation.name}<br/>
            <strong>Status:</strong> <span className={`status-badge ${activeCalculation.status}`}>
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
        <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#f8f9fa', borderRadius: '8px' }}>
          <h2>Calculation Summary</h2>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '15px' }}>
            <div>
              <strong>Molecule:</strong> {(parameters as any).name|| 'Unknown'}
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
              <strong>Completed:</strong> {new Date(completedAt).toLocaleString()}
            </div>
            <div>
              <strong>Convergence:</strong> {results.converged ? '✅ Converged' : '❌ Not Converged'}
            </div>
          </div>
        </section>

        {/* Energy Results */}
        <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#e8f6f3', borderRadius: '8px' }}>
          <h2>Energy Results</h2>
          <div style={{ fontSize: '18px' }}>
            <strong>SCF Energy:</strong> <code>{results.scf_energy?.toFixed(8) || 'N/A'} hartree</code>
          </div>
        </section>

        {/* Orbital Information */}
        <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#fef9e7', borderRadius: '8px' }}>
          <h2>Molecular Orbitals</h2>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '15px' }}>
            <div>
              <strong>HOMO Index:</strong> <code>{results.homo_index}</code>
            </div>
            <div>
              <strong>LUMO Index:</strong> <code>{results.lumo_index}</code>
            </div>
            <div>
              <strong>Occupied Orbitals:</strong> <code>{results.num_occupied_orbitals}</code>
            </div>
            <div>
              <strong>Virtual Orbitals:</strong> <code>{results.num_virtual_orbitals}</code>
            </div>
          </div>
        </section>

        {/* TDDFT Results Section */}
        {parameters.calculation_method === 'TDDFT' && results.excitation_energies && (
          <>
            {/* Excited States Summary */}
            <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#e8f0ff', borderRadius: '8px' }}>
              <h2>Excited States Summary</h2>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '15px' }}>
                <div>
                  <strong>Number of States:</strong> <code>{results.excitation_energies.length}</code>
                </div>
                <div>
                  <strong>TDDFT Method:</strong> <code>{(parameters as any).tddft_method || 'TDDFT'}</code>
                </div>
                <div>
                  <strong>Lowest Excitation:</strong> <code>{results.excitation_energies[0]?.toFixed(4) || 'N/A'} eV</code>
                </div>
                <div>
                  <strong>UV-Vis Range:</strong> <code>{results.excitation_wavelengths?.[0]?.toFixed(0)} nm</code>
                </div>
              </div>
            </section>

            {/* Excitation Energies Table */}
            <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#fff8dc', borderRadius: '8px' }}>
              <h2>Excitation Energies and Transitions</h2>
              <div style={{ overflowX: 'auto' }}>
                <table style={{ 
                  width: '100%', 
                  borderCollapse: 'collapse', 
                  backgroundColor: 'white',
                  borderRadius: '4px',
                  overflow: 'hidden'
                }}>
                  <thead>
                    <tr style={{ backgroundColor: '#f0f8ff' }}>
                      <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'left' }}>State</th>
                      <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>Energy (eV)</th>
                      <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>Wavelength (nm)</th>
                      <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>Osc. Strength</th>
                      <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'left' }}>Transition Type</th>
                    </tr>
                  </thead>
                  <tbody>
                    {results.excitation_energies.map((energy: number, index: number) => {
                      const wavelength = results.excitation_wavelengths?.[index];
                      const oscStrength = results.oscillator_strengths?.[index];
                      const transition = results.major_transitions?.[index];
                      
                      return (
                        <tr key={index} style={{ backgroundColor: index % 2 === 0 ? '#fafafa' : 'white' }}>
                          <td style={{ padding: '10px', border: '1px solid #ddd' }}>S{index + 1}</td>
                          <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                            {energy.toFixed(4)}
                          </td>
                          <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                            {wavelength ? wavelength.toFixed(1) : 'N/A'}
                          </td>
                          <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                            {oscStrength !== undefined ? oscStrength.toFixed(6) : 'N/A'}
                          </td>
                          <td style={{ padding: '10px', border: '1px solid #ddd' }}>
                            {transition?.dominant_transition || 'Unknown'}
                          </td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
            </section>

            {/* UV-Vis Spectrum Visualization */}
            <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#f0fff0', borderRadius: '8px' }}>
              <h2>UV-Vis Spectrum (Simulated)</h2>
              <div style={{ 
                height: '300px', 
                border: '1px solid #ddd', 
                borderRadius: '4px', 
                backgroundColor: 'white',
                position: 'relative',
                overflow: 'hidden'
              }}>
                <svg width="100%" height="100%" viewBox="0 0 800 300" style={{ display: 'block' }}>
                  {/* Background Grid */}
                  <defs>
                    <pattern id="grid" width="40" height="30" patternUnits="userSpaceOnUse">
                      <path d="M 40 0 L 0 0 0 30" fill="none" stroke="#f0f0f0" strokeWidth="1"/>
                    </pattern>
                  </defs>
                  <rect width="800" height="300" fill="url(#grid)"/>
                  
                  {/* Spectrum Bars */}
                  {results.excitation_wavelengths?.map((wavelength: number, index: number) => {
                    if (!wavelength || wavelength < 200 || wavelength > 800) return null;
                    
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
                  })}
                  
                  {/* Axis Labels */}
                  <text x="20" y="295" fontSize="12" fill="#666">200nm</text>
                  <text x="400" y="295" fontSize="12" fill="#666">500nm</text>
                  <text x="780" y="295" fontSize="12" fill="#666" textAnchor="end">800nm</text>
                  <text x="10" y="15" fontSize="12" fill="#666" transform="rotate(-90, 10, 15)" textAnchor="end">Intensity</text>
                </svg>
              </div>
              <div style={{ marginTop: '10px', fontSize: '14px', color: '#666', textAlign: 'center' }}>
                UV-Vis absorption spectrum showing calculated transitions. 
                Colors represent approximate wavelength regions.
              </div>
            </section>

            {/* Transition Dipole Moments */}
            {results.transition_dipoles && results.transition_dipoles.length > 0 && (
              <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#fff0f5', borderRadius: '8px' }}>
                <h2>Transition Dipole Moments</h2>
                <div style={{ overflowX: 'auto' }}>
                  <table style={{ 
                    width: '100%', 
                    borderCollapse: 'collapse', 
                    backgroundColor: 'white',
                    borderRadius: '4px',
                    overflow: 'hidden'
                  }}>
                    <thead>
                      <tr style={{ backgroundColor: '#f0f8ff' }}>
                        <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'left' }}>State</th>
                        <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>μx (a.u.)</th>
                        <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>μy (a.u.)</th>
                        <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>μz (a.u.)</th>
                        <th style={{ padding: '12px', border: '1px solid #ddd', textAlign: 'right' }}>|μ| (a.u.)</th>
                      </tr>
                    </thead>
                    <tbody>
                      {results.transition_dipoles.map((dipole: any, index: number) => {
                        const magnitude = Math.sqrt(dipole.x * dipole.x + dipole.y * dipole.y + dipole.z * dipole.z);
                        
                        return (
                          <tr key={index} style={{ backgroundColor: index % 2 === 0 ? '#fafafa' : 'white' }}>
                            <td style={{ padding: '10px', border: '1px solid #ddd' }}>S{index + 1}</td>
                            <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                              {dipole.x.toFixed(6)}
                            </td>
                            <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                              {dipole.y.toFixed(6)}
                            </td>
                            <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                              {dipole.z.toFixed(6)}
                            </td>
                            <td style={{ padding: '10px', border: '1px solid #ddd', textAlign: 'right', fontFamily: 'monospace' }}>
                              <strong>{magnitude.toFixed(6)}</strong>
                            </td>
                          </tr>
                        );
                      })}
                    </tbody>
                  </table>
                </div>
              </section>
            )}
          </>
        )}

        {/* Checkpoint File Information */}
        <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#f0f3ff', borderRadius: '8px' }}>
          <h2>Checkpoint File Information</h2>
          <div>
            <div style={{ marginBottom: '10px' }}>
              <strong>Working Directory:</strong> <code>{results.working_directory}</code>
            </div>
            <div style={{ marginBottom: '10px' }}>
              <strong>Checkpoint File:</strong> <code>{results.checkpoint_file}</code>
            </div>
            <div style={{ marginBottom: '15px' }}>
              <strong>File Status:</strong> {results.checkpoint_exists ? '✅ File exists' : '❌ File not found'}
            </div>
            {results.checkpoint_exists && (
              <div style={{ 
                padding: '15px', 
                backgroundColor: '#e8f5e8', 
                borderRadius: '4px',
                border: '1px solid #4caf50'
              }}>
                <h4 style={{ margin: '0 0 10px 0', color: '#2e7d32' }}>📁 ファイルへのアクセス方法</h4>
                <div style={{ fontSize: '14px', lineHeight: '1.4' }}>
                  <p><strong>Finder:</strong> {results.working_directory} をFinderで開く</p>
                  <p><strong>Terminal:</strong> <code>cd {results.working_directory}</code></p>
                  <p><strong>チェックポイントファイル:</strong> <code>calculation.chk</code></p>
                  <p style={{ marginTop: '10px', fontStyle: 'italic', color: '#666' }}>
                    ※ このディレクトリには分子軌道データや波動関数情報が保存されています
                  </p>
                </div>
              </div>
            )}
          </div>
        </section>

        {/* Molecular Structure */}
        <section style={{ marginBottom: '30px', padding: '20px', backgroundColor: '#fff5f5', borderRadius: '8px' }}>
          <h2>Optimized Geometry</h2>
          <div>
            <strong>Number of Atoms:</strong> {results.atom_count}
          </div>
          <div style={{ marginTop: '15px' }}>
            <strong>XYZ Coordinates:</strong>
            <pre style={{ 
              backgroundColor: '#f8f8f8', 
              padding: '15px', 
              borderRadius: '4px', 
              overflow: 'auto',
              fontSize: '14px',
              fontFamily: 'monospace',
              border: '1px solid #ddd',
              marginTop: '10px'
            }}>
              {results.optimized_geometry}
            </pre>
          </div>
        </section>

        {/* Calculation Parameters */}
        <section style={{ padding: '20px', backgroundColor: '#f5f5f5', borderRadius: '8px' }}>
          <h2>Calculation Parameters</h2>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '15px' }}>
            <div>
              <strong>Max SCF Cycles:</strong> {results.max_cycle}
            </div>
            <div>
              <strong>CPU Cores:</strong> {parameters.cpu_cores || 'Default'}
            </div>
            <div>
              <strong>Memory:</strong> {parameters.memory_mb ? `${parameters.memory_mb} MB` : 'Default'}
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