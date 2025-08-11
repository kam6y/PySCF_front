import { useEffect, useState } from "react";
import { CalculationInstance } from "../types/api-types";
import { useCalculationPolling } from "../hooks/useCalculationPolling";

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

  const { startPolling, stopPolling } = useCalculationPolling({
    calculationId: activeCalculation?.id || null,
    onUpdate: onCalculationUpdate
  });

  useEffect(() => {
    setError(detailsError);
  }, [detailsError]);

  useEffect(() => {
    // If activeCalculation is running, start polling for its status
    if (activeCalculation && activeCalculation.status === 'running') {
      startPolling();
    } else {
      stopPolling();
    }
  }, [activeCalculation, startPolling, stopPolling]);

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