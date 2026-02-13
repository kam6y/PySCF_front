import styles from '../../pages/CalculationResultsPage.module.css';
import {
  CalculationInstance,
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';

interface CalculationSettingsSummaryProps {
  activeCalculation: CalculationInstance;
  processedData: {
    results: CalculationResults;
    parameters: CalculationParameters;
    shouldShowTDDFTSection: boolean;
    shouldShowCASSection: boolean;
  };
}

export const CalculationSettingsSummary = ({
  activeCalculation,
  processedData,
}: CalculationSettingsSummaryProps) => {
  const { results, parameters } = processedData;
  const completedAt = activeCalculation.updatedAt;
  const computeDeviceLabel = results.gpu_enabled ? 'GPU' : 'CPU';

  return (
<section className={styles.calculationSettings}>
  <h2 className={styles.primaryHeader}>Calculation Settings</h2>

  <div className={styles.settingsGrid}>
    {/* Basic Information Category */}
    <div className={styles.categorySection}>
      <h3 className={styles.categoryTitle}>Basic Information</h3>
      <div className={styles.infoGrid}>
        <div className={styles.infoRow}>
          <span className={styles.label}>Calculation Name:</span>
          <span className={styles.value}>{activeCalculation.name}</span>
        </div>
        <div className={styles.infoRow}>
          <span className={styles.label}>Calculation Method:</span>
          <span className={styles.value}>
            {parameters.calculation_method}
          </span>
        </div>
        <div className={styles.infoRow}>
          <span className={styles.label}>Compute Device:</span>
          <span className={styles.value}>{computeDeviceLabel}</span>
        </div>
        <div className={styles.infoRow}>
          <span className={styles.label}>Completed At:</span>
          <span className={styles.value}>
            {new Date(completedAt).toLocaleString()}
          </span>
        </div>
        <div className={styles.infoRow}>
          <span className={styles.label}>Convergence:</span>
          <span className={styles.value}>
            {results.converged ? 'Converged ✓' : 'Not Converged'}
          </span>
        </div>
      </div>
    </div>

    {/* Molecular Configuration Category */}
    <div className={styles.categorySection}>
      <h3 className={styles.categoryTitle}>Molecular Configuration</h3>
      <div className={styles.infoGrid}>
        <div className={styles.infoRow}>
          <span className={styles.label}>Charge:</span>
          <span className={styles.value}>{results.charge}</span>
        </div>
        <div className={styles.infoRow}>
          <span className={styles.label}>Spin (2S):</span>
          <span className={styles.value}>{results.spin}</span>
        </div>
        {results.total_electrons !== undefined &&
          results.total_electrons !== null && (
            <div className={styles.infoRow}>
              <span className={styles.label}>Total Electrons:</span>
              <span className={styles.value}>
                {results.total_electrons}
              </span>
            </div>
          )}
      </div>
    </div>

    {/* Basis Set Configuration Category */}
    <div className={styles.categorySection}>
      <h3 className={styles.categoryTitle}>Basis Set Configuration</h3>
      <div className={styles.infoGrid}>
        <div className={styles.infoRow}>
          <span className={styles.label}>Basis Set:</span>
          <span className={styles.value}>{results.basis}</span>
        </div>
        {results.xc_functional && (
          <div className={styles.infoRow}>
            <span className={styles.label}>XC Functional:</span>
            <span className={styles.value}>
              {results.xc_functional}
            </span>
          </div>
        )}
        {results.num_basis_functions !== undefined &&
          results.num_basis_functions !== null && (
            <div className={styles.infoRow}>
              <span className={styles.label}>
                Number of Basis Functions:
              </span>
              <span className={styles.value}>
                {results.num_basis_functions}
              </span>
            </div>
          )}
        {results.num_primitive_gaussians !== undefined &&
          results.num_primitive_gaussians !== null && (
            <div className={styles.infoRow}>
              <span className={styles.label}>
                Number of Primitive Gaussians:
              </span>
              <span className={styles.value}>
                {results.num_primitive_gaussians}
              </span>
            </div>
          )}
      </div>
    </div>

    {/* Solvation Effects Category (Conditional) */}
    {parameters.solvent && parameters.solvent !== '-' && (
      <div className={styles.categorySection}>
        <h3 className={styles.categoryTitle}>Solvation Effects</h3>
        <div className={styles.infoGrid}>
          <div className={styles.infoRow}>
            <span className={styles.label}>Solvation Method:</span>
            <span className={styles.value}>
              {parameters.solvent_method || 'none'}
            </span>
          </div>
          <div className={styles.infoRow}>
            <span className={styles.label}>Solvent:</span>
            <span className={styles.value}>{parameters.solvent}</span>
          </div>
        </div>
      </div>
    )}

    {/* TDDFT Configuration Category (Conditional) */}
    {processedData.shouldShowTDDFTSection && (
      <div className={styles.categorySection}>
        <h3 className={styles.categoryTitle}>TDDFT Configuration</h3>
        <div className={styles.infoGrid}>
          {(parameters as any).tddft_nstates !== undefined &&
            (parameters as any).tddft_nstates !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Number of Excited States:
                </span>
                <span className={styles.value}>
                  {(parameters as any).tddft_nstates}
                </span>
              </div>
            )}
          {(parameters as any).tddft_method && (
            <div className={styles.infoRow}>
              <span className={styles.label}>TDDFT Method:</span>
              <span className={styles.value}>
                {(parameters as any).tddft_method}
              </span>
            </div>
          )}
        </div>
      </div>
    )}

    {/* CASCI/CASSCF Configuration Category (Conditional) */}
    {processedData.shouldShowCASSection && (
      <div className={styles.categorySection}>
        <h3 className={styles.categoryTitle}>
          {parameters.calculation_method === 'CASSCF'
            ? 'CASSCF'
            : 'CASCI'}{' '}
          Configuration
        </h3>
        <div className={styles.infoGrid}>
          {(parameters as any).ncas !== undefined &&
            (parameters as any).ncas !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Active Space Orbitals (ncas):
                </span>
                <span className={styles.value}>
                  {(parameters as any).ncas}
                </span>
              </div>
            )}
          {(parameters as any).nelecas !== undefined &&
            (parameters as any).nelecas !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Active Space Electrons (nelecas):
                </span>
                <span className={styles.value}>
                  {(parameters as any).nelecas}
                </span>
              </div>
            )}
          {parameters.calculation_method === 'CASSCF' &&
            (parameters as any).max_cycle_macro !== undefined &&
            (parameters as any).max_cycle_macro !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Max Macro Iterations:
                </span>
                <span className={styles.value}>
                  {(parameters as any).max_cycle_macro}
                </span>
              </div>
            )}
          {(parameters as any).max_cycle_micro !== undefined &&
            (parameters as any).max_cycle_micro !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Max CI Micro Iterations:
                </span>
                <span className={styles.value}>
                  {(parameters as any).max_cycle_micro}
                </span>
              </div>
            )}
          {(parameters as any).natorb !== undefined &&
            (parameters as any).natorb !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Natural Orbital Transform:
                </span>
                <span className={styles.value}>
                  {(parameters as any).natorb ? 'Enabled' : 'Disabled'}
                </span>
              </div>
            )}
          {(parameters as any).conv_tol !== undefined &&
            (parameters as any).conv_tol !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Energy Convergence Tolerance:
                </span>
                <span className={styles.value}>
                  {(parameters as any).conv_tol}
                </span>
              </div>
            )}
          {parameters.calculation_method === 'CASSCF' &&
            (parameters as any).conv_tol_grad !== undefined &&
            (parameters as any).conv_tol_grad !== null && (
              <div className={styles.infoRow}>
                <span className={styles.label}>
                  Gradient Convergence Tolerance:
                </span>
                <span className={styles.value}>
                  {(parameters as any).conv_tol_grad}
                </span>
              </div>
            )}
        </div>
      </div>
    )}
  </div>
</section>
  );
};
