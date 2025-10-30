import React, { useState, useEffect, useCallback } from 'react';
import styles from './CalculationResultsPage.module.css';
import { CalculationInstance } from '../types/api-types';
import { MolecularOrbitalViewer } from '../components/MolecularOrbitalViewer';
import { MolecularOrbitalEnergyDiagram } from '../components/MolecularOrbitalEnergyDiagram';
import { IRSpectrumViewer } from '../components/IRSpectrumViewer';
import { CIAnalysisViewer } from '../components/CIAnalysisViewer';
import { MullikenChargeViewer } from '../components/MullikenChargeViewer';
import { LazyViewer } from '../components/LazyViewer';
import { useProcessedCalculationResults } from '../hooks/useProcessedCalculationResults';

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
  const [isMullikenListOpen, setIsMullikenListOpen] = useState(false);

  useEffect(() => {
    setError(detailsError);
  }, [detailsError]);

  // Optimize event handlers with useCallback
  const handleSetError = useCallback((errorMessage: string) => {
    setError(errorMessage);
  }, []);

  const handleOrbitalSelect = useCallback((orbitalIndex: number) => {
    setSelectedOrbitalIndex(orbitalIndex);
  }, []);

  // Process and memoize calculation results data
  const processedData = useProcessedCalculationResults(activeCalculation);

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

  // Early return if processedData is null (already handled by loading/error states above)
  if (!processedData) {
    return null;
  }

  const { results, parameters, primaryEnergyLabel, primaryEnergyValue } =
    processedData;
  const completedAt = activeCalculation.updatedAt;

  return (
    <div className="page-container">
      <div className="page-content">
        <h1>Quantum Chemistry Calculation Results</h1>

        {/* ========================================
            1️⃣ OVERVIEW SECTION - Integrated Summary and Energy
            ======================================== */}
        <section
          className={`${styles.calculationSection} ${styles.overviewSection}`}
        >
          <h2 className={styles.primaryHeader}>Calculation Overview</h2>

          {/* Basic Information */}
          <div className={styles.overviewGrid}>
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
              <strong>Spin (2S):</strong> {results.spin}
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

          {/* Primary Energy Result - Prominently Displayed */}
          <div className={styles.primaryEnergyResult}>
            <div className={styles.energyLabel}>Primary Energy Result</div>
            <div className={styles.energyValue}>
              <strong>{primaryEnergyLabel}:</strong>{' '}
              <code>{primaryEnergyValue}</code>
            </div>
          </div>
        </section>

        {/* ========================================
            2️⃣ MOLECULAR STRUCTURE SECTION
            ======================================== */}
        <section
          className={`${styles.calculationSection} ${styles.structureSection}`}
        >
          {(() => {
            // Check if geometry optimization was performed
            const optimizeGeometry =
              (parameters as any).optimize_geometry !== false;

            if (!optimizeGeometry) {
              return (
                <>
                  <h2 className={styles.primaryHeader}>
                    Initial Molecular Structure
                  </h2>
                  <div className={styles.sectionDescription}>
                    ℹ️ No geometry optimization performed - using initial
                    structure
                  </div>
                </>
              );
            }

            // Determine optimization method based on calculation method
            switch (parameters.calculation_method) {
              case 'MP2':
                return (
                  <>
                    <h2 className={styles.primaryHeader}>
                      MP2-Optimized Molecular Structure
                    </h2>
                    <div className={styles.sectionDescription}>
                      ℹ️ Geometry optimized using MP2 method
                    </div>
                  </>
                );
              case 'HF':
                return (
                  <>
                    <h2 className={styles.primaryHeader}>
                      HF-Optimized Molecular Structure
                    </h2>
                    <div className={styles.sectionDescription}>
                      ℹ️ Geometry optimized using Hartree-Fock method
                    </div>
                  </>
                );
              case 'CCSD':
              case 'CCSD_T':
                return (
                  <>
                    <h2 className={styles.primaryHeader}>
                      Initial Molecular Structure
                    </h2>
                    <div className={styles.sectionDescription}>
                      ℹ️ CCSD calculations use initial geometry (no
                      optimization)
                    </div>
                  </>
                );
              case 'CASCI':
              case 'CASSCF':
                return (
                  <>
                    <h2 className={styles.primaryHeader}>
                      Initial Molecular Structure
                    </h2>
                    <div className={styles.sectionDescription}>
                      ℹ️ CASCI/CASSCF calculations use initial geometry (no
                      optimization)
                    </div>
                  </>
                );
              case 'TDDFT':
                return (
                  <>
                    <h2 className={styles.primaryHeader}>
                      Initial Molecular Structure
                    </h2>
                    <div className={styles.sectionDescription}>
                      ℹ️ TDDFT calculations use initial geometry (no
                      optimization)
                    </div>
                  </>
                );
              default:
                // DFT and other methods
                return (
                  <>
                    <h2 className={styles.primaryHeader}>
                      DFT-Optimized Molecular Structure
                    </h2>
                    <div className={styles.sectionDescription}>
                      ℹ️ Geometry optimized using DFT method
                    </div>
                  </>
                );
            }
          })()}
          <div className={styles.molecularStructureInfo}>
            <strong>Number of Atoms:</strong> {results.atom_count}
          </div>
          <div className={styles.xyzCoordinatesContainer}>
            <strong>XYZ Coordinates:</strong>
            <pre className={styles.xyzCoordinates}>
              {results.optimized_geometry}
            </pre>
          </div>
        </section>

        {/* ========================================
            3️⃣ ELECTRONIC PROPERTIES SECTION - New Unified Section
            ======================================== */}
        {processedData.shouldShowElectronicProperties && (
          <section
            className={`${styles.calculationSection} ${styles.electronicPropertiesSection}`}
          >
            <h2 className={styles.primaryHeader}>Electronic Properties</h2>

            {/* 3D Charge Distribution Visualization - Now first */}
            {results.mulliken_charges &&
              results.mulliken_charges.length > 0 &&
              results.optimized_geometry && (
                <div className={styles.propertySubsection}>
                  <h3>3D Charge Distribution Visualization</h3>
                  <div className={styles.sectionDescription}>
                    Interactive 3D visualization of the electrostatic potential
                    on the molecular surface. The surface color represents
                    charge distribution based on Mulliken population analysis.
                  </div>
                  <LazyViewer>
                    <MullikenChargeViewer
                      key={activeCalculation.id}
                      xyzData={results.optimized_geometry}
                      mullikenCharges={results.mulliken_charges}
                    />
                  </LazyViewer>
                </div>
              )}

            {/* Mulliken Charge List - Now second, with toggle */}
            {results.mulliken_charges &&
              results.mulliken_charges.length > 0 && (
                <div className={styles.propertySubsection}>
                  <div
                    className={styles.toggleHeader}
                    onClick={() => setIsMullikenListOpen(!isMullikenListOpen)}
                  >
                    <h3>Mulliken Charge List</h3>
                    <span
                      className={`${styles.toggleIcon} ${isMullikenListOpen ? styles.rotated : ''}`}
                    >
                      <svg
                        width="20"
                        height="20"
                        viewBox="0 0 16 16"
                        fill="none"
                        xmlns="http://www.w3.org/2000/svg"
                      >
                        <path
                          d="M4 6L8 10L12 6"
                          stroke="currentColor"
                          strokeWidth="1.5"
                          strokeLinecap="round"
                          strokeLinejoin="round"
                        />
                      </svg>
                    </span>
                  </div>
                  {isMullikenListOpen && (
                    <>
                      <div className={styles.sectionDescription}>
                        Partial charges of each atom by Mulliken population
                        analysis. Positive values indicate electron deficiency
                        (positive charge), negative values indicate electron
                        excess (negative charge).
                      </div>
                      <div className={styles.tableContainer}>
                        <table className={styles.mullikenChargeTable}>
                          <thead>
                            <tr>
                              <th>Atom Number</th>
                              <th>Element</th>
                              <th>Mulliken Charge (e)</th>
                              <th>Charge Character</th>
                            </tr>
                          </thead>
                          <tbody>
                            {results.mulliken_charges.map(
                              (chargeData: any, index: number) => {
                                const isPositive = chargeData.charge > 0;

                                return (
                                  <tr key={index}>
                                    <td>{chargeData.atom_index + 1}</td>
                                    <td style={{ fontWeight: 'bold' }}>
                                      {chargeData.element}
                                    </td>
                                    <td
                                      className={`${styles.chargeValueCell} ${
                                        isPositive
                                          ? styles.chargeValueCellPositive
                                          : styles.chargeValueCellNegative
                                      }`}
                                    >
                                      {chargeData.charge > 0 ? '+' : ''}
                                      {chargeData.charge.toFixed(4)}
                                    </td>
                                    <td
                                      className={
                                        isPositive
                                          ? styles.chargeCharacterPositive
                                          : styles.chargeCharacterNegative
                                      }
                                    >
                                      {isPositive
                                        ? 'Positive (δ+)'
                                        : 'Negative (δ−)'}
                                    </td>
                                  </tr>
                                );
                              }
                            )}
                          </tbody>
                        </table>
                      </div>
                      <div className={styles.chargeSummary}>
                        <strong>Total Charge:</strong>{' '}
                        <code>
                          {results.mulliken_charges
                            .reduce(
                              (sum: number, charge: any) => sum + charge.charge,
                              0
                            )
                            .toFixed(4)}{' '}
                          e
                        </code>{' '}
                        (Molecular Charge: <code>{results.charge || 0}</code> e)
                      </div>
                    </>
                  )}
                </div>
              )}

            {/* Mulliken Spin Density Analysis - moved from CASCI/CASSCF section */}
            {(results as any).mulliken_spin_analysis &&
              (results as any).mulliken_spin_analysis.available && (
                <div className={styles.propertySubsection}>
                  <h3>Mulliken Spin Density Analysis</h3>
                  <div className={styles.sectionDescription}>
                    Mulliken atomic spin densities show the distribution of
                    unpaired electron density across atoms in open-shell
                    systems.
                  </div>

                  <div className={styles.spinSummary}>
                    <div className={styles.spinStats}>
                      <div>
                        <strong>Total Spin Density:</strong>{' '}
                        <code>
                          {(
                            results as any
                          ).mulliken_spin_analysis.total_spin_density?.toFixed(
                            3
                          ) || 'N/A'}
                        </code>
                      </div>
                      <div>
                        <strong>Expected Spin:</strong>{' '}
                        <code>
                          {(results as any).mulliken_spin_analysis
                            .expected_spin || 'N/A'}
                        </code>
                      </div>
                      <div>
                        <strong>Total Absolute Spin:</strong>{' '}
                        <code>
                          {(
                            results as any
                          ).mulliken_spin_analysis.total_absolute_spin_density?.toFixed(
                            3
                          ) || 'N/A'}
                        </code>
                      </div>
                    </div>
                  </div>

                  {(results as any).mulliken_spin_analysis
                    .atomic_spin_densities &&
                    (results as any).mulliken_spin_analysis
                      .atomic_spin_densities.length > 0 && (
                      <div className={styles.atomicSpinSection}>
                        <h4>Atomic Spin Densities</h4>
                        <div className={styles.tableContainer}>
                          <table className={styles.dataTable}>
                            <thead>
                              <tr>
                                <th>Atom</th>
                                <th>Element</th>
                                <th className={styles.rightAlign}>
                                  Spin Density
                                </th>
                                <th className={styles.rightAlign}>
                                  |Spin Density|
                                </th>
                              </tr>
                            </thead>
                            <tbody>
                              {(
                                results as any
                              ).mulliken_spin_analysis.atomic_spin_densities.map(
                                (atom: any, index: number) => (
                                  <tr key={index}>
                                    <td>{atom.atom_index + 1}</td>
                                    <td>{atom.element}</td>
                                    <td className={styles.rightAlign}>
                                      <code
                                        style={{
                                          color:
                                            atom.spin_density > 0
                                              ? '#2e7d32'
                                              : atom.spin_density < 0
                                                ? '#d32f2f'
                                                : 'inherit',
                                        }}
                                      >
                                        {atom.spin_density.toFixed(3)}
                                      </code>
                                    </td>
                                    <td className={styles.rightAlign}>
                                      <code>
                                        {atom.abs_spin_density.toFixed(3)}
                                      </code>
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
          </section>
        )}

        {/* ========================================
            4️⃣ METHOD-SPECIFIC ADVANCED RESULTS
            ======================================== */}
        {/* CASCI/CASSCF Results Section */}
        {processedData.shouldShowCASSection && (
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
        )}

        {/* ========================================
            5️⃣ VIBRATIONAL ANALYSIS SECTION - Unified
            ======================================== */}
        {results.frequency_analysis_performed && (
          <section
            className={`${styles.calculationSection} ${styles.vibrationalSection}`}
          >
            <h2 className={styles.primaryHeader}>Vibrational Analysis</h2>

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

            {/* IR Spectrum - Integrated into Vibrational Analysis */}
            {results.vibrational_frequencies &&
              results.vibrational_frequencies.length > 0 && (
                <div className={styles.irSpectrumSubsection}>
                  <h3>Infrared (IR) Spectrum</h3>
                  <div className={styles.sectionDescription}>
                    Theoretical infrared spectrum generated from vibrational
                    frequency calculations with scale factor corrections and
                    Lorentzian broadening for realistic peak shapes.
                  </div>
                  <LazyViewer>
                    <IRSpectrumViewer
                      calculationId={activeCalculation.id}
                      onError={handleSetError}
                    />
                  </LazyViewer>
                </div>
              )}
          </section>
        )}

        {/* ========================================
            6️⃣ MOLECULAR ORBITALS SECTION - Unified
            ======================================== */}
        <section
          className={`${styles.calculationSection} ${styles.molecularOrbitalsSection}`}
        >
          <h2 className={styles.primaryHeader}>Molecular Orbitals</h2>

          {/* Basic Orbital Information */}
          <div className={styles.orbitalBasicInfo}>
            <h3>Orbital Information</h3>
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
          </div>

          {/* Molecular Orbital Energy Diagram */}
          <div className={styles.orbitalEnergyDiagram}>
            <h3>Energy Level Diagram</h3>
            <div className={styles.sectionDescription}>
              Energy levels of molecular orbitals are illustrated. Click on
              orbitals to view details in 3D visualization below.
            </div>
            <LazyViewer>
              <MolecularOrbitalEnergyDiagram
                key={`energy-${activeCalculation.id}`}
                calculationId={activeCalculation.id}
                selectedOrbitalIndex={selectedOrbitalIndex}
                onOrbitalSelect={handleOrbitalSelect}
                onError={handleSetError}
              />
            </LazyViewer>
          </div>

          {/* Molecular Orbital 3D Visualization */}
          <div className={styles.orbitalVisualization}>
            <h3>3D Orbital Visualization</h3>
            <div className={styles.sectionDescription}>
              Interactive 3D visualization of molecular orbitals. Select
              orbitals from the energy diagram above or use the controls to view
              their shapes and spatial distributions.
            </div>
            <LazyViewer>
              <MolecularOrbitalViewer
                key={activeCalculation.id}
                calculationId={activeCalculation.id}
                onError={handleSetError}
              />
            </LazyViewer>
          </div>
        </section>

        {/* TDDFT Results Section */}
        {processedData.shouldShowTDDFTSection &&
          results.excitation_energies && (
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
                          <div className={styles.ntoNoData}>
                            No significant NTO pairs found for this excited
                            state.
                          </div>
                        )}
                        <div className={styles.ntoCount}>
                          Total NTO pairs analyzed:{' '}
                          {stateData.total_nto_pairs || 0}
                        </div>
                      </div>
                    )
                  )}
                  <div className={styles.ntoHelpBox}>
                    <h4 className={styles.ntoHelpTitle}>
                      💡 How to Read NTO Analysis
                    </h4>
                    <ul className={styles.ntoHelpList}>
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
        {processedData.shouldShowCCSDSection && (
          <section
            className={`${styles.calculationSection} ${styles.ccsdSection}`}
          >
            <h2 className={styles.secondaryHeader}>CCSD Advanced Results</h2>
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
              <div className={styles.sectionDescription}>
                ℹ️ Frozen core approximation was used in this calculation
              </div>
            )}
          </section>
        )}

        {/* ========================================
            7️⃣ TECHNICAL DETAILS SECTION - Bottom
            ======================================== */}
        <section
          className={`${styles.calculationSection} ${styles.technicalDetailsSection}`}
        >
          <h2 className={styles.primaryHeader}>Technical Details</h2>

          {/* Calculation Parameters */}
          <div className={styles.technicalSubsection}>
            <h3>Calculation Parameters</h3>
            <div className={styles.parametersGrid}>
              <div>
                <strong>Max SCF Cycles:</strong> {results.max_cycle}
              </div>
              <div>
                <strong>CPU Cores:</strong> {parameters.cpu_cores || 'Default'}
              </div>
              <div>
                <strong>Memory:</strong>{' '}
                {parameters.memory_mb
                  ? `${parameters.memory_mb} MB`
                  : 'Default'}
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
          </div>

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
                    ※ This directory contains molecular orbital data and wave
                    function information
                  </p>
                </div>
              </div>
            )}
          </div>
        </section>
      </div>
    </div>
  );
};
