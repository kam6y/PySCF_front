import styles from '../../pages/CalculationResultsPage.module.css';
import {
  CalculationInstance,
  CalculationParameters,
  CalculationResults,
} from '../../types/api-types';
import { MullikenChargeViewer } from '../MullikenChargeViewer';
import { LazyViewer } from '../LazyViewer';

interface ElectronicPropertiesSectionProps {
  results: CalculationResults;
  parameters: CalculationParameters;
  activeCalculation: Pick<CalculationInstance, 'id'>;
}

export const ElectronicPropertiesSection = ({
  results,
  parameters,
  activeCalculation,
}: ElectronicPropertiesSectionProps) => {
  return (
<section
  className={`${styles.calculationSection} ${styles.electronicPropertiesSection}`}
>
  <h2 className={styles.primaryHeader}>Electronic Properties</h2>

  {/* Dipole Moment - First subsection */}
  {results.dipole_moment_total_debye != null && (
    <div className={styles.propertySubsection}>
      <h3>Dipole Moment</h3>
      <div className={styles.sectionDescription}>
        Electric dipole moment quantifies the separation of positive
        and negative charges in the molecule.
      </div>
      <div className={styles.dipoleMomentContainer}>
        <div className={styles.dipoleComponents}>
          <table className={styles.dipoleTable}>
            <thead>
              <tr>
                <th>Component</th>
                <th className={styles.rightAlign}>Value (Debye)</th>
                <th className={styles.rightAlign}>Value (a.u.)</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>
                  <strong>μx</strong>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_x_debye?.toFixed(4)}
                  </code>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_x_au?.toFixed(4)}
                  </code>
                </td>
              </tr>
              <tr>
                <td>
                  <strong>μy</strong>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_y_debye?.toFixed(4)}
                  </code>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_y_au?.toFixed(4)}
                  </code>
                </td>
              </tr>
              <tr>
                <td>
                  <strong>μz</strong>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_z_debye?.toFixed(4)}
                  </code>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_z_au?.toFixed(4)}
                  </code>
                </td>
              </tr>
              <tr className={styles.totalRow}>
                <td>
                  <strong>|μ| (Total)</strong>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    <strong>
                      {results.dipole_moment_total_debye.toFixed(4)}
                    </strong>
                  </code>
                </td>
                <td className={styles.rightAlign}>
                  <code>
                    {results.dipole_moment_total_au?.toFixed(4)}
                  </code>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>
  )}

  {/* Flex layout for Mulliken Charge List (left) and 3D Visualization (right) */}
  {results.mulliken_charges &&
    results.mulliken_charges.length > 0 && (
      <div className={styles.electronicPropertiesFlexWrapper}>
        {/* Left Column: Mulliken Charge List */}
        <div className={styles.mullikenChargeListColumn}>
          <h3>Mulliken Charge List</h3>
          <div className={styles.sectionDescription}>
            Partial charges of each atom by Mulliken population
            analysis.
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
          <div className={styles.mullikenChargeTableWrapper}>
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
        </div>

        {/* Right Column: 3D Charge Distribution Visualization */}
        {(results.optimized_geometry || parameters.xyz) && (
          <div className={styles.chargeVisualizationColumn}>
            <h3>3D Charge Distribution Visualization</h3>
            <div className={styles.sectionDescription}>
              Interactive 3D visualization of the electrostatic
              potential on the molecular surface.
            </div>
            <LazyViewer>
              <MullikenChargeViewer
                key={activeCalculation.id}
                xyzData={results.optimized_geometry || parameters.xyz}
                mullikenCharges={results.mulliken_charges}
              />
            </LazyViewer>
          </div>
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
  );
};
