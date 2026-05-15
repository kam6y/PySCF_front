import React from 'react';
import styles from '../../pages/CalculationResultsPage.module.css';
import { MolecularOrbitalViewer } from '../MolecularOrbitalViewer';
import { MolecularOrbitalEnergyDiagram } from '../MolecularOrbitalEnergyDiagram';
import { LazyViewer } from '../LazyViewer';

interface MolecularOrbitalsSectionProps {
  calculationId: string;
  selectedOrbitalIndex: number | null;
  onOrbitalSelect: (orbitalIndex: number) => void;
  onError: (error: string) => void;
}

export const MolecularOrbitalsSection =
  React.memo<MolecularOrbitalsSectionProps>(
    ({ calculationId, selectedOrbitalIndex, onOrbitalSelect, onError }) => {
      return (
        <section
          className={`${styles.calculationSection} ${styles.molecularOrbitalsSection}`}
        >
          <h2 className={styles.primaryHeader}>Molecular Orbitals</h2>

          {/* Flex wrapper for horizontal layout */}
          <div className={styles.orbitalsFlexWrapper}>
            {/* Molecular Orbital Energy Diagram */}
            <div className={styles.orbitalEnergyDiagram}>
              <h3>Energy Level Diagram</h3>
              <div className={styles.sectionDescription}>
                Energy levels of molecular orbitals are illustrated.
              </div>
              <LazyViewer>
                <MolecularOrbitalEnergyDiagram
                  key={`energy-${calculationId}`}
                  calculationId={calculationId}
                  selectedOrbitalIndex={selectedOrbitalIndex}
                  onOrbitalSelect={onOrbitalSelect}
                  onError={onError}
                />
              </LazyViewer>
            </div>

            {/* Molecular Orbital 3D Visualization */}
            <div className={styles.orbitalVisualization}>
              <h3>3D Orbital Visualization</h3>
              <div className={styles.sectionDescription}>
                Interactive 3D visualization of molecular orbitals.
              </div>
              <LazyViewer>
                <MolecularOrbitalViewer
                  key={calculationId}
                  calculationId={calculationId}
                  selectedOrbitalIndex={selectedOrbitalIndex}
                  onOrbitalSelect={onOrbitalSelect}
                  onError={onError}
                />
              </LazyViewer>
            </div>
          </div>
        </section>
      );
    }
  );
