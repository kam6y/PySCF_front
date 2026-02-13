import React from 'react';
import { MoleculeViewer } from './MoleculeViewer';
import { StyleControls } from './StyleControls';
import { StyleSpec } from '../../types/3dmol';
import styles from './MoleculeViewerSection.module.css';

interface MoleculeViewerSectionProps {
  hasValidMolecule: boolean;
  xyzData?: string | null;
  currentStyle?: StyleSpec | null;
  onStyleChange: (style: StyleSpec) => void;
  showAxes: boolean;
  onShowAxesChange: (show: boolean) => void;
  showCoordinates: boolean;
  onShowCoordinatesChange: (show: boolean) => void;
  useAtomicRadii: boolean;
  onUseAtomicRadiiChange: (use: boolean) => void;
}

export const MoleculeViewerSection = React.memo<MoleculeViewerSectionProps>(({
  hasValidMolecule,
  xyzData,
  currentStyle,
  onStyleChange,
  showAxes,
  onShowAxesChange,
  showCoordinates,
  onShowCoordinatesChange,
  useAtomicRadii,
  onUseAtomicRadiiChange,
}) => {
  return (
    <div className={styles.mainContent}>
      <div className={styles.leftColumn}>
        <MoleculeViewer
          xyzData={xyzData}
          currentStyle={currentStyle}
          showAxes={showAxes}
          showCoordinates={showCoordinates}
          width={'100%'}
          height={'100%'}
          backgroundColor="white"
        />
        {!hasValidMolecule && (
          <div className={styles.viewerPlaceholder}>
            <div className={styles.placeholderContent}>
              <h3>No Molecule Loaded</h3>
              <p>
                Enter a molecular structure in the right panel to see the 3D
                visualization
              </p>
            </div>
          </div>
        )}
      </div>
      <section className={styles.rightColumn}>
        <StyleControls
          onStyleChange={onStyleChange}
          showAxes={showAxes}
          onShowAxesChange={onShowAxesChange}
          showCoordinates={showCoordinates}
          onShowCoordinatesChange={onShowCoordinatesChange}
          useAtomicRadii={useAtomicRadii}
          onUseAtomicRadiiChange={onUseAtomicRadiiChange}
        />
      </section>
    </div>
  );
});
