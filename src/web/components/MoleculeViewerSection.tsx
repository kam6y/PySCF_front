import { RefObject } from 'react';
import { MoleculeViewer, MoleculeViewerRef } from './MoleculeViewer';
import { StyleControls } from './StyleControls';
import { StyleSpec } from '../../types/3dmol';

interface MoleculeViewerSectionProps {
  moleculeViewerRef: RefObject<MoleculeViewerRef | null>;
  hasValidMolecule: boolean;
  onStyleChange: (style: StyleSpec) => void;
  showAxes: boolean;
  onShowAxesChange: (show: boolean) => void;
  showCoordinates: boolean;
  onShowCoordinatesChange: (show: boolean) => void;
  useAtomicRadii: boolean;
  onUseAtomicRadiiChange: (use: boolean) => void;
}

export const MoleculeViewerSection = ({
  moleculeViewerRef,
  hasValidMolecule,
  onStyleChange,
  showAxes,
  onShowAxesChange,
  showCoordinates,
  onShowCoordinatesChange,
  useAtomicRadii,
  onUseAtomicRadiiChange,
}: MoleculeViewerSectionProps) => {
  return (
    <div className="main-content">
      <div className="left-column">
        <MoleculeViewer
          ref={moleculeViewerRef}
          width={'100%'}
          height={'100%'}
          backgroundColor="white"
        />
        {!hasValidMolecule && (
          <div className="viewer-placeholder">
            <div className="placeholder-content">
              <h3>No Molecule Loaded</h3>
              <p>
                Enter a molecular structure in the right panel to see the 3D
                visualization
              </p>
            </div>
          </div>
        )}
      </div>
      <div className="right-column">
        <section className="visualization-section">
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
    </div>
  );
};
