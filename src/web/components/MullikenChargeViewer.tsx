import React, { useEffect, useRef, useState, useCallback } from 'react';
import * as $3Dmol from '3dmol';
import { GLViewer } from '../../types/3dmol';
import styles from './MullikenChargeViewer.module.css';

interface MullikenChargeData {
  atom_index?: number;
  element?: string;
  charge?: number;
}

interface MullikenChargeViewerProps {
  xyzData: string;
  mullikenCharges: MullikenChargeData[];
  width?: string;
  height?: string;
}

export const MullikenChargeViewer: React.FC<MullikenChargeViewerProps> = ({
  xyzData,
  mullikenCharges,
  width = '100%',
  height = '500px',
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<GLViewer | null>(null);
  const resizeObserverRef = useRef<ResizeObserver | null>(null);

  const opacity = 0.7; // Fixed opacity
  const [chargeRange, setChargeRange] = useState(0.5);
  const surfaceType = 'VDW'; // Fixed to Van der Waals
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  // Calculate charge range from data
  const calculateChargeRange = useCallback(() => {
    if (!mullikenCharges || mullikenCharges.length === 0) return 0.5;
    const charges = mullikenCharges
      .map(c => c.charge)
      .filter((charge): charge is number => charge !== undefined);
    if (charges.length === 0) return 0.5;
    const maxAbs = Math.max(
      Math.abs(Math.min(...charges)),
      Math.abs(Math.max(...charges))
    );
    return Math.max(maxAbs, 0.1); // Ensure minimum range
  }, [mullikenCharges]);

  // Initialize viewer
  useEffect(() => {
    if (!containerRef.current) return;

    try {
      const viewer = $3Dmol.createViewer(containerRef.current, {
        backgroundColor: 'white',
      });
      viewerRef.current = viewer;

      // Setup resize observer
      resizeObserverRef.current = new ResizeObserver(() => {
        if (viewerRef.current) {
          viewerRef.current.resize();
          viewerRef.current.render();
        }
      });
      resizeObserverRef.current.observe(containerRef.current);

      // Set initial charge range
      setChargeRange(calculateChargeRange());
    } catch (err) {
      console.error('Failed to initialize 3Dmol viewer:', err);
      setError('Failed to initialize molecular viewer');
    }

    return () => {
      if (resizeObserverRef.current) {
        resizeObserverRef.current.disconnect();
      }
      if (viewerRef.current) {
        viewerRef.current.clear();
      }
    };
  }, [calculateChargeRange]);

  // Update visualization when data or settings change
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !xyzData || !mullikenCharges) return;

    setIsLoading(true);
    setError(null);

    try {
      // Clear previous content
      viewer.removeAllModels();
      viewer.removeAllSurfaces();

      // Add model from XYZ data
      const model = viewer.addModel(xyzData, 'xyz');

      // Map Mulliken charges to atoms
      const atoms = model.atoms;
      if (atoms && atoms.length > 0) {
        atoms.forEach((atom, index) => {
          const chargeData = mullikenCharges.find(c => c.atom_index === index);
          if (chargeData && chargeData.charge !== undefined) {
            // Store charge as a property on the atom
            if (!atom.properties) {
              atom.properties = {};
            }
            atom.properties.charge = chargeData.charge;
          }
        });
      }

      // Set molecular structure style (ball and stick)
      viewer.setStyle(
        {},
        {
          stick: { radius: 0.15, colorscheme: 'default' },
          sphere: { radius: 0.3, colorscheme: 'default' },
        }
      );

      // Add surface with charge-based coloring
      // Use string literals as 3Dmol.js expects
      const surfaceTypeValue = surfaceType === 'VDW' ? 'VDW' : 'SAS';

      viewer.addSurface(
        surfaceTypeValue as any,
        {
          opacity: opacity,
          colorscheme: {
            prop: 'charge',
            gradient: 'rwb', // Red (negative) - White (neutral) - Blue (positive)
            min: -chargeRange,
            max: chargeRange,
          },
        },
        {},
        {}
      );

      viewer.zoomTo();
      viewer.render();

      setIsLoading(false);
    } catch (err) {
      console.error('Failed to render charge distribution:', err);
      setError(
        `Failed to render visualization: ${err instanceof Error ? err.message : String(err)}`
      );
      setIsLoading(false);
    }
  }, [xyzData, mullikenCharges, chargeRange]);

  if (error) {
    return (
      <div className={styles.errorContainer}>
        <div className={styles.errorIcon}>❌</div>
        <div className={styles.errorMessage}>{error}</div>
      </div>
    );
  }

  return (
    <div className={styles.chargeViewerContainer}>
      <div className={styles.viewerWrapper}>
        <div
          ref={containerRef}
          className={styles.viewer3D}
          style={{ width, height }}
        />
        {isLoading && (
          <div className={styles.loadingOverlay}>
            <div className={styles.spinner}></div>
            <div>Generating charge distribution...</div>
          </div>
        )}
      </div>

      {/* Vertical Charge Range Slider */}
      <div className={styles.verticalSliderContainer}>
        <div className={styles.sliderLabel}>
          <div>Charge</div>
          <div>Range</div>
        </div>
        <div className={styles.sliderValue}>±{chargeRange.toFixed(2)} e</div>
        <input
          type="range"
          min="0.01"
          max="1.0"
          step="0.01"
          value={chargeRange}
          onChange={e => setChargeRange(parseFloat(e.target.value))}
          className={styles.verticalSlider}
        />
      </div>
    </div>
  );
};
