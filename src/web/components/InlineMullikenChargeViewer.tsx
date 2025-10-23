// src/web/components/InlineMullikenChargeViewer.tsx

import React, { useEffect, useRef, useState, useCallback } from 'react';
import * as $3Dmol from '3dmol';
import { GLViewer } from '../../types/3dmol';
import { useGetCalculationDetails } from '../hooks/useCalculationQueries';
import styles from './InlineMullikenChargeViewer.module.css';

interface InlineMullikenChargeViewerProps {
  calculation_id: string;
  height?: number;
  onError?: (error: string) => void;
}

interface MullikenChargeData {
  atom_index?: number;
  element?: string;
  charge?: number;
}

export const InlineMullikenChargeViewer: React.FC<
  InlineMullikenChargeViewerProps
> = ({ calculation_id, height = 400, onError }) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<GLViewer | null>(null);
  const resizeObserverRef = useRef<ResizeObserver | null>(null);

  const [chargeRange, setChargeRange] = useState(0.5);
  const [isViewerReady, setIsViewerReady] = useState(false);

  // Fetch calculation details to get XYZ and Mulliken charges
  const {
    data: calculationData,
    isLoading: isLoadingCalculation,
    error: calculationError,
  } = useGetCalculationDetails(calculation_id);

  // Extract XYZ and Mulliken charges from calculation data
  const xyzData =
    calculationData?.calculation?.results?.optimized_geometry ||
    calculationData?.calculation?.parameters?.xyz;
  const mullikenCharges: MullikenChargeData[] =
    calculationData?.calculation?.results?.mulliken_charges || [];

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

  // Initialize 3Dmol viewer
  useEffect(() => {
    if (!containerRef.current) return;

    try {
      const viewer = $3Dmol.createViewer(containerRef.current, {
        backgroundColor: 'white',
      });
      viewerRef.current = viewer;
      setIsViewerReady(true);

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
      onError?.('Failed to initialize molecular viewer');
    }

    return () => {
      if (resizeObserverRef.current) {
        resizeObserverRef.current.disconnect();
      }
      if (viewerRef.current) {
        viewerRef.current.clear();
      }
    };
  }, [calculateChargeRange, onError]);

  // Render molecule with Mulliken charges
  const renderMoleculeWithCharges = useCallback(() => {
    if (!viewerRef.current || !xyzData || !mullikenCharges.length) return;

    const viewer = viewerRef.current;
    viewer.clear();

    try {
      // Add molecular structure
      const model = viewer.addModel(xyzData, 'xyz');

      // Map Mulliken charges to atoms
      const atoms = model.atoms;
      if (atoms && atoms.length > 0) {
        atoms.forEach((atom, index) => {
          const chargeData = mullikenCharges.find(
            c => (c.atom_index ?? index) === index
          );
          if (chargeData && chargeData.charge !== undefined) {
            // Store charge as a property on the atom
            if (!atom.properties) {
              atom.properties = {};
            }
            atom.properties.charge = chargeData.charge;
          }
        });
      }

      // Add ball-and-stick style
      viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { radius: 0.4 } });

      // Add Van der Waals surface colored by Mulliken charges
      viewer.addSurface(
        'VDW' as any,
        {
          opacity: 0.7,
          colorscheme: {
            prop: 'charge',
            gradient: 'rwb',
            min: -chargeRange,
            max: chargeRange,
          },
        },
        {},
        {}
      );

      // Add labels for charges
      mullikenCharges.forEach((chargeData, index) => {
        if (chargeData.charge !== undefined) {
          const charge = chargeData.charge;
          const chargeStr = charge.toFixed(3);
          const color =
            charge > 0 ? '#dc2626' : charge < 0 ? '#2563eb' : '#6b7280';

          viewer.addLabel(chargeStr, {
            position: { serial: index + 1 },
            backgroundColor: 'rgba(255, 255, 255, 0.8)',
            fontColor: color,
            fontSize: 10,
            showBackground: true,
            backgroundOpacity: 0.8,
          });
        }
      });

      viewer.zoomTo();
      viewer.render();
    } catch (err) {
      console.error('Failed to render molecule with charges:', err);
      onError?.(
        `Failed to render molecule: ${err instanceof Error ? err.message : String(err)}`
      );
    }
  }, [xyzData, mullikenCharges, chargeRange, onError]);

  // Update visualization when data or charge range changes
  useEffect(() => {
    if (isViewerReady && xyzData && mullikenCharges.length > 0) {
      renderMoleculeWithCharges();
    }
  }, [isViewerReady, xyzData, mullikenCharges, renderMoleculeWithCharges]);

  // Handle errors
  useEffect(() => {
    if (calculationError) {
      const errorMsg =
        calculationError instanceof Error
          ? calculationError.message
          : 'Failed to load calculation data';
      onError?.(errorMsg);
    }
  }, [calculationError, onError]);

  if (isLoadingCalculation) {
    return (
      <div className={styles.container}>
        <div className={styles.loadingContainer}>
          <div className={styles.spinner}></div>
          <div>Loading Mulliken charge data...</div>
        </div>
      </div>
    );
  }

  if (calculationError) {
    return (
      <div className={styles.container}>
        <div className={styles.errorContainer}>
          <div>❌ Failed to load calculation data</div>
          <div className={styles.errorDetail}>
            {calculationError instanceof Error
              ? calculationError.message
              : 'Unknown error'}
          </div>
        </div>
      </div>
    );
  }

  if (!xyzData || !mullikenCharges.length) {
    return (
      <div className={styles.container}>
        <div className={styles.errorContainer}>
          <div>❌ No Mulliken charge data available</div>
          <div className={styles.errorDetail}>
            This calculation does not contain Mulliken population analysis data.
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className={styles.container}>
      {/* Header */}
      <div className={styles.header}>
        <div className={styles.title}>Mulliken Charge Distribution</div>
        <div className={styles.metadata}>
          <span>Atoms: {mullikenCharges.length}</span>
          <span className={styles.separator}>|</span>
          <span>Charge range: ±{chargeRange.toFixed(3)} e</span>
        </div>
      </div>

      {/* 3D Viewer */}
      <div className={styles.viewerContainer}>
        <div
          ref={containerRef}
          style={{
            width: '100%',
            height: `${height}px`,
            backgroundColor: 'white',
          }}
        />
      </div>

      {/* Charge range slider */}
      <div className={styles.controls}>
        <div className={styles.controlRow}>
          <label>Charge Color Range: ±{chargeRange.toFixed(3)} e</label>
          <input
            type="range"
            min="0.01"
            max="1.0"
            step="0.01"
            value={chargeRange}
            onChange={e => setChargeRange(parseFloat(e.target.value))}
          />
        </div>
      </div>
    </div>
  );
};
