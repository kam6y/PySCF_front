import React, {
  useRef,
  useEffect,
  useImperativeHandle,
  forwardRef,
} from 'react';
import * as $3Dmol from '3dmol';
import {
  GLViewer,
  GLModel,
  StyleSpec,
  Label,
  AtomSpec,
} from '../../types/3dmol';
import { getAtomicRadius, VAN_DER_WAALS_RADII } from '../data/atomicRadii';
import styles from './MoleculeViewer.module.css';

export interface MoleculeViewerProps {
  width?: number | string;
  height?: number | string;
  backgroundColor?: string;
  className?: string;
}

export interface MoleculeViewerRef {
  loadXYZ: (xyzData: string) => void;
  setStyle: (style: StyleSpec) => void;
  clearModels: () => void;
  zoomToFit: () => void;
  takeScreenshot: () => string;
  showAxes: (shouldShow: boolean) => void;
  showAtomCoordinates: (shouldShow: boolean) => void;
}

export const MoleculeViewer = forwardRef<
  MoleculeViewerRef,
  MoleculeViewerProps
>(
  (
    {
      width = '600px',
      height = '600px',
      backgroundColor = 'white',
      className = '',
    },
    ref
  ) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const viewerRef = useRef<GLViewer | null>(null);
    const areAxesVisibleRef = useRef(true);
    const areCoordinatesVisibleRef = useRef(false);

    /**
     * 軸やラベルなどのオーバーレイをすべて更新する統一関数
     */
    const updateOverlays = (viewer: GLViewer) => {
      // 既存のシェイプとラベルをすべてクリア
      viewer.removeAllShapes();
      viewer.removeAllLabels();

      const model = viewer.getModel();
      if (!model || !model.atoms || model.atoms.length === 0) {
        viewer.render();
        return; // モデルがない場合は何も描画しない
      }

      // --- 軸の描画 (有効な場合) ---
      if (areAxesVisibleRef.current) {
        let maxDist = 0;
        const atoms = model.atoms;
        let xmin = Infinity,
          xmax = -Infinity,
          ymin = Infinity,
          ymax = -Infinity,
          zmin = Infinity,
          zmax = -Infinity;

        for (const atom of atoms) {
          if (
            atom.x === undefined ||
            atom.y === undefined ||
            atom.z === undefined
          )
            continue;
          xmin = Math.min(xmin, atom.x);
          xmax = Math.max(xmax, atom.x);
          ymin = Math.min(ymin, atom.y);
          ymax = Math.max(ymax, atom.y);
          zmin = Math.min(zmin, atom.z);
          zmax = Math.max(zmax, atom.z);
        }
        const dx = xmax - xmin;
        const dy = ymax - ymin;
        const dz = zmax - zmin;
        maxDist = Math.sqrt(dx * dx + dy * dy + dz * dz);

        const axisLength = maxDist > 0 ? maxDist * 0.35 : 3.0;
        const radius = axisLength * 0.05;
        const labelOffset = radius * 4;

        viewer.addArrow({
          start: { x: 0, y: 0, z: 0 },
          end: { x: axisLength, y: 0, z: 0 },
          radius,
          color: 'red',
        });
        viewer.addLabel('X', {
          position: { x: axisLength + labelOffset, y: 0, z: 0 },
          fontColor: 'red',
          fontSize: 32,
          bold: true,
          backgroundOpacity: 0,
        });

        viewer.addArrow({
          start: { x: 0, y: 0, z: 0 },
          end: { x: 0, y: axisLength, z: 0 },
          radius,
          color: 'green',
        });
        viewer.addLabel('Y', {
          position: { x: 0, y: axisLength + labelOffset, z: 0 },
          fontColor: 'green',
          fontSize: 32,
          bold: true,
          backgroundOpacity: 0,
        });

        viewer.addArrow({
          start: { x: 0, y: 0, z: 0 },
          end: { x: 0, y: 0, z: axisLength },
          radius,
          color: 'blue',
        });
        viewer.addLabel('Z', {
          position: { x: 0, y: 0, z: axisLength + labelOffset },
          fontColor: 'blue',
          fontSize: 32,
          bold: true,
          backgroundOpacity: 0,
        });
      }

      // --- 原子座標の描画 (有効な場合) ---
      if (areCoordinatesVisibleRef.current) {
        model.atoms.forEach(atom => {
          if (
            atom.x === undefined ||
            atom.y === undefined ||
            atom.z === undefined
          )
            return;
          const text = `(${atom.x.toFixed(4)}, ${atom.y.toFixed(4)}, ${atom.z.toFixed(4)})`;
          viewer.addLabel(text, {
            position: { x: atom.x, y: atom.y, z: atom.z },
            fontColor: '#333333',
            fontSize: 16,
            inFront: true,
            backgroundColor: 'white',
            backgroundOpacity: 0.6,
          });
        });
      }

      viewer.render();
    };

    useEffect(() => {
      if (!containerRef.current) return;
      let viewer: GLViewer | null = null;
      try {
        viewer = $3Dmol.createViewer(containerRef.current, { backgroundColor });
        viewerRef.current = viewer;
        const resizeObserver = new ResizeObserver(() => {
          viewer?.resize();
        });
        resizeObserver.observe(containerRef.current);
        return () => {
          resizeObserver.disconnect();
        };
      } catch (error) {
        console.error('Failed to initialize 3Dmol viewer:', error);
      }
    }, [backgroundColor]);

    useImperativeHandle(ref, () => ({
      loadXYZ: (xyzData: string) => {
        const viewer = viewerRef.current;
        if (!viewer) return;
        try {
          viewer.removeAllModels();
          viewer.addModel(xyzData, 'xyz');

          // Apply default ball-and-stick style
          viewer.setStyle(
            {},
            {
              stick: { radius: 0.15, colorscheme: 'default' },
              sphere: { radius: 0.3, colorscheme: 'default' },
            }
          );

          updateOverlays(viewer);
          viewer.zoomTo();
        } catch (error) {
          console.error('Failed to load XYZ data:', error);
        }
      },

      clearModels: () => {
        const viewer = viewerRef.current;
        if (!viewer) return;
        viewer.removeAllModels();
        viewer.removeAllShapes();
        viewer.removeAllLabels();
        viewer.render();
      },

      showAxes: (shouldShow: boolean) => {
        areAxesVisibleRef.current = shouldShow;
        const viewer = viewerRef.current;
        if (viewer) {
          updateOverlays(viewer);
        }
      },

      showAtomCoordinates: (shouldShow: boolean) => {
        areCoordinatesVisibleRef.current = shouldShow;
        const viewer = viewerRef.current;
        if (viewer) {
          updateOverlays(viewer);
        }
      },

      setStyle: (style: StyleSpec) => {
        const viewer = viewerRef.current;
        if (!viewer) return;

        const useAtomicRadii = (style as any)._useAtomicRadii;
        const baseAtomRadius = (style as any)._baseAtomRadius || 0.3;

        if (useAtomicRadii) {
          // Apply atomic radii: set styles for each element separately
          const model = viewer.getModel();
          if (model && model.atoms) {
            // First clear existing styles
            viewer.setStyle({}, {});

            // Group atoms by element
            const elementGroups: Record<string, AtomSpec[]> = {};
            for (const atom of model.atoms) {
              if (atom.elem) {
                if (!elementGroups[atom.elem]) {
                  elementGroups[atom.elem] = [];
                }
                elementGroups[atom.elem].push(atom);
              }
            }

            // Apply styles for each element with appropriate radius
            for (const [element, atoms] of Object.entries(elementGroups)) {
              const atomRadius = getAtomicRadius(
                element,
                VAN_DER_WAALS_RADII['H'],
                baseAtomRadius
              );
              const elementSelector = { elem: element };

              // Create element-specific style
              const elementStyle = { ...style };
              if (elementStyle.sphere) {
                elementStyle.sphere.radius = atomRadius;
              }

              // Remove metadata
              delete (elementStyle as any)._useAtomicRadii;
              delete (elementStyle as any)._baseAtomRadius;

              viewer.addStyle(elementSelector, elementStyle);
            }
          } else {
            // Fallback to uniform style if no model loaded
            const fallbackStyle = { ...style };
            delete (fallbackStyle as any)._useAtomicRadii;
            delete (fallbackStyle as any)._baseAtomRadius;
            viewer.setStyle({}, fallbackStyle);
          }
        } else {
          // Apply uniform style to all atoms
          const uniformStyle = { ...style };
          delete (uniformStyle as any)._useAtomicRadii;
          delete (uniformStyle as any)._baseAtomRadius;
          viewer.setStyle({}, uniformStyle);
        }

        viewer.render();
      },

      zoomToFit: () => {
        viewerRef.current?.zoomTo();
        viewerRef.current?.render();
      },

      takeScreenshot: () => {
        if (viewerRef.current && containerRef.current) {
          const w = containerRef.current.clientWidth;
          const h = containerRef.current.clientHeight;
          return viewerRef.current.screenshot(w, h);
        }
        return '';
      },
    }));

    return (
      <div
        className={`${styles.moleculeViewer} ${className}`}
        style={{
          width,
          height,
        }}
      >
        <div ref={containerRef} className={styles.moleculeViewerContainer} />
      </div>
    );
  }
);
