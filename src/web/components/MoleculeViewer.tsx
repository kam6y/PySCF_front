import { useRef, useEffect } from 'react';
import * as $3Dmol from '3dmol';
import {
  GLViewer,
  GLModel,
  StyleSpec,
  Label,
  AtomSpec,
  ExtendedStyleSpec,
} from '../../types/3dmol';
import { getAtomicRadius, VAN_DER_WAALS_RADII } from '../data/atomicRadii';
import type { components } from '../types/generated-api';
import styles from './MoleculeViewer.module.css';

type AtomDisplacement = components['schemas']['AtomDisplacement'];

export interface MoleculeViewerProps {
  width?: number | string;
  height?: number | string;
  backgroundColor?: string;
  className?: string;
  vibrationMode?: AtomDisplacement[] | null;
  animationAmplitude?: number;
  xyzData?: string | null;
  currentStyle?: StyleSpec | null;
  showAxes?: boolean;
  showCoordinates?: boolean;
  showAtomNumbers?: boolean;
}

export const MoleculeViewer = ({
  width = '600px',
  height = '600px',
  backgroundColor = 'white',
  className = '',
  vibrationMode = null,
  animationAmplitude = 0.3,
  xyzData = null,
  currentStyle = null,
  showAxes = false,
  showCoordinates = false,
  showAtomNumbers = false,
}: MoleculeViewerProps) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<GLViewer | null>(null);
  const modelRef = useRef<GLModel | null>(null);
  const areAxesVisibleRef = useRef(true);
  const areCoordinatesVisibleRef = useRef(false);
  const areAtomNumbersVisibleRef = useRef(false);
  const animationIntervalRef = useRef<number | null>(null);
  const animationTimeoutRef = useRef<number | null>(null);
  const basePositionsRef = useRef<Array<{ x: number; y: number; z: number }>>(
    []
  );

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

    // --- 原子番号・座標ラベルの描画 ---
    const showCoords = areCoordinatesVisibleRef.current;
    const showNumbers = areAtomNumbersVisibleRef.current;
    if (showCoords || showNumbers) {
      model.atoms.forEach((atom, index) => {
        if (
          atom.x === undefined ||
          atom.y === undefined ||
          atom.z === undefined
        )
          return;

        let text: string;
        if (showNumbers && showCoords) {
          text = `(${index + 1}, ${atom.x.toFixed(4)}, ${atom.y.toFixed(4)}, ${atom.z.toFixed(4)})`;
        } else if (showNumbers) {
          text = (index + 1).toString();
        } else {
          text = `(${atom.x.toFixed(4)}, ${atom.y.toFixed(4)}, ${atom.z.toFixed(4)})`;
        }

        viewer.addLabel(text, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          fontColor: 'black',
          fontSize: 12,
          inFront: true,
          backgroundColor: 'white',
          backgroundOpacity: 0.8,
        });
      });
    }

    viewer.render();
  };

  useEffect(() => {
    if (!containerRef.current) {
      return;
    }

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
        // Clean up animation on unmount
        if (animationIntervalRef.current !== null) {
          clearInterval(animationIntervalRef.current);
        }
      };
    } catch (error) {
      console.error('Failed to initialize 3Dmol viewer:', error);
    }
  }, [backgroundColor]);

  // Load XYZ data when xyzData prop changes
  useEffect(() => {
    if (!xyzData || !viewerRef.current) return;

    try {
      const viewer = viewerRef.current;
      viewer.removeAllModels();
      const model = viewer.addModel(xyzData, 'xyz');
      modelRef.current = model;

      // Apply default ball-and-stick style
      viewer.setStyle(
        {},
        {
          stick: { radius: 0.15, colorscheme: 'default' },
          sphere: { radius: 0.3, colorscheme: 'default' },
        }
      );

      // Store base positions immediately after loading
      if (model && model.atoms) {
        basePositionsRef.current = model.atoms.map(atom => ({
          x: atom.x ?? 0,
          y: atom.y ?? 0,
          z: atom.z ?? 0,
        }));
      }

      updateOverlays(viewer);
      viewer.zoomTo();
      viewer.render();
    } catch (error) {
      console.error('Failed to load XYZ data from props:', error);
    }
  }, [xyzData]);

  // Handle vibration mode changes via props
  useEffect(() => {
    if (vibrationMode && vibrationMode.length > 0) {
      let retryCount = 0;
      const MAX_RETRIES = 50; // 最大50回 (5秒) まで再試行

      // Start animation when vibration mode is provided
      const startAnimation = () => {
        if (!viewerRef.current) {
          retryCount++;
          if (retryCount >= MAX_RETRIES) {
            console.error(
              'MoleculeViewer - viewer failed to initialize after',
              MAX_RETRIES,
              'retries'
            );
            return;
          }
          animationTimeoutRef.current = window.setTimeout(startAnimation, 100);
          return;
        }

        const model = modelRef.current;

        if (!model || !model.atoms || model.atoms.length === 0) {
          retryCount++;
          if (retryCount >= MAX_RETRIES) {
            console.error(
              'MoleculeViewer - model failed to load after',
              MAX_RETRIES,
              'retries'
            );
            return;
          }
          animationTimeoutRef.current = window.setTimeout(startAnimation, 100);
          return;
        }

        // Stop any existing animation
        if (animationIntervalRef.current !== null) {
          clearInterval(animationIntervalRef.current);
          animationIntervalRef.current = null;
        }

        // Store base positions
        basePositionsRef.current = model.atoms.map(atom => ({
          x: atom.x ?? 0,
          y: atom.y ?? 0,
          z: atom.z ?? 0,
        }));

        // Create a map of atom index to displacement
        const displacementMap = new Map<number, AtomDisplacement>();
        vibrationMode.forEach(disp => {
          displacementMap.set(disp.atom_index, disp);
        });

        // Animation parameters
        const fps = 30;
        const period = 2000; // ms for one full cycle
        let startTime = Date.now();

        const animateFrame = () => {
          const viewer = viewerRef.current;
          const currentModel = modelRef.current;
          if (!viewer || !currentModel || !currentModel.atoms) {
            console.warn(
              'AnimateFrame early return - viewer:',
              !!viewer,
              'model:',
              !!currentModel,
              'atoms:',
              currentModel?.atoms?.length
            );
            return;
          }

          const elapsed = Date.now() - startTime;
          const phase = (elapsed % period) / period; // 0 to 1
          const displacement =
            Math.sin(phase * 2 * Math.PI) * animationAmplitude;

          currentModel.atoms.forEach((atom, index) => {
            const basePos = basePositionsRef.current[index];
            const disp = displacementMap.get(index);

            if (basePos && disp) {
              atom.x = basePos.x + (disp.dx ?? 0) * displacement;
              atom.y = basePos.y + (disp.dy ?? 0) * displacement;
              atom.z = basePos.z + (disp.dz ?? 0) * displacement;
            }
          });

          // 座標変更後、スタイルを再設定してジオメトリを再計算
          viewer.setStyle(
            {},
            {
              stick: { radius: 0.15, colorscheme: 'default' },
              sphere: { radius: 0.3, colorscheme: 'default' },
            }
          );
          viewer.render();
        };

        // Start animation loop
        animationIntervalRef.current = window.setInterval(
          animateFrame,
          1000 / fps
        );
      };

      startAnimation();
    } else {
      // Stop animation when vibration mode is cleared
      if (animationIntervalRef.current !== null) {
        clearInterval(animationIntervalRef.current);
        animationIntervalRef.current = null;

        // Restore base positions
        const viewer = viewerRef.current;
        const model = modelRef.current;
        if (
          viewer &&
          model &&
          model.atoms &&
          basePositionsRef.current.length > 0
        ) {
          model.atoms.forEach((atom, index) => {
            const basePos = basePositionsRef.current[index];
            if (basePos) {
              atom.x = basePos.x;
              atom.y = basePos.y;
              atom.z = basePos.z;
            }
          });
          viewer.render();
        }
      }
    }

    return () => {
      // Clean up animation and timeouts when vibrationMode changes
      if (animationIntervalRef.current !== null) {
        clearInterval(animationIntervalRef.current);
        animationIntervalRef.current = null;
      }
      if (animationTimeoutRef.current !== null) {
        clearTimeout(animationTimeoutRef.current);
        animationTimeoutRef.current = null;
      }

      // Restore atom positions to equilibrium when switching modes
      const viewer = viewerRef.current;
      const model = modelRef.current;
      if (
        viewer &&
        model &&
        model.atoms &&
        basePositionsRef.current.length > 0
      ) {
        model.atoms.forEach((atom, index) => {
          const basePos = basePositionsRef.current[index];
          if (basePos) {
            atom.x = basePos.x;
            atom.y = basePos.y;
            atom.z = basePos.z;
          }
        });
        viewer.setStyle(
          {},
          {
            stick: { radius: 0.15, colorscheme: 'default' },
            sphere: { radius: 0.3, colorscheme: 'default' },
          }
        );
        viewer.render();
      }
    };
  }, [vibrationMode, animationAmplitude]);

  // Apply style when currentStyle prop changes
  useEffect(() => {
    if (!currentStyle || !viewerRef.current) return;

    const viewer = viewerRef.current;
    const style = currentStyle as ExtendedStyleSpec;
    const useAtomicRadii = style._useAtomicRadii;
    const baseAtomRadius = style._baseAtomRadius || 0.3;

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
          const elementStyle = { ...style } as ExtendedStyleSpec;
          if (elementStyle.sphere) {
            elementStyle.sphere.radius = atomRadius;
          }

          // Remove metadata
          delete elementStyle._useAtomicRadii;
          delete elementStyle._baseAtomRadius;

          viewer.addStyle(elementSelector, elementStyle);
        }
      } else {
        // Fallback to uniform style if no model loaded
        const fallbackStyle = { ...style } as ExtendedStyleSpec;
        delete fallbackStyle._useAtomicRadii;
        delete fallbackStyle._baseAtomRadius;
        viewer.setStyle({}, fallbackStyle);
      }
    } else {
      // Apply uniform style to all atoms
      const uniformStyle = { ...style } as ExtendedStyleSpec;
      delete uniformStyle._useAtomicRadii;
      delete uniformStyle._baseAtomRadius;
      viewer.setStyle({}, uniformStyle);
    }

    viewer.render();
  }, [currentStyle]);

  // Update axes visibility when showAxes prop changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;

    areAxesVisibleRef.current = showAxes ?? false;
    updateOverlays(viewer);
  }, [showAxes]);

  // Update coordinates visibility when showCoordinates prop changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;

    areCoordinatesVisibleRef.current = showCoordinates ?? false;
    updateOverlays(viewer);
  }, [showCoordinates]);

  // Update atom numbers visibility when showAtomNumbers prop changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;

    areAtomNumbersVisibleRef.current = showAtomNumbers ?? false;
    updateOverlays(viewer);
  }, [showAtomNumbers]);

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
};
