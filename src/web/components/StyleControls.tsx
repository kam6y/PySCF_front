import React, { useState, useEffect } from 'react';
import { StyleSpec } from '../../types/3dmol';
import { VAN_DER_WAALS_RADII } from '../data/atomicRadii';
import styles from './StyleControls.module.css';

export interface StyleControlsProps {
  onStyleChange: (style: StyleSpec) => void;
  className?: string;
  showAxes: boolean;
  onShowAxesChange: (show: boolean) => void;
  showCoordinates: boolean;
  onShowCoordinatesChange: (show: boolean) => void;
  useAtomicRadii?: boolean;
  onUseAtomicRadiiChange?: (use: boolean) => void;
}

export type VisualizationStyle = 'stick' | 'sphere' | 'ball-and-stick' | 'line';

export interface StyleOption {
  id: VisualizationStyle;
  label: string;
  description: string;
}

const styleOptions: StyleOption[] = [
  {
    id: 'stick',
    label: 'Stick',
    description: 'Show bonds as sticks',
  },
  {
    id: 'sphere',
    label: 'Space-filling',
    description: 'Show atoms as spheres',
  },
  {
    id: 'ball-and-stick',
    label: 'Ball & Stick',
    description: 'Combination of spheres and sticks',
  },
  {
    id: 'line',
    label: 'Wireframe',
    description: 'Show bonds as thin lines',
  },
];

export const StyleControls: React.FC<StyleControlsProps> = ({
  onStyleChange,
  className = '',
  showAxes,
  onShowAxesChange,
  showCoordinates,
  onShowCoordinatesChange,
  useAtomicRadii = false,
  onUseAtomicRadiiChange,
}) => {
  const [selectedStyle, setSelectedStyle] =
    useState<VisualizationStyle>('ball-and-stick');
  const [atomRadius, setAtomRadius] = useState(0.3);
  const [bondRadius, setBondRadius] = useState(0.15);

  const generateStyleSpec = (style: VisualizationStyle): StyleSpec => {
    switch (style) {
      case 'stick':
        return {
          stick: {
            radius: bondRadius,
            colorscheme: 'default',
          },
        };

      case 'sphere':
        return {
          sphere: {
            radius: atomRadius,
            colorscheme: 'default',
          },
        };

      case 'ball-and-stick':
        return {
          stick: {
            radius: bondRadius,
            colorscheme: 'default',
          },
          sphere: {
            radius: atomRadius,
            colorscheme: 'default',
          },
        };

      case 'line':
        return {
          line: {
            linewidth: 2,
          },
        };

      default:
        return {
          stick: { radius: 0.2 },
          sphere: { radius: 0.3 },
        };
    }
  };

  useEffect(() => {
    const styleSpec = generateStyleSpec(selectedStyle);
    // Add metadata to indicate if atomic radii should be used
    if (useAtomicRadii) {
      (styleSpec as any)._useAtomicRadii = true;
      (styleSpec as any)._baseAtomRadius = atomRadius;
    } else {
      (styleSpec as any)._useAtomicRadii = false;
    }
    onStyleChange(styleSpec);
  }, [selectedStyle, atomRadius, bondRadius, onStyleChange, useAtomicRadii]);

  return (
    <div className={`${styles.styleControls} ${className}`}>
      <div className={styles.styleControlsHeader}>
        <h3 className={styles.sectionTitle}>Visualization Style</h3>
      </div>

      <div className={styles.styleOptions}>
        {styleOptions.map(option => (
          <label
            key={option.id}
            className={`${styles.styleOption} ${selectedStyle === option.id ? styles.selected : ''}`}
          >
            <input
              type="radio"
              name="visualization-style"
              value={option.id}
              checked={selectedStyle === option.id}
              onChange={() => setSelectedStyle(option.id)}
              className={styles.styleRadio}
            />
            <div className={styles.optionContent}>
              <div className={styles.optionLabel}>{option.label}</div>
            </div>
          </label>
        ))}
      </div>

      {(selectedStyle === 'sphere' || selectedStyle === 'ball-and-stick') && (
        <div className={styles.sizeControlSection}>
          <div
            className={`${styles.toggleSwitch} ${styles.toggleSwitchWithMargin}`}
          >
            <span className={styles.toggleLabel}>Use Atomic Radii</span>
            <label className={styles.switch}>
              <input
                type="checkbox"
                checked={useAtomicRadii}
                onChange={e => onUseAtomicRadiiChange?.(e.target.checked)}
              />
              <span className={styles.slider}></span>
            </label>
          </div>
          <div className={styles.sliderControl}>
            <label className={styles.sliderLabel}>
              {useAtomicRadii ? 'Base Size: ' : 'Atom Size: '}
              {atomRadius.toFixed(2)}
            </label>
            <input
              type="range"
              min="0.1"
              max="1.0"
              step="0.05"
              value={atomRadius}
              onChange={e => setAtomRadius(parseFloat(e.target.value))}
              className={styles.sizeSlider}
            />
          </div>
        </div>
      )}

      {(selectedStyle === 'stick' || selectedStyle === 'ball-and-stick') && (
        <div className={styles.sizeControlSection}>
          <div className={styles.sliderControl}>
            <label className={styles.sliderLabel}>
              Bond Size: {bondRadius.toFixed(2)}
            </label>
            <input
              type="range"
              min="0.05"
              max="0.5"
              step="0.05"
              value={bondRadius}
              onChange={e => setBondRadius(parseFloat(e.target.value))}
              className={styles.sizeSlider}
            />
          </div>
        </div>
      )}

      <div className={styles.toggleSwitchSection}>
        <div className={styles.toggleSwitch}>
          <span className={styles.toggleLabel}>Show XYZ Axes</span>
          <label className="switch">
            <input
              type="checkbox"
              checked={showAxes}
              onChange={e => onShowAxesChange(e.target.checked)}
            />
            <span className="slider"></span>
          </label>
        </div>
        <div className={styles.toggleSwitch} style={{ marginTop: '12px' }}>
          <span className={styles.toggleLabel}>Show Atom Coordinates</span>
          <label className="switch">
            <input
              type="checkbox"
              checked={showCoordinates}
              onChange={e => onShowCoordinatesChange(e.target.checked)}
            />
            <span className="slider"></span>
          </label>
        </div>
      </div>
    </div>
  );
};
