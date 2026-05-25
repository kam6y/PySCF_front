import React, { useState, useEffect, useCallback } from 'react';
import { StyleSpec } from '../../types/3dmol';
import styles from './StyleControls.module.css';

export interface StyleControlsProps {
  onStyleChange: (style: StyleSpec) => void;
  className?: string;
  showAxes: boolean;
  onShowAxesChange: (show: boolean) => void;
  showCoordinates: boolean;
  onShowCoordinatesChange: (show: boolean) => void;
  showAtomNumbers: boolean;
  onShowAtomNumbersChange: (show: boolean) => void;
  useAtomicRadii?: boolean;
  onUseAtomicRadiiChange?: (use: boolean) => void;
}

export type VisualizationStyle = 'stick' | 'sphere';

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
];

export const StyleControls: React.FC<StyleControlsProps> = ({
  onStyleChange,
  className = '',
  showAxes,
  onShowAxesChange,
  showCoordinates,
  onShowCoordinatesChange,
  showAtomNumbers,
  onShowAtomNumbersChange,
  useAtomicRadii = false,
  onUseAtomicRadiiChange,
}) => {
  const [selectedStyles, setSelectedStyles] = useState<Set<VisualizationStyle>>(
    () => new Set<VisualizationStyle>(['stick', 'sphere'])
  );
  const [atomRadius, setAtomRadius] = useState(0.3);
  const [bondRadius, setBondRadius] = useState(0.15);

  const toggleStyle = useCallback((id: VisualizationStyle) => {
    setSelectedStyles(prev => {
      const next = new Set(prev);
      if (next.has(id)) {
        // Prevent unchecking the last selected style
        if (next.size === 1) return prev;
        next.delete(id);
      } else {
        next.add(id);
      }
      return next;
    });
  }, []);

  const generateStyleSpec = useCallback(
    (selected: Set<VisualizationStyle>): StyleSpec => {
      const spec: StyleSpec = {};
      if (selected.has('stick')) {
        spec.stick = {
          radius: bondRadius,
          colorscheme: 'default',
        };
      }
      if (selected.has('sphere')) {
        spec.sphere = {
          radius: atomRadius,
          colorscheme: 'default',
        };
      }
      return spec;
    },
    [atomRadius, bondRadius]
  );

  useEffect(() => {
    const styleSpec = generateStyleSpec(selectedStyles);
    // Add metadata to indicate if atomic radii should be used
    if (useAtomicRadii) {
      (styleSpec as any)._useAtomicRadii = true;
      (styleSpec as any)._baseAtomRadius = atomRadius;
    } else {
      (styleSpec as any)._useAtomicRadii = false;
    }
    onStyleChange(styleSpec);
  }, [
    selectedStyles,
    generateStyleSpec,
    onStyleChange,
    useAtomicRadii,
    atomRadius,
  ]);

  return (
    <div className={`${styles.styleControls} ${className}`}>
      <div className={styles.styleControlsHeader}>
        <h3 className={styles.sectionTitle}>Visualization Style</h3>
      </div>

      <div className={styles.styleOptions}>
        {styleOptions.map(option => {
          const isSelected = selectedStyles.has(option.id);
          return (
            <label
              key={option.id}
              className={`${styles.styleOption} ${isSelected ? styles.selected : ''}`}
            >
              <input
                type="checkbox"
                value={option.id}
                checked={isSelected}
                onChange={() => toggleStyle(option.id)}
                className={styles.styleRadio}
              />
              <div className={styles.optionContent}>
                <div className={styles.optionLabel}>{option.label}</div>
              </div>
            </label>
          );
        })}
      </div>

      {selectedStyles.has('sphere') && (
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

      {selectedStyles.has('stick') && (
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
        <div className={styles.toggleSwitch} style={{ marginTop: '12px' }}>
          <span className={styles.toggleLabel}>Show Atom Numbers</span>
          <label className="switch">
            <input
              type="checkbox"
              checked={showAtomNumbers}
              onChange={e => onShowAtomNumbersChange(e.target.checked)}
            />
            <span className="slider"></span>
          </label>
        </div>
      </div>
    </div>
  );
};
