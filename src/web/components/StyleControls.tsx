import React, { useState, useEffect } from "react";
import { StyleSpec } from "../../types/3dmol";

export interface StyleControlsProps {
  onStyleChange: (style: StyleSpec) => void;
  className?: string;
}

export type VisualizationStyle = 
  | "stick"
  | "sphere" 
  | "ball-and-stick"
  | "line"
  | "cartoon";

export interface StyleOption {
  id: VisualizationStyle;
  label: string;
  description: string;
}

const styleOptions: StyleOption[] = [
  {
    id: "stick",
    label: "Stick",
    description: "Show bonds as sticks"
  },
  {
    id: "sphere",
    label: "Space-filling",
    description: "Show atoms as spheres"
  },
  {
    id: "ball-and-stick",
    label: "Ball & Stick",
    description: "Combination of spheres and sticks"
  },
  {
    id: "line",
    label: "Wireframe",
    description: "Show bonds as thin lines"
  },
  {
    id: "cartoon",
    label: "Cartoon",
    description: "Simplified representation for biomolecules"
  }
];

export const StyleControls: React.FC<StyleControlsProps> = ({ onStyleChange, className = "" }) => {
  const [selectedStyle, setSelectedStyle] = useState<VisualizationStyle>("ball-and-stick");
  const [atomRadius, setAtomRadius] = useState(0.3);
  const [bondRadius, setBondRadius] = useState(0.15);

  const generateStyleSpec = (style: VisualizationStyle): StyleSpec => {
    switch (style) {
      case "stick":
        return {
          stick: {
            radius: bondRadius,
            colorscheme: "element"
          }
        };
      
      case "sphere":
        return {
          sphere: {
            radius: atomRadius,
            colorscheme: "element"
          }
        };
      
      case "ball-and-stick":
        return {
          stick: {
            radius: bondRadius,
            colorscheme: "element"
          },
          sphere: {
            radius: atomRadius,
            colorscheme: "element"
          }
        };
      
      case "line":
        return {
          line: {
            linewidth: 2
          }
        };
      
      case "cartoon":
        return {
          cartoon: {
            colorscheme: "element"
          }
        };
      
      default:
        return {
          stick: { radius: 0.2 },
          sphere: { radius: 0.3 }
        };
    }
  };

  // Update 3D viewer as a side effect after state is updated
  useEffect(() => {
    const styleSpec = generateStyleSpec(selectedStyle);
    onStyleChange(styleSpec);
  }, [selectedStyle, atomRadius, bondRadius, onStyleChange]); // dependency array

  return (
    <div className={`style-controls ${className}`}>
      <div className="style-controls-header">
        <h3 className="section-title">Visualization Style</h3>
      </div>

      {/* Style Selection */}
      <div className="style-options">
        {styleOptions.map((option) => (
          <label 
            key={option.id}
            className={`style-option ${selectedStyle === option.id ? 'selected' : ''}`}
          >
            <input
              type="radio"
              name="visualization-style"
              value={option.id}
              checked={selectedStyle === option.id}
              onChange={() => setSelectedStyle(option.id)}
              className="style-radio"
            />
            <div className="option-content">
              <div className="option-label">{option.label}</div>
              <div className="option-description">{option.description}</div>
            </div>
          </label>
        ))}
      </div>

      {/* Size Controls */}
      {(selectedStyle === "sphere" || selectedStyle === "ball-and-stick") && (
        <div className="size-control-section">
          <div className="slider-control">
            <label className="slider-label">
              Atom Size: {atomRadius.toFixed(2)}
            </label>
            <input
              type="range"
              min="0.1"
              max="1.0"
              step="0.05"
              value={atomRadius}
              onChange={(e) => setAtomRadius(parseFloat(e.target.value))}
              className="size-slider"
            />
          </div>
        </div>
      )}

      {(selectedStyle === "stick" || selectedStyle === "ball-and-stick") && (
        <div className="size-control-section">
          <div className="slider-control">
            <label className="slider-label">
              Bond Size: {bondRadius.toFixed(2)}
            </label>
            <input
              type="range"
              min="0.05"
              max="0.5"
              step="0.05"
              value={bondRadius}
              onChange={(e) => setBondRadius(parseFloat(e.target.value))}
              className="size-slider"
            />
          </div>
        </div>
      )}

      {/* Quick Actions */}
      <div className="quick-actions">
        <h4 className="quick-actions-title">Quick Actions</h4>
        <div className="quick-action-buttons">
          <button
            onClick={() => {
              setSelectedStyle("ball-and-stick");
              setAtomRadius(0.3);
              setBondRadius(0.15);
            }}
            className="quick-action-btn"
          >
            Default
          </button>
          <button
            onClick={() => {
              setAtomRadius(0.6);
              setBondRadius(0.4);
            }}
            className="quick-action-btn"
          >
            Large
          </button>
          <button
            onClick={() => {
              setAtomRadius(0.2);
              setBondRadius(0.1);
            }}
            className="quick-action-btn"
          >
            Small
          </button>
        </div>
      </div>
    </div>
  );
};