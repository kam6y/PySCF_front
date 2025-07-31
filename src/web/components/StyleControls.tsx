import React, { useState } from "react";
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
    description: "Simplified representation"
  }
];

export const StyleControls: React.FC<StyleControlsProps> = ({ onStyleChange, className = "" }) => {
  const [selectedStyle, setSelectedStyle] = useState<VisualizationStyle>("ball-and-stick");
  const [atomRadius, setAtomRadius] = useState(0.3);
  const [bondRadius, setBondRadius] = useState(0.2);

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
            radius: bondRadius * 0.8,
            colorscheme: "element"
          },
          sphere: {
            radius: atomRadius * 0.8,
            colorscheme: "element"
          }
        };
      
      case "line":
        return {
          line: {
            linewidth: 2,
            color: "black"
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
          stick: { radius: bondRadius },
          sphere: { radius: atomRadius }
        };
    }
  };

  const handleStyleChange = (style: VisualizationStyle) => {
    setSelectedStyle(style);
    const styleSpec = generateStyleSpec(style);
    onStyleChange(styleSpec);
  };

  const handleRadiusChange = () => {
    const styleSpec = generateStyleSpec(selectedStyle);
    onStyleChange(styleSpec);
  };

  return (
    <div className={`style-controls ${className}`}>
      <div className="style-controls-header" style={{ marginBottom: "15px" }}>
        <h3 style={{ margin: "0 0 10px 0", fontSize: "16px", fontWeight: "bold" }}>
          Visualization Style
        </h3>
      </div>

      {/* Style Selection */}
      <div className="style-options" style={{ marginBottom: "15px" }}>
        {styleOptions.map((option) => (
          <label 
            key={option.id}
            style={{ 
              display: "flex", 
              alignItems: "center", 
              marginBottom: "8px",
              cursor: "pointer",
              padding: "8px",
              backgroundColor: selectedStyle === option.id ? "#e3f2fd" : "transparent",
              borderRadius: "4px",
              border: selectedStyle === option.id ? "1px solid #2196f3" : "1px solid transparent"
            }}
          >
            <input
              type="radio"
              name="visualization-style"
              value={option.id}
              checked={selectedStyle === option.id}
              onChange={() => handleStyleChange(option.id)}
              style={{ marginRight: "8px" }}
            />
            <div>
              <div style={{ fontWeight: "bold", fontSize: "14px" }}>
                {option.label}
              </div>
              <div style={{ fontSize: "12px", color: "#666" }}>
                {option.description}
              </div>
            </div>
          </label>
        ))}
      </div>

      {/* Size Controls */}
      {(selectedStyle === "sphere" || selectedStyle === "ball-and-stick") && (
        <div className="size-controls" style={{ marginBottom: "15px" }}>
          <div style={{ marginBottom: "10px" }}>
            <label style={{ display: "block", fontSize: "14px", fontWeight: "bold", marginBottom: "5px" }}>
              Atom Size: {atomRadius.toFixed(2)}
            </label>
            <input
              type="range"
              min="0.1"
              max="1.0"
              step="0.1"
              value={atomRadius}
              onChange={(e) => {
                setAtomRadius(parseFloat(e.target.value));
                handleRadiusChange();
              }}
              style={{ width: "100%" }}
            />
          </div>
        </div>
      )}

      {(selectedStyle === "stick" || selectedStyle === "ball-and-stick") && (
        <div className="size-controls" style={{ marginBottom: "15px" }}>
          <div style={{ marginBottom: "10px" }}>
            <label style={{ display: "block", fontSize: "14px", fontWeight: "bold", marginBottom: "5px" }}>
              Bond Size: {bondRadius.toFixed(2)}
            </label>
            <input
              type="range"
              min="0.05"
              max="0.5"
              step="0.05"
              value={bondRadius}
              onChange={(e) => {
                setBondRadius(parseFloat(e.target.value));
                handleRadiusChange();
              }}
              style={{ width: "100%" }}
            />
          </div>
        </div>
      )}

      {/* Quick Actions */}
      <div className="quick-actions">
        <h4 style={{ margin: "0 0 8px 0", fontSize: "14px", fontWeight: "bold" }}>
          Quick Actions
        </h4>
        <div style={{ display: "flex", gap: "8px", flexWrap: "wrap" }}>
          <button
            onClick={() => handleStyleChange("ball-and-stick")}
            style={{
              padding: "6px 12px",
              fontSize: "12px",
              backgroundColor: "#f8f9fa",
              border: "1px solid #ccc",
              borderRadius: "4px",
              cursor: "pointer"
            }}
          >
            Default
          </button>
          <button
            onClick={() => {
              setAtomRadius(0.6);
              setBondRadius(0.3);
              handleRadiusChange();
            }}
            style={{
              padding: "6px 12px",
              fontSize: "12px",
              backgroundColor: "#f8f9fa",
              border: "1px solid #ccc",
              borderRadius: "4px",
              cursor: "pointer"
            }}
          >
            Large
          </button>
          <button
            onClick={() => {
              setAtomRadius(0.2);
              setBondRadius(0.1);
              handleRadiusChange();
            }}
            style={{
              padding: "6px 12px",
              fontSize: "12px",
              backgroundColor: "#f8f9fa",
              border: "1px solid #ccc",
              borderRadius: "4px",
              cursor: "pointer"
            }}
          >
            Small
          </button>
        </div>
      </div>
    </div>
  );
};