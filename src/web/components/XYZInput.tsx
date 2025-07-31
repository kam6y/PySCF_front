import React, { useState, useCallback } from "react";
import { parseXYZData, getSampleXYZData, XYZValidationResult } from "../utils/xyzParser";

export interface XYZInputProps {
  onXYZChange: (xyzData: string, isValid: boolean) => void;
  className?: string;
}

export const XYZInput: React.FC<XYZInputProps> = ({ onXYZChange, className = "" }) => {
  const [xyzInput, setXyzInput] = useState("");
  const [validationResult, setValidationResult] = useState<XYZValidationResult | null>(null);
  const [isValidating, setIsValidating] = useState(false);

  const sampleData = getSampleXYZData();

  const validateXYZ = useCallback((input: string) => {
    if (!input.trim()) {
      setValidationResult(null);
      onXYZChange("", false);
      return;
    }

    setIsValidating(true);
    
    // Add a small delay to avoid excessive validation during typing
    const timeoutId = setTimeout(() => {
      const result = parseXYZData(input);
      setValidationResult(result);
      onXYZChange(input, result.isValid);
      setIsValidating(false);
    }, 300);

    return () => clearTimeout(timeoutId);
  }, [onXYZChange]);

  const handleInputChange = (event: React.ChangeEvent<HTMLTextAreaElement>) => {
    const value = event.target.value;
    setXyzInput(value);
    validateXYZ(value);
  };

  const handleSampleSelect = (event: React.ChangeEvent<HTMLSelectElement>) => {
    const sampleKey = event.target.value;
    if (sampleKey && sampleData[sampleKey]) {
      const sampleXYZ = sampleData[sampleKey];
      setXyzInput(sampleXYZ);
      validateXYZ(sampleXYZ);
    }
  };

  const handleClear = () => {
    setXyzInput("");
    setValidationResult(null);
    onXYZChange("", false);
  };

  const getValidationStatusStyle = () => {
    if (!validationResult && !isValidating) return {};
    
    if (isValidating) {
      return { borderColor: "#ffa500", backgroundColor: "#fff9e6" };
    }
    
    return validationResult?.isValid 
      ? { borderColor: "#28a745", backgroundColor: "#f8fff8" }
      : { borderColor: "#dc3545", backgroundColor: "#fff8f8" };
  };

  return (
    <div className={`xyz-input-container ${className}`}>
      <div className="xyz-input-header" style={{ marginBottom: "10px", display: "flex", alignItems: "center", gap: "10px" }}>
        <label htmlFor="xyz-textarea" style={{ fontWeight: "bold", fontSize: "14px" }}>
          XYZ Coordinates:
        </label>
        
        <select 
          onChange={handleSampleSelect}
          style={{
            padding: "4px 8px",
            fontSize: "12px",
            border: "1px solid #ccc",
            borderRadius: "4px"
          }}
        >
          <option value="">Load Sample...</option>
          <option value="water">Water (H₂O)</option>
          <option value="methane">Methane (CH₄)</option>
          <option value="benzene">Benzene (C₆H₆)</option>
        </select>

        <button 
          onClick={handleClear}
          style={{
            padding: "4px 8px",
            fontSize: "12px",
            backgroundColor: "#f8f9fa",
            border: "1px solid #ccc",
            borderRadius: "4px",
            cursor: "pointer"
          }}
        >
          Clear
        </button>
      </div>

      <textarea
        id="xyz-textarea"
        value={xyzInput}
        onChange={handleInputChange}
        placeholder={`Enter XYZ coordinates in the format:
<number of atoms>
<comment line (optional)>
<element> <x> <y> <z>
<element> <x> <y> <z>
...

Example:
3
Water molecule
O   0.000000   0.000000   0.119262
H   0.000000   0.763239  -0.477047
H   0.000000  -0.763239  -0.477047`}
        style={{
          width: "100%",
          height: "200px",
          padding: "10px",
          fontSize: "12px",
          fontFamily: "monospace",
          border: "1px solid #ccc",
          borderRadius: "4px",
          resize: "vertical",
          ...getValidationStatusStyle()
        }}
      />

      {/* Validation Status */}
      <div className="validation-status" style={{ marginTop: "8px" }}>
        {isValidating && (
          <div style={{ color: "#ffa500", fontSize: "12px" }}>
            ⏳ Validating...
          </div>
        )}
        
        {validationResult && !isValidating && (
          <div style={{ fontSize: "12px" }}>
            {validationResult.isValid ? (
              <div style={{ color: "#28a745" }}>
                ✅ Valid XYZ data - {validationResult.data?.numAtoms} atoms detected
              </div>
            ) : (
              <div style={{ color: "#dc3545" }}>
                ❌ {validationResult.error}
              </div>
            )}
          </div>
        )}
      </div>

      {/* Format Help */}
      <details style={{ marginTop: "10px", fontSize: "12px", color: "#666" }}>
        <summary style={{ cursor: "pointer", fontWeight: "bold" }}>
          XYZ Format Help
        </summary>
        <div style={{ marginTop: "8px", padding: "8px", backgroundColor: "#f8f9fa", borderRadius: "4px" }}>
          <p><strong>XYZ format structure:</strong></p>
          <ol style={{ paddingLeft: "20px", margin: "8px 0" }}>
            <li>First line: Number of atoms (integer)</li>
            <li>Second line: Comment or description (optional)</li>
            <li>Following lines: Element symbol followed by X, Y, Z coordinates</li>
          </ol>
          <p><strong>Example:</strong></p>
          <pre style={{ backgroundColor: "#fff", padding: "8px", borderRadius: "4px", fontSize: "11px" }}>
{`3
Water molecule
O   0.000000   0.000000   0.119262
H   0.000000   0.763239  -0.477047
H   0.000000  -0.763239  -0.477047`}
          </pre>
        </div>
      </details>
    </div>
  );
};