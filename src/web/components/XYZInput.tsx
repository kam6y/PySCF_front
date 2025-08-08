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
      <textarea
        id="xyz-textarea"
        value={xyzInput}
        onChange={handleInputChange}
        placeholder={`C -1.151228 0.497765 -1.393083
N -2.046769 -0.227518 -1.197056
C -3.273621 0.358425 -0.274278
C -3.383621 1.488582 -0.274278
C -2.238126 0.144473 -1.132090
C -2.317317 -0.752840 -0.250109
...`}
        className="xyz-textarea"
      />

      {/* Validation Status */}
      <div className="xyz-validation-status">
        {isValidating && (
          <div className="validation-message validating">
            ⏳ Validating...
          </div>
        )}
        
        {validationResult && !isValidating && (
          <>
            {validationResult.isValid ? (
              <div className="validation-message valid">
                ✅ Valid XYZ data - {validationResult.data?.numAtoms} atoms detected
              </div>
            ) : (
              <div className="validation-message invalid">
                ❌ {validationResult.error}
              </div>
            )}
          </>
        )}
      </div>

      {/* Format Help Link */}
      <div className="xyz-format-help">
        <a href="#" className="format-help-link" onClick={(e) => {
          e.preventDefault();
          const details = e.currentTarget.parentElement?.querySelector('details');
          if (details) {
            details.open = !details.open;
          }
        }}>
          ▶ XYZ Format Help
        </a>
        
        <details className="format-help-details">
          <summary className="format-help-summary"></summary>
          <div className="format-help-content">
            <p><strong>XYZ format structure:</strong></p>
            <ol>
              <li>First line: Number of atoms (integer)</li>
              <li>Second line: Comment or description (optional)</li>
              <li>Following lines: Element symbol followed by X, Y, Z coordinates</li>
            </ol>
            <p><strong>Example:</strong></p>
            <pre className="format-example">
{`3
Water molecule
O   0.000000   0.000000   0.119262
H   0.000000   0.763239  -0.477047
H   0.000000  -0.763239  -0.477047`}
            </pre>
          </div>
        </details>
      </div>
    </div>
  );
};