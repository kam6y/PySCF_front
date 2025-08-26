import React, { useState, useCallback, useEffect } from 'react';
import {
  parseXYZData,
  getSampleXYZData,
  XYZValidationResult,
} from '../utils/xyzParser';
import styles from './XYZInput.module.css';

export interface XYZInputProps {
  onXYZChange: (xyzData: string, isValid: boolean) => void;
  className?: string;
  value?: string; // Allow external control of the input value
}

export const XYZInput: React.FC<XYZInputProps> = ({
  onXYZChange,
  className = '',
  value,
}) => {
  const [xyzInput, setXyzInput] = useState('');
  const [validationResult, setValidationResult] =
    useState<XYZValidationResult | null>(null);
  const [isValidating, setIsValidating] = useState(false);

  // Update internal state when external value changes
  useEffect(() => {
    if (value !== undefined && value !== xyzInput) {
      setXyzInput(value);
      // Validate the external value without triggering onChange
      if (!value.trim()) {
        setValidationResult(null);
        return;
      }

      setIsValidating(true);
      const timeoutId = setTimeout(() => {
        const result = parseXYZData(value);
        setValidationResult(result);
        setIsValidating(false);
      }, 100); // Shorter delay for external updates

      return () => clearTimeout(timeoutId);
    }
  }, [value]);

  const sampleData = getSampleXYZData();

  const validateXYZ = useCallback(
    (input: string) => {
      if (!input.trim()) {
        setValidationResult(null);
        onXYZChange('', false);
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
    },
    [onXYZChange]
  );

  const handleInputChange = (event: React.ChangeEvent<HTMLTextAreaElement>) => {
    const value = event.target.value;
    setXyzInput(value);
    validateXYZ(value);
  };

  return (
    <div className={`${styles.xyzInputContainer} ${className}`}>
      <textarea
        id="xyz-textarea"
        value={xyzInput}
        onChange={handleInputChange}
        placeholder={`3
Water molecule
O   0.000000   0.000000   0.119262
H   0.000000   0.763239  -0.477047
H   0.000000  -0.763239  -0.477047
`}
        className={styles.xyzTextarea}
      />

      {/* Validation Status */}
      <div className={styles.xyzValidationStatus}>
        {isValidating && (
          <div className={`${styles.validationMessage} ${styles.validating}`}>⏳ Validating...</div>
        )}

        {validationResult && !isValidating && (
          <>
            {validationResult.isValid ? (
              <div className={`${styles.validationMessage} ${styles.valid}`}>
                ✅ Valid XYZ data - {validationResult.data?.numAtoms} atoms
                detected
              </div>
            ) : (
              <div className={`${styles.validationMessage} ${styles.invalid}`}>
                ❌ {validationResult.error}
              </div>
            )}
          </>
        )}
      </div>

      {/* Format Help Link */}
      <div className={styles.xyzFormatHelp}>
        <a
          href="#"
          className={styles.formatHelpLink}
          onClick={e => {
            e.preventDefault();
            const details =
              e.currentTarget.parentElement?.querySelector('details');
            if (details) {
              details.open = !details.open;
            }
          }}
        >
          ▶ XYZ Format Help
        </a>

        <details className={styles.formatHelpDetails}>
          <summary className={styles.formatHelpSummary}></summary>
          <div className={styles.formatHelpContent}>
            <p>
              <strong>XYZ format structure:</strong>
            </p>
            <ol>
              <li>First line: Number of atoms (integer)</li>
              <li>Second line: Comment or description (optional)</li>
              <li>
                Following lines: Element symbol followed by X, Y, Z coordinates
              </li>
            </ol>
            <p>
              <strong>Example:</strong>
            </p>
            <pre className={styles.formatExample}>
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
