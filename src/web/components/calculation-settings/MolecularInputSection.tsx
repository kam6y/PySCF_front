import React from 'react';
import styles from '../../pages/CalculationSettingsPage.module.css';
import { XYZInput } from '../XYZInput';
import { INPUT_PLACEHOLDERS } from '../../constants/calculationDefaults';

interface MolecularInputSectionProps {
  inputMethod: string;
  pubchemInput: string;
  isConverting: boolean;
  convertError: string | null;
  calculationStatus: string;
  xyzInputValue: string;
  onInputMethodChange: (method: string) => void;
  onPubchemInputChange: (value: string) => void;
  onXYZConvert: () => void;
  onXYZChange: (xyzData: string, isValid: boolean) => void;
  onConvertErrorClear: () => void;
}

export const MolecularInputSection = React.memo<MolecularInputSectionProps>(
  ({
    inputMethod,
    pubchemInput,
    isConverting,
    convertError,
    calculationStatus,
    xyzInputValue,
    onInputMethodChange,
    onPubchemInputChange,
    onXYZConvert,
    onXYZChange,
    onConvertErrorClear,
  }) => {
    return (
      <section className={styles.molecularInputSection}>
        <h3>Molecular Structure Input</h3>
        <div className={styles.inputMethodSelection}>
          <div className={styles.radioOptions}>
            <label className={styles.radioOption}>
              <input
                type="radio"
                name="inputMethod"
                value="pubchem"
                checked={inputMethod === 'pubchem'}
                onChange={e => onInputMethodChange(e.target.value)}
                disabled={calculationStatus === 'running'}
              />
              <span className={styles.radioText}>Get from PubChem Name/CID</span>
            </label>
            <label className={styles.radioOption}>
              <input
                type="radio"
                name="inputMethod"
                value="smiles"
                checked={inputMethod === 'smiles'}
                onChange={e => onInputMethodChange(e.target.value)}
                disabled={calculationStatus === 'running'}
              />
              <span className={styles.radioText}>Get from SMILES</span>
            </label>
          </div>
        </div>
        <div className={styles.pubchemInputSection}>
          <div className={styles.inputWithButton}>
            <input
              type="text"
              value={pubchemInput}
              placeholder={
                INPUT_PLACEHOLDERS[inputMethod] ?? INPUT_PLACEHOLDERS.pubchem
              }
              onChange={e => {
                onPubchemInputChange(e.target.value);
                if (convertError) {
                  onConvertErrorClear();
                }
              }}
              className={styles.pubchemInput}
              disabled={calculationStatus === 'running'}
            />
            <button
              onClick={onXYZConvert}
              className={styles.convertButton}
              disabled={
                isConverting ||
                !pubchemInput.trim() ||
                calculationStatus === 'running'
              }
            >
              {isConverting ? 'Converting...' : 'Convert to XYZ'}
            </button>
          </div>
          {isConverting && (
            <div
              className={`${styles.validationMessage} ${styles.validating} ${styles.convertingMessage}`}
            >
              Converting...
            </div>
          )}
          {convertError && (
            <div
              className={`${styles.validationMessage} ${styles.invalid} ${styles.errorMessage}`}
            >
              ❌ {convertError}
            </div>
          )}
        </div>
        <div className={styles.xyzDirectInput}>
          <h4 className={styles.subsectionTitle}>Direct XYZ Input/Edit</h4>
          <XYZInput onXYZChange={onXYZChange} value={xyzInputValue} />
        </div>
      </section>
    );
  }
);
