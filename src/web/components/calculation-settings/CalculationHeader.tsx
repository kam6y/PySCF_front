import React from 'react';
import styles from '../../pages/CalculationSettingsPage.module.css';
import { CALCULATION_BUTTON_TEXT } from '../../constants/calculationDefaults';
import {
  CalculationInstance,
  CalculationParameters,
  QuantumCalculationRequest,
} from '../../types/api-types';
import { DistributiveKeyOf } from '../../hooks/useCalculationForm';
import { isCalculationEditable } from '../../utils/calculationStatus';

interface CalculationHeaderProps {
  activeCalculation?: CalculationInstance;
  localName: string;
  isEditingName: boolean;
  calculationStatus: CalculationInstance['status'];
  showCpuSettings: boolean;
  hasValidMolecule: boolean;
  calculationError: string | null;
  params: CalculationParameters;
  onNameChange: (e: React.ChangeEvent<HTMLInputElement>) => void;
  onNameBlur: () => void;
  onNameKeyDown: (e: React.KeyboardEvent<HTMLInputElement>) => void;
  onClearName: () => void;
  onParamChange: (
    field: DistributiveKeyOf<QuantumCalculationRequest>,
    value: string | number | boolean
  ) => void;
  onStartCalculation: () => void;
  onCalculationPause: (id: string) => Promise<void>;
  onCalculationResume: (id: string) => Promise<void>;
}

export const CalculationHeader = React.memo<CalculationHeaderProps>(
  ({
    activeCalculation,
    localName,
    isEditingName,
    calculationStatus,
    showCpuSettings,
    hasValidMolecule,
    calculationError,
    params,
    onNameChange,
    onNameBlur,
    onNameKeyDown,
    onClearName,
    onParamChange,
    onStartCalculation,
    onCalculationPause,
    onCalculationResume,
  }) => {
    return (
      <>
        <div className={styles.headerTitleBar}>
          <div className={styles.moleculeNameSection}>
            <input
              type="text"
              placeholder="Molecule name..."
              value={localName}
              onBlur={onNameBlur}
              onKeyDown={onNameKeyDown}
              onChange={onNameChange}
              className={styles.moleculeNameInput}
              disabled={!isCalculationEditable(calculationStatus)}
              data-editing={isEditingName}
            />
            {localName && (
              <button
                onClick={onClearName}
                className={styles.clearMoleculeName}
                aria-label="Clear molecule name"
              >
                {' '}
                ×{' '}
              </button>
            )}
          </div>
          <div className={styles.computationSettings}>
            {showCpuSettings && (
              <>
                <div className={styles.cpuSetting}>
                  <label>CPU Cores</label>
                  <div className={styles.cpuInputContainer}>
                    <input
                      type="number"
                      value={params.cpu_cores || 1}
                      onChange={e =>
                        onParamChange(
                          'cpu_cores',
                          Math.max(1, Number(e.target.value))
                        )
                      }
                      min="1"
                      max="32"
                      className={styles.cpuCoresInput}
                      disabled={!isCalculationEditable(calculationStatus)}
                    />
                    <div className={styles.spinnerArrows}>
                      <button
                        type="button"
                        className={`${styles.spinnerBtn} ${styles.up}`}
                        onClick={() =>
                          onParamChange(
                            'cpu_cores',
                            Math.min(32, (params.cpu_cores || 1) + 1)
                          )
                        }
                        disabled={!isCalculationEditable(calculationStatus)}
                      >
                        ▲
                      </button>
                      <button
                        type="button"
                        className={`${styles.spinnerBtn} ${styles.down}`}
                        onClick={() =>
                          onParamChange(
                            'cpu_cores',
                            Math.max(1, (params.cpu_cores || 1) - 1)
                          )
                        }
                        disabled={!isCalculationEditable(calculationStatus)}
                      >
                        ▼
                      </button>
                    </div>
                  </div>
                </div>
                <div className={styles.memorySetting}>
                  <label>Memory Usage</label>
                  <div className={styles.memoryInputContainer}>
                    <input
                      type="number"
                      value={params.memory_mb || 2000}
                      onChange={e =>
                        onParamChange(
                          'memory_mb',
                          Math.max(128, Number(e.target.value))
                        )
                      }
                      min="128"
                      className={styles.memoryValueInput}
                      disabled={!isCalculationEditable(calculationStatus)}
                    />
                    <span className={styles.memoryUnit}>MB</span>
                  </div>
                </div>
              </>
            )}
            {calculationStatus === 'running' || calculationStatus === 'pausing' ? (
              <button
                className={`${styles.pauseBtn} ${
                  calculationStatus === 'pausing' ? styles.pausing : ''
                }`}
                onClick={() => {
                  if (!activeCalculation) return;
                  void onCalculationPause(activeCalculation.id);
                }}
                disabled={calculationStatus === 'pausing'}
              >
                <span className={styles.pauseIcon}>⏸</span>
                {calculationStatus === 'pausing' ? 'Pausing...' : 'Pause'}
              </button>
            ) : calculationStatus === 'paused' ? (
              <button
                className={styles.resumeBtn}
                onClick={() => {
                  if (!activeCalculation) return;
                  void onCalculationResume(activeCalculation.id);
                }}
              >
                <span className={styles.resumeIcon}>▶</span>
                Resume
              </button>
            ) : calculationStatus === 'waiting' ? (
              <button
                className={`${styles.startCalculationBtn} ${styles.waiting}`}
                disabled
              >
                Waiting...
              </button>
            ) : (
              <button
                className={`${styles.startCalculationBtn} ${
                  calculationStatus === 'completed'
                    ? styles.completed
                    : calculationStatus === 'error'
                      ? styles.error
                      : styles.pending
                }`}
                onClick={onStartCalculation}
                disabled={!hasValidMolecule || calculationStatus === 'completed'}
              >
                {CALCULATION_BUTTON_TEXT[calculationStatus] ?? '+ Start Calc'}
              </button>
            )}
          </div>
        </div>
        {calculationError && (
          <div className={styles.calculationError}>❌ {calculationError}</div>
        )}
      </>
    );
  }
);
