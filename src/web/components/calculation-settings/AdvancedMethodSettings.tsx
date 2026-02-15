import React from 'react';
import styles from '../../pages/CalculationSettingsPage.module.css';
import {
  CalculationParameters,
  CalculationStatus,
  QuantumCalculationRequest,
  SupportedParametersResponseData,
} from '../../types/api-types';
import { DistributiveKeyOf } from '../../hooks/useCalculationForm';
import { isCalculationEditable } from '../../utils/calculationStatus';

interface AdvancedMethodSettingsProps {
  params: CalculationParameters;
  calculationStatus: CalculationStatus;
  isLoadingParams: boolean;
  paramsError: unknown;
  supportedParams?: SupportedParametersResponseData;
  onParamChange: (
    field: DistributiveKeyOf<QuantumCalculationRequest>,
    value: string | number | boolean
  ) => void;
}

export const AdvancedMethodSettings = React.memo<AdvancedMethodSettingsProps>(
  ({
    params,
    calculationStatus,
    isLoadingParams,
    paramsError,
    supportedParams,
    onParamChange,
  }) => {
    return (
      <>
        {(params.calculation_method === 'CASCI' ||
          params.calculation_method === 'CASSCF') && (
          <section className={styles.calculationSettingsSection}>
            <div className={styles.settingRow}>
              <label>Number of Active Orbitals</label>
              <input
                type="number"
                value={(params as any).ncas}
                onChange={e =>
                  onParamChange(
                    'ncas',
                    Math.max(1, Math.min(20, Number(e.target.value)))
                  )
                }
                min={1}
                max={20}
                step={1}
                className={`${styles.numberInput} ${styles.withSpinner}`}
                disabled={!isCalculationEditable(calculationStatus) || isLoadingParams}
              />
            </div>
            <div className={styles.settingRow}>
              <label>Number of Active Electrons</label>
              <input
                type="number"
                value={(params as any).nelecas}
                onChange={e =>
                  onParamChange(
                    'nelecas',
                    Math.max(1, Math.min(40, Number(e.target.value)))
                  )
                }
                min={1}
                max={40}
                step={1}
                className={`${styles.numberInput} ${styles.withSpinner}`}
                disabled={!isCalculationEditable(calculationStatus) || isLoadingParams}
              />
            </div>
            {params.calculation_method === 'CASSCF' && (
              <>
                <div className={styles.settingRow}>
                  <label>Energy Convergence Tolerance</label>
                  <select
                    value={(params as any).conv_tol ?? 1e-6}
                    onChange={e =>
                      onParamChange('conv_tol', parseFloat(e.target.value))
                    }
                    disabled={!isCalculationEditable(calculationStatus)}
                  >
                    <option value={1e-5}>1e-5 (loose)</option>
                    <option value={1e-6}>1e-6 (normal)</option>
                    <option value={1e-7}>1e-7 (tight)</option>
                    <option value={1e-8}>1e-8 (very tight)</option>
                  </select>
                </div>
                <div className={styles.settingRow}>
                  <label>Gradient Convergence Tolerance</label>
                  <select
                    value={(params as any).conv_tol_grad ?? 1e-4}
                    onChange={e =>
                      onParamChange('conv_tol_grad', parseFloat(e.target.value))
                    }
                    disabled={!isCalculationEditable(calculationStatus)}
                  >
                    <option value={1e-3}>1e-3 (loose)</option>
                    <option value={1e-4}>1e-4 (normal)</option>
                    <option value={1e-5}>1e-5 (tight)</option>
                    <option value={1e-6}>1e-6 (very tight)</option>
                  </select>
                </div>
                <div className={styles.settingRow}>
                  <label>CASSCF Max Macro Iterations</label>
                  <input
                    type="number"
                    value={(params as any).max_cycle_macro}
                    onChange={e =>
                      onParamChange(
                        'max_cycle_macro',
                        Math.max(1, Math.min(200, Number(e.target.value)))
                      )
                    }
                    min={1}
                    max={200}
                    step={1}
                    className={`${styles.numberInput} ${styles.withSpinner}`}
                    disabled={!isCalculationEditable(calculationStatus)}
                  />
                </div>
              </>
            )}
            <div className={styles.settingRow}>
              <label>CI Max Micro Iterations</label>
              <input
                type="number"
                value={
                  (params as any).max_cycle_micro !== undefined
                    ? (params as any).max_cycle_micro
                    : 3
                }
                onChange={e =>
                  onParamChange(
                    'max_cycle_micro',
                    Math.max(1, Math.min(100, Number(e.target.value)))
                  )
                }
                min={1}
                max={100}
                step={1}
                className={`${styles.numberInput} ${styles.withSpinner}`}
                disabled={!isCalculationEditable(calculationStatus)}
              />
            </div>
            <div className={styles.settingRow}>
              <label>
                <input
                  type="checkbox"
                  checked={(params as any).natorb !== false}
                  onChange={e => onParamChange('natorb', e.target.checked)}
                  disabled={!isCalculationEditable(calculationStatus)}
                />
                Transform to Natural Orbitals in Active Space
              </label>
            </div>
          </section>
        )}

        {params.calculation_method === 'TDDFT' && (
          <section className={styles.calculationSettingsSection}>
            <div className={styles.settingRow}>
              <label>Number of Excited States</label>
              <input
                type="number"
                value={(params as any).tddft_nstates}
                onChange={e =>
                  onParamChange(
                    'tddft_nstates',
                    Math.max(1, Math.min(50, Number(e.target.value)))
                  )
                }
                min={1}
                max={50}
                step={1}
                className={`${styles.numberInput} ${styles.withSpinner}`}
                disabled={!isCalculationEditable(calculationStatus) || isLoadingParams}
              />
            </div>
            <div className={styles.settingRow}>
              <label>TDDFT Method</label>
              <select
                value={(params as any).tddft_method}
                onChange={e => onParamChange('tddft_method', e.target.value)}
                disabled={!isCalculationEditable(calculationStatus) || isLoadingParams}
              >
                {isLoadingParams ? (
                  <option value="">Loading...</option>
                ) : paramsError ? (
                  <option value="">Error loading TDDFT methods</option>
                ) : (
                  supportedParams?.tddft_methods?.map(method => (
                    <option key={method} value={method}>
                      {method === 'TDDFT'
                        ? 'Full TDDFT'
                        : method === 'TDA'
                          ? 'Tamm-Dancoff Approximation (TDA)'
                          : method}
                    </option>
                  ))
                )}
              </select>
            </div>
            <div className={styles.settingRow}>
              <label>
                <input
                  type="checkbox"
                  checked={(params as any).tddft_analyze_nto || false}
                  onChange={e =>
                    onParamChange('tddft_analyze_nto', e.target.checked)
                  }
                  disabled={!isCalculationEditable(calculationStatus)}
                />
                Natural Transition Orbital Analysis
              </label>
            </div>
          </section>
        )}

        {(params.calculation_method === 'CCSD' ||
          params.calculation_method === 'CCSD_T') && (
          <section className={styles.calculationSettingsSection}>
            <div className={styles.settingRow}>
              <label>
                <input
                  type="checkbox"
                  checked={(params as any).frozen_core !== false}
                  onChange={e => onParamChange('frozen_core', e.target.checked)}
                  disabled={!isCalculationEditable(calculationStatus)}
                />
                Use Frozen Core Approximation
              </label>
              <div className={styles.frozenCoreHelp}>
                Freeze core orbitals to reduce computational cost (recommended)
              </div>
            </div>
          </section>
        )}
      </>
    );
  }
);
