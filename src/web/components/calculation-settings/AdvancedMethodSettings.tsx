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
    const isMethodSettingDisabled = !isCalculationEditable(calculationStatus);

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
                disabled={
                  !isCalculationEditable(calculationStatus) || isLoadingParams
                }
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
                disabled={
                  !isCalculationEditable(calculationStatus) || isLoadingParams
                }
              />
            </div>
            {params.calculation_method === 'CASSCF' && (
              <>
                <div className={styles.settingRow}>
                  <label>Energy Convergence Tolerance</label>
                  <input
                    type="number"
                    step={1e-6}
                    value={(params as any).conv_tol ?? 1e-6}
                    onChange={e => {
                      const val = parseFloat(e.target.value);
                      if (!isNaN(val) && val > 0) {
                        onParamChange('conv_tol', val);
                      }
                    }}
                    placeholder="1e-6"
                    className={styles.numberInput}
                    disabled={!isCalculationEditable(calculationStatus)}
                  />
                </div>
                <div className={styles.settingRow}>
                  <label>Gradient Convergence Tolerance</label>
                  <input
                    type="number"
                    step={1e-4}
                    value={(params as any).conv_tol_grad ?? 1e-4}
                    onChange={e => {
                      const val = parseFloat(e.target.value);
                      if (!isNaN(val) && val > 0) {
                        onParamChange('conv_tol_grad', val);
                      }
                    }}
                    placeholder="1e-4"
                    className={styles.numberInput}
                    disabled={!isCalculationEditable(calculationStatus)}
                  />
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
              <div
                className={`${styles.toggleSwitch} ${
                  isMethodSettingDisabled ? styles.toggleDisabled : ''
                }`}
              >
                <label className={styles.toggleLabel} htmlFor="toggle-natorb">
                  Transform to Natural Orbitals in Active Space
                </label>
                <label className={styles.switch}>
                  <input
                    id="toggle-natorb"
                    type="checkbox"
                    checked={(params as any).natorb !== false}
                    onChange={e => onParamChange('natorb', e.target.checked)}
                    disabled={isMethodSettingDisabled}
                    aria-label="Transform to Natural Orbitals in Active Space"
                  />
                  <span className={styles.slider}></span>
                </label>
              </div>
            </div>
          </section>
        )}

        {['DFT', 'HF', 'MP2'].includes(params.calculation_method) &&
          (params as any).optimize_geometry !== false && (
            <section className={styles.calculationSettingsSection}>
              <div className={styles.settingRow}>
                <label>Max Optimization Steps</label>
                <input
                  type="number"
                  value={
                    (params as any).geomopt_maxsteps !== undefined
                      ? (params as any).geomopt_maxsteps
                      : 100
                  }
                  onChange={e =>
                    onParamChange(
                      'geomopt_maxsteps',
                      Math.max(1, Math.min(1000, Number(e.target.value)))
                    )
                  }
                  min={1}
                  max={1000}
                  step={1}
                  className={`${styles.numberInput} ${styles.withSpinner}`}
                  disabled={
                    !isCalculationEditable(calculationStatus) || isLoadingParams
                  }
                />
              </div>
              <div className={styles.settingRow}>
                <label>Energy Convergence (Hartree)</label>
                <input
                  type="number"
                  step={1e-6}
                  value={(params as any).geomopt_conv_energy ?? 1e-6}
                  onChange={e => {
                    const val = parseFloat(e.target.value);
                    if (!isNaN(val) && val > 0) {
                      onParamChange('geomopt_conv_energy', val);
                    }
                  }}
                  placeholder="1e-6"
                  className={styles.numberInput}
                  disabled={
                    !isCalculationEditable(calculationStatus) || isLoadingParams
                  }
                />
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
                disabled={
                  !isCalculationEditable(calculationStatus) || isLoadingParams
                }
              />
            </div>
            <div className={styles.settingRow}>
              <label>TDDFT Method</label>
              <select
                value={(params as any).tddft_method}
                onChange={e => onParamChange('tddft_method', e.target.value)}
                disabled={
                  !isCalculationEditable(calculationStatus) || isLoadingParams
                }
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
              <div
                className={`${styles.toggleSwitch} ${
                  isMethodSettingDisabled ? styles.toggleDisabled : ''
                }`}
              >
                <label
                  className={styles.toggleLabel}
                  htmlFor="toggle-nto-analysis"
                >
                  Natural Transition Orbital Analysis
                </label>
                <label className={styles.switch}>
                  <input
                    id="toggle-nto-analysis"
                    type="checkbox"
                    checked={(params as any).tddft_analyze_nto || false}
                    onChange={e =>
                      onParamChange('tddft_analyze_nto', e.target.checked)
                    }
                    disabled={isMethodSettingDisabled}
                    aria-label="Natural Transition Orbital Analysis"
                  />
                  <span className={styles.slider}></span>
                </label>
              </div>
            </div>
          </section>
        )}

        {(params.calculation_method === 'CCSD' ||
          params.calculation_method === 'CCSD_T') && (
          <section className={styles.calculationSettingsSection}>
            <div className={styles.settingRow}>
              <div
                className={`${styles.toggleSwitch} ${
                  isMethodSettingDisabled ? styles.toggleDisabled : ''
                }`}
              >
                <label
                  className={styles.toggleLabel}
                  htmlFor="toggle-frozen-core"
                >
                  Use Frozen Core Approximation
                </label>
                <label className={styles.switch}>
                  <input
                    id="toggle-frozen-core"
                    type="checkbox"
                    checked={(params as any).frozen_core !== false}
                    onChange={e =>
                      onParamChange('frozen_core', e.target.checked)
                    }
                    disabled={isMethodSettingDisabled}
                    aria-label="Use Frozen Core Approximation"
                  />
                  <span className={styles.slider}></span>
                </label>
              </div>
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
