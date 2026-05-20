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
import { clampIntegerInput } from '../../utils/numberInput';

const CHARGE_INPUT = {
  min: -10,
  max: 10,
  fallback: 0,
} as const;

const SPIN_INPUT = {
  min: 0,
  max: 10,
  fallback: 0,
} as const;

interface BasicSettingsSectionProps {
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

export const BasicSettingsSection = React.memo<BasicSettingsSectionProps>(
  ({
    params,
    calculationStatus,
    isLoadingParams,
    paramsError,
    supportedParams,
    onParamChange,
  }) => {
    return (
      <section className={styles.calculationSettingsSection}>
        <div className={styles.settingRow}>
          <label>Calculation Method</label>
          <select
            value={params.calculation_method}
            onChange={e => onParamChange('calculation_method', e.target.value)}
            disabled={
              !isCalculationEditable(calculationStatus) || isLoadingParams
            }
          >
            {isLoadingParams ? (
              <option value="">Loading...</option>
            ) : paramsError ? (
              <option value="">Error loading methods</option>
            ) : (
              supportedParams?.calculation_methods?.map(method => (
                <option key={method} value={method}>
                  {method === 'CCSD_T' ? 'CCSD(T)' : method}
                </option>
              ))
            )}
          </select>
        </div>
        <div className={styles.settingRow}>
          <label>Basis Function</label>
          <select
            value={params.basis_function}
            onChange={e => onParamChange('basis_function', e.target.value)}
            disabled={
              !isCalculationEditable(calculationStatus) || isLoadingParams
            }
          >
            {isLoadingParams ? (
              <option value="">Loading...</option>
            ) : paramsError ? (
              <option value="">Error loading basis functions</option>
            ) : (
              supportedParams?.basis_functions &&
              Object.entries(supportedParams.basis_functions).map(
                ([group, functions]) => (
                  <optgroup key={group} label={group}>
                    {functions.map(func => (
                      <option key={func} value={func}>
                        {func}
                      </option>
                    ))}
                  </optgroup>
                )
              )
            )}
          </select>
        </div>
        <div className={styles.settingRow}>
          <label>Exchange Functional</label>
          <select
            value={params.exchange_correlation || ''}
            onChange={e =>
              onParamChange('exchange_correlation', e.target.value)
            }
            disabled={
              !(
                params.calculation_method === 'DFT' ||
                params.calculation_method === 'TDDFT'
              ) ||
              !isCalculationEditable(calculationStatus) ||
              isLoadingParams
            }
          >
            {isLoadingParams ? (
              <option value="">Loading...</option>
            ) : paramsError ? (
              <option value="">Error loading functionals</option>
            ) : (
              supportedParams?.exchange_correlation &&
              Object.entries(supportedParams.exchange_correlation).map(
                ([group, functionals]) => (
                  <optgroup key={group} label={group}>
                    {functionals.map(functional => (
                      <option key={functional} value={functional}>
                        {functional === 'wB97XD' ? 'ωB97X-D' : functional}
                      </option>
                    ))}
                  </optgroup>
                )
              )
            )}
          </select>
        </div>
        <div className={styles.settingRow}>
          <label>Charge</label>
          <input
            type="number"
            value={params.charges ?? CHARGE_INPUT.fallback}
            onChange={e =>
              onParamChange(
                'charges',
                clampIntegerInput(e.target.value, CHARGE_INPUT)
              )
            }
            min={CHARGE_INPUT.min}
            max={CHARGE_INPUT.max}
            step={1}
            className={`${styles.numberInput} ${styles.withSpinner}`}
            disabled={!isCalculationEditable(calculationStatus)}
          />
        </div>
        <div className={styles.settingRow}>
          <label>Spin (2S)</label>
          <input
            type="number"
            value={params.spin ?? SPIN_INPUT.fallback}
            onChange={e =>
              onParamChange('spin', clampIntegerInput(e.target.value, SPIN_INPUT))
            }
            min={SPIN_INPUT.min}
            max={SPIN_INPUT.max}
            step={1}
            className={`${styles.numberInput} ${styles.withSpinner}`}
            disabled={!isCalculationEditable(calculationStatus)}
          />
        </div>
      </section>
    );
  }
);
