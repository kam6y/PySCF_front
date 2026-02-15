import React from 'react';
import styles from '../../pages/CalculationSettingsPage.module.css';
import { isCustomDielectricConstant } from '../../constants/calculationDefaults';
import {
  CalculationParameters,
  CalculationStatus,
  QuantumCalculationRequest,
  SupportedParametersResponseData,
} from '../../types/api-types';
import { DistributiveKeyOf } from '../../hooks/useCalculationForm';
import { isCalculationEditable } from '../../utils/calculationStatus';

interface SolventSettingsSectionProps {
  params: CalculationParameters;
  calculationStatus: CalculationStatus;
  isLoadingParams: boolean;
  paramsError: unknown;
  supportedParams?: SupportedParametersResponseData;
  solventDisplayValue: string;
  customDielectricValue: string;
  onParamChange: (
    field: DistributiveKeyOf<QuantumCalculationRequest>,
    value: string | number | boolean
  ) => void;
}

export const SolventSettingsSection = React.memo<SolventSettingsSectionProps>(
  ({
    params,
    calculationStatus,
    isLoadingParams,
    paramsError,
    supportedParams,
    solventDisplayValue,
    customDielectricValue,
    onParamChange,
  }) => {
    return (
      <section className={styles.calculationSettingsSection}>
        <div className={styles.settingRow}>
          <label>Solvent Effect Method</label>
          <select
            value={params.solvent_method || 'none'}
            onChange={e => onParamChange('solvent_method', e.target.value)}
            disabled={!isCalculationEditable(calculationStatus)}
          >
            {isLoadingParams ? (
              <option value="">Loading...</option>
            ) : paramsError ? (
              <option value="">Error loading solvent methods</option>
            ) : (
              supportedParams?.solvent_methods
                ?.map(method => {
                  if (method === 'none') {
                    return (
                      <option key={method} value={method}>
                        None
                      </option>
                    );
                  }
                  if (
                    ['ief-pcm', 'c-pcm', 'cosmo', 'ssvpe'].includes(method) ||
                    method === 'ddcosmo'
                  ) {
                    return null;
                  }
                  return (
                    <option key={method} value={method}>
                      {method}
                    </option>
                  );
                })
                .filter(Boolean)
            )}
            {!isLoadingParams &&
              !paramsError &&
              supportedParams?.solvent_methods && (
                <>
                  <optgroup label="PCM Methods">
                    {['ief-pcm', 'c-pcm', 'cosmo', 'ssvpe'].map(
                      method =>
                        supportedParams.solvent_methods.includes(method) && (
                          <option key={method} value={method}>
                            {method === 'ief-pcm'
                              ? 'IEF-PCM'
                              : method === 'c-pcm'
                                ? 'C-PCM'
                                : method === 'cosmo'
                                  ? 'COSMO'
                                  : method === 'ssvpe'
                                    ? 'SS(V)PE'
                                    : method}
                          </option>
                        )
                    )}
                  </optgroup>
                  {supportedParams.solvent_methods.includes('ddcosmo') && (
                    <optgroup label="ddCOSMO Method">
                      <option value="ddcosmo">ddCOSMO</option>
                    </optgroup>
                  )}
                </>
              )}
          </select>
        </div>
        <div className={styles.settingRow}>
          <label>Solvent(dielectric constant)</label>
          <select
            value={solventDisplayValue}
            onChange={e => onParamChange('solvent', e.target.value)}
            disabled={
              params.solvent_method === 'none' ||
              !isCalculationEditable(calculationStatus)
            }
          >
            {isLoadingParams ? (
              <option value="">Loading...</option>
            ) : paramsError ? (
              <option value="">Error loading solvents</option>
            ) : (
              supportedParams?.solvents &&
              Object.entries(supportedParams.solvents).map(([group, solvents]) => (
                <optgroup key={group} label={group}>
                  {solvents.map(solvent => (
                    <option key={solvent.value} value={solvent.value}>
                      {solvent.display}
                    </option>
                  ))}
                </optgroup>
              ))
            )}
            {!isLoadingParams && !paramsError && (
              <option value="custom">Custom (Enter dielectric constant below)</option>
            )}
          </select>
        </div>
        {(params.solvent === 'custom' ||
          isCustomDielectricConstant(params.solvent)) &&
          params.solvent_method !== 'none' && (
            <div className={styles.settingRow}>
              <label>Custom Dielectric Constant</label>
              <input
                type="number"
                min="0"
                value={customDielectricValue}
                onChange={e =>
                  onParamChange('solvent', e.target.value || '78.36')
                }
                className={styles.numberInput}
                disabled={!isCalculationEditable(calculationStatus)}
              />
            </div>
          )}
      </section>
    );
  }
);
