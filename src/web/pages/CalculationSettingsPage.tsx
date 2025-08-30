// src/web/pages/CalculationSettingsPage.tsx

import { useRef, useState, useEffect, useCallback } from 'react';

import styles from './CalculationSettingsPage.module.css';
import { MoleculeViewerRef } from '../components/MoleculeViewer';
import { XYZInput } from '../components/XYZInput';
import { MoleculeViewerSection } from '../components/MoleculeViewerSection';
import { StyleSpec } from '../../types/3dmol';
import {
  QuantumCalculationRequest,
  CalculationInstance,
  PubChemSearchResponseData,
  SMILESConvertResponseData,
} from '../types/api-types';
import { searchPubChem, convertSmilesToXyz } from '../apiClient';
import { useSupportedParameters } from '../hooks/useCalculationQueries';

interface CalculationSettingsPageProps {
  activeCalculation?: CalculationInstance;
  onCalculationUpdate: (updatedCalculation: CalculationInstance) => void;
  onStartCalculation: (
    params: QuantumCalculationRequest
  ) => Promise<CalculationInstance>;
  onCalculationRename: (id: string, newName: string) => Promise<void>;
  createNewCalculationFromExisting: (
    originalCalc: CalculationInstance,
    newParams: QuantumCalculationRequest
  ) => void;
}

export const CalculationSettingsPage = ({
  activeCalculation,
  onCalculationUpdate,
  onStartCalculation,
  onCalculationRename,
  createNewCalculationFromExisting,
}: CalculationSettingsPageProps) => {
  const moleculeViewerRef = useRef<MoleculeViewerRef>(null);
  const previousCalculationIdRef = useRef<string | null>(null);
  const currentStyleRef = useRef<StyleSpec | null>(null);
  const [calculationError, setCalculationError] = useState<string | null>(null);
  const [inputMethod, setInputMethod] = useState('pubchem');
  const [pubchemInput, setPubchemInput] = useState('');
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);
  const [localName, setLocalName] = useState('');
  const [isEditingName, setIsEditingName] = useState(false);
  const [showAxes, setShowAxes] = useState(false);
  const [showCoordinates, setShowCoordinates] = useState(false);
  const [useAtomicRadii, setUseAtomicRadii] = useState(false);

  // サポートされているパラメータを取得
  const {
    data: supportedParams,
    isLoading: isLoadingParams,
    error: paramsError,
  } = useSupportedParameters();

  useEffect(() => {
    const currentCalculationId = activeCalculation?.id || null;
    const previousCalculationId = previousCalculationIdRef.current;
    const isNewCalculation = currentCalculationId !== previousCalculationId;

    if (activeCalculation) {
      if (isNewCalculation || !isEditingName) {
        setLocalName(
          activeCalculation.name ||
            activeCalculation.parameters?.molecule_name ||
            ''
        );
      }

      const xyz = activeCalculation.parameters?.xyz;
      if (xyz && xyz.trim() !== '') {
        moleculeViewerRef.current?.loadXYZ(xyz);
        // Apply the current style after loading the molecule
        if (currentStyleRef.current) {
          setTimeout(() => {
            moleculeViewerRef.current?.setStyle(currentStyleRef.current!);
          }, 100); // Small delay to ensure the molecule is fully loaded
        }
      } else {
        moleculeViewerRef.current?.clearModels();
      }
    } else {
      setLocalName('');
      moleculeViewerRef.current?.clearModels();
      setIsEditingName(false);
    }

    previousCalculationIdRef.current = currentCalculationId;
  }, [activeCalculation, isEditingName]);

  useEffect(() => {
    moleculeViewerRef.current?.showAxes(showAxes);
  }, [showAxes, activeCalculation]);

  useEffect(() => {
    moleculeViewerRef.current?.showAtomCoordinates(showCoordinates);
  }, [showCoordinates, activeCalculation]);

  useEffect(() => {
    // Re-apply the current style when atomic radii setting changes
    const hasValidMolecule = !!(
      activeCalculation?.parameters?.xyz &&
      activeCalculation.parameters.xyz.trim() !== ''
    );
    if (currentStyleRef.current && hasValidMolecule) {
      const style = { ...currentStyleRef.current };
      if (useAtomicRadii) {
        (style as any)._useAtomicRadii = true;
        (style as any)._baseAtomRadius = 0.3;
      } else {
        (style as any)._useAtomicRadii = false;
      }
      moleculeViewerRef.current?.setStyle(style);
    }
  }, [useAtomicRadii, activeCalculation?.parameters?.xyz]);

  const handleParamChange = useCallback(
    (
      field: keyof QuantumCalculationRequest,
      value: string | number | boolean
    ) => {
      if (!activeCalculation || !onCalculationUpdate) return;

      if (field === 'name') {
        if (activeCalculation.id.startsWith('new-calculation-')) {
          const stringValue = String(value);
          const updatedParams = {
            ...activeCalculation.parameters,
            [field]: stringValue,
          };
          onCalculationUpdate({
            ...activeCalculation,
            name: stringValue,
            parameters: updatedParams,
          });
        }
        return;
      }

      const isCompleted =
        activeCalculation.status === 'completed' ||
        activeCalculation.status === 'error';
      const isParamChange = field !== 'xyz';
      const currentParams = activeCalculation.parameters;

      let processedValue = value;
      if (field === 'solvent' && value === 'custom') {
        processedValue = '78.36';
      }

      // Adjust defaults when calculation method changes
      let adjustedParams = { ...currentParams };
      if (field === 'calculation_method') {
        if (value === 'CCSD' || value === 'CCSD_T') {
          // CCSD defaults: correlation-consistent basis and higher memory
          adjustedParams.basis_function = 'cc-pVDZ';
          adjustedParams.memory_mb = adjustedParams.memory_mb || 4000;
        } else if (value === 'MP2') {
          // MP2 defaults: higher memory
          adjustedParams.memory_mb = adjustedParams.memory_mb || 3000;
        } else if (value === 'TDDFT') {
          // TDDFT defaults
          adjustedParams.basis_function =
            adjustedParams.basis_function || '6-31G(d)';
          adjustedParams.memory_mb = adjustedParams.memory_mb || 2000;
          // Initialize TDDFT parameters if not already set
          (adjustedParams as any).tddft_nstates =
            (adjustedParams as any).tddft_nstates || 10;
          (adjustedParams as any).tddft_method =
            (adjustedParams as any).tddft_method || 'TDDFT';
          (adjustedParams as any).tddft_analyze_nto =
            (adjustedParams as any).tddft_analyze_nto !== undefined
              ? (adjustedParams as any).tddft_analyze_nto
              : false;
        } else {
          // DFT/HF defaults
          adjustedParams.basis_function =
            adjustedParams.basis_function || '6-31G(d)';
          adjustedParams.memory_mb = adjustedParams.memory_mb || 2000;
        }
      }

      const safeParams: QuantumCalculationRequest & { frozen_core?: boolean } =
        {
          xyz: adjustedParams.xyz || '',
          calculation_method: adjustedParams.calculation_method || 'DFT',
          basis_function: adjustedParams.basis_function || '6-31G(d)',
          exchange_correlation: adjustedParams.exchange_correlation || 'B3LYP',
          charges: adjustedParams.charges || 0,
          spin_multiplicity: adjustedParams.spin_multiplicity || 0,
          solvent_method: adjustedParams.solvent_method || 'none',
          solvent: adjustedParams.solvent || '-',
          name:
            (adjustedParams as any).name ||
            (adjustedParams as any).molecule_name ||
            'Unnamed Calculation',
          cpu_cores: adjustedParams.cpu_cores || undefined,
          memory_mb: adjustedParams.memory_mb || undefined,
          tddft_nstates:
            (adjustedParams as any).tddft_nstates !== undefined
              ? (adjustedParams as any).tddft_nstates
              : 10,
          tddft_method: (adjustedParams as any).tddft_method || 'TDDFT',
          tddft_analyze_nto:
            (adjustedParams as any).tddft_analyze_nto !== undefined
              ? (adjustedParams as any).tddft_analyze_nto
              : false,
          frozen_core: (adjustedParams as any).frozen_core !== false, // Default to true
        };

      if (isCompleted && isParamChange) {
        const newParams = { ...safeParams, [field]: processedValue };
        createNewCalculationFromExisting(activeCalculation, newParams);
      } else {
        const updatedParams = { ...safeParams, [field]: processedValue };
        onCalculationUpdate({
          ...activeCalculation,
          parameters: updatedParams,
        });
      }
    },
    [activeCalculation, onCalculationUpdate, createNewCalculationFromExisting]
  );

  // Calculate hasValidMolecule before early return
  const hasValidMolecule = !!(
    activeCalculation?.parameters?.xyz &&
    activeCalculation.parameters.xyz.trim() !== ''
  );

  const handleStyleChange = useCallback(
    (style: StyleSpec) => {
      currentStyleRef.current = style;
      if (hasValidMolecule) {
        moleculeViewerRef.current?.setStyle(style);
      }
    },
    [hasValidMolecule]
  );

  const handleXYZChange = useCallback(
    (xyzData: string, isValid: boolean) => {
      if (isValid && activeCalculation) {
        const isCompleted =
          activeCalculation.status === 'completed' ||
          activeCalculation.status === 'error';
        const currentParams = activeCalculation.parameters;

        const safeParams: QuantumCalculationRequest & {
          frozen_core?: boolean;
        } = {
          xyz: xyzData,
          calculation_method: currentParams.calculation_method || 'DFT',
          basis_function: currentParams.basis_function || '6-31G(d)',
          exchange_correlation: currentParams.exchange_correlation || 'B3LYP',
          charges: currentParams.charges || 0,
          spin_multiplicity: currentParams.spin_multiplicity || 0,
          solvent_method: currentParams.solvent_method || 'none',
          solvent: currentParams.solvent || '-',
          name:
            (currentParams as any).name ||
            (currentParams as any).molecule_name ||
            'Unnamed Calculation',
          cpu_cores: currentParams.cpu_cores || undefined,
          memory_mb: currentParams.memory_mb || undefined,
          tddft_nstates:
            (currentParams as any).tddft_nstates !== undefined
              ? (currentParams as any).tddft_nstates
              : 10,
          tddft_method: (currentParams as any).tddft_method || 'TDDFT',
          tddft_analyze_nto:
            (currentParams as any).tddft_analyze_nto !== undefined
              ? (currentParams as any).tddft_analyze_nto
              : false,
          frozen_core: (currentParams as any).frozen_core !== false, // Default to true
        };

        if (isCompleted) {
          createNewCalculationFromExisting(activeCalculation, safeParams);
        } else {
          const updatedParams = { ...currentParams, xyz: xyzData };
          onCalculationUpdate({
            ...activeCalculation,
            parameters: updatedParams,
          });
        }
      }
    },
    [activeCalculation, onCalculationUpdate, createNewCalculationFromExisting]
  );

  const handleNameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setLocalName(e.target.value);
    setIsEditingName(true);
  };

  const handleNameBlur = async () => {
    const newName = localName.trim();

    if (!activeCalculation || !newName || activeCalculation.name === newName) {
      if (activeCalculation) setLocalName(activeCalculation.name);
      setIsEditingName(false);
      return;
    }

    if (activeCalculation.id.startsWith('new-calculation-')) {
      try {
        onCalculationUpdate({
          ...activeCalculation,
          name: newName,
          parameters: {
            ...activeCalculation.parameters,
            molecule_name: newName,
          },
        });
      } catch (error) {
        console.error('Error updating new calculation name:', error);
        setLocalName(activeCalculation.name || ''); // Revert on error
      }
      setIsEditingName(false);
      return;
    }

    if (
      window.confirm(
        `Are you sure you want to rename this calculation to "${newName}"?`
      )
    ) {
      try {
        await onCalculationRename(activeCalculation.id, newName);
      } catch (error) {
        console.error('Error renaming calculation:', error);
        setLocalName(activeCalculation.name || ''); // Revert on error
        alert(
          `Failed to rename calculation: ${error instanceof Error ? error.message : 'Unknown error'}`
        );
      }
    } else {
      setLocalName(activeCalculation.name || '');
    }
    setIsEditingName(false);
  };

  const handleNameKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter') {
      e.currentTarget.blur();
    }
  };

  const handleStartCalculation = async () => {
    if (
      !activeCalculation ||
      !activeCalculation.parameters?.xyz ||
      !activeCalculation.parameters.xyz.trim()
    ) {
      setCalculationError('A valid molecular structure is required.');
      return;
    }

    const moleculeName = localName.trim();
    if (!moleculeName) {
      setCalculationError('A molecule name is required.');
      return;
    }

    setCalculationError(null);

    const currentParams = activeCalculation.parameters;
    const finalParams: QuantumCalculationRequest & { frozen_core?: boolean } = {
      xyz: currentParams.xyz || '',
      calculation_method: currentParams.calculation_method || 'DFT',
      basis_function: currentParams.basis_function || '6-31G(d)',
      exchange_correlation: currentParams.exchange_correlation || 'B3LYP',
      charges: currentParams.charges || 0,
      spin_multiplicity: currentParams.spin_multiplicity || 0,
      solvent_method: currentParams.solvent_method || 'none',
      solvent: currentParams.solvent || '-',
      name: moleculeName,
      cpu_cores: currentParams.cpu_cores || undefined,
      memory_mb: currentParams.memory_mb || undefined,
      tddft_nstates:
        (currentParams as any).tddft_nstates !== undefined
          ? (currentParams as any).tddft_nstates
          : 10,
      tddft_method: (currentParams as any).tddft_method || 'TDDFT',
      tddft_analyze_nto:
        (currentParams as any).tddft_analyze_nto !== undefined
          ? (currentParams as any).tddft_analyze_nto
          : false,
      frozen_core: (currentParams as any).frozen_core !== false, // Default to true
    };

    try {
      const runningCalculation = await onStartCalculation(finalParams);
      onCalculationUpdate(runningCalculation);
    } catch (error) {
      setCalculationError(
        error instanceof Error ? error.message : 'An unknown error occurred.'
      );
      if (activeCalculation) {
        onCalculationUpdate({ ...activeCalculation, status: 'error' });
      }
    }
  };

  if (!activeCalculation || !activeCalculation.parameters) {
    return (
      <div className={styles.calculationSettingsContainers}>
        <div className={styles.loadingContainer}>
          <h2>Welcome PySCF_front</h2>
          <p>
            Please start a new calculation or select a calculation from the
            sidebar.
          </p>
        </div>
      </div>
    );
  }

  const { parameters: params, status: calculationStatus } = activeCalculation;
  const xyzInputValue = params.xyz || '';

  const handleXYZConvert = async () => {
    if (!pubchemInput.trim()) return;
    setIsConverting(true);
    setConvertError(null);

    try {
      let data: PubChemSearchResponseData | SMILESConvertResponseData;
      let moleculeName = localName;

      if (inputMethod === 'smiles') {
        data = await convertSmilesToXyz(pubchemInput.trim());
        moleculeName = pubchemInput.trim();
      } else {
        const searchType = /^\d+$/.test(pubchemInput.trim()) ? 'cid' : 'name';
        data = await searchPubChem(pubchemInput.trim(), searchType);
        if ('compound_info' in data && data.compound_info?.iupac_name) {
          moleculeName = data.compound_info.iupac_name;
        }
      }

      setLocalName(moleculeName);

      const isCompleted =
        activeCalculation.status === 'completed' ||
        activeCalculation.status === 'error';

      const safeParams: QuantumCalculationRequest & { frozen_core?: boolean } =
        {
          xyz: data.xyz,
          calculation_method: params.calculation_method || 'DFT',
          basis_function: params.basis_function || '6-31G(d)',
          exchange_correlation: params.exchange_correlation || 'B3LYP',
          charges: params.charges || 0,
          spin_multiplicity: params.spin_multiplicity || 0,
          solvent_method: params.solvent_method || 'none',
          solvent: params.solvent || '-',
          name: moleculeName,
          cpu_cores: params.cpu_cores || undefined,
          memory_mb: params.memory_mb || undefined,
          tddft_nstates:
            (params as any).tddft_nstates !== undefined
              ? (params as any).tddft_nstates
              : 10,
          tddft_method: (params as any).tddft_method || 'TDDFT',
          tddft_analyze_nto:
            (params as any).tddft_analyze_nto !== undefined
              ? (params as any).tddft_analyze_nto
              : false,
          frozen_core: (params as any).frozen_core !== false, // Default to true
        };

      if (isCompleted) {
        createNewCalculationFromExisting(activeCalculation, safeParams);
      } else {
        const updatedParams = {
          ...params,
          xyz: data.xyz,
          molecule_name: moleculeName,
        };
        onCalculationUpdate({
          ...activeCalculation,
          name: moleculeName,
          parameters: updatedParams,
        });
      }
    } catch (error) {
      setConvertError(
        error instanceof Error
          ? error.message
          : 'An unknown error occurred during conversion.'
      );
    } finally {
      setIsConverting(false);
    }
  };

  const getInputPlaceholder = () => {
    return inputMethod === 'smiles'
      ? 'e.g., CCO for ethanol'
      : 'e.g., aspirin, or 2244';
  };

  const getCalculationButtonText = () => {
    switch (calculationStatus) {
      case 'running':
        return 'Running...';
      case 'completed':
        return 'Completed!';
      case 'error':
        return 'Error';
      default:
        return '+ Start Calc';
    }
  };

  const isCustomDielectricConstant = (
    solventValue: string | undefined
  ): boolean => {
    if (!solventValue || solventValue === '-') return false;

    const predefinedSolvents = [
      'water',
      'dimethylsulfoxide',
      'n,n-dimethylformamide',
      'nitromethane',
      'methanol',
      'ethanol',
      'acetone',
      'dichloroethane',
      'dichloromethane',
      'tetrahydrofuran',
      'chlorobenzene',
      'chloroform',
      'diethylether',
      'toluene',
      'benzene',
      '1,4-dioxane',
      'cyclohexane',
      'custom',
    ];

    if (predefinedSolvents.includes(solventValue.toLowerCase())) return false;

    const numValue = parseFloat(solventValue);
    return !isNaN(numValue) && numValue > 0;
  };

  const getSolventDisplayValue = (): string => {
    const solventValue = params.solvent || '-';
    return isCustomDielectricConstant(solventValue) ? 'custom' : solventValue;
  };

  const getCustomDielectricValue = (): string => {
    const solventValue = params.solvent || '';
    return isCustomDielectricConstant(solventValue) ? solventValue : '';
  };

  return (
    <div className={styles.calculationSettingsContainers}>
      <div className={styles.calculationSettingsContainer}>
        <div className={styles.calculationHeader}>
          <div className={styles.headerTitleBar}>
            <div className={styles.moleculeNameSection}>
              <input
                type="text"
                placeholder="Molecule name..."
                value={localName}
                onBlur={handleNameBlur}
                onKeyDown={handleNameKeyDown}
                onChange={handleNameChange}
                className={styles.moleculeNameInput}
                disabled={calculationStatus === 'running'}
              />
              {localName && (
                <button
                  onClick={() => setLocalName('')}
                  className={styles.clearMoleculeName}
                  aria-label="Clear molecule name"
                >
                  {' '}
                  ×{' '}
                </button>
              )}
            </div>
            <div className={styles.computationSettings}>
              <div className={styles.cpuSetting}>
                <label>CPU Cores</label>
                <div className={styles.cpuInputContainer}>
                  <input
                    type="number"
                    value={params.cpu_cores || 1}
                    onChange={e =>
                      handleParamChange(
                        'cpu_cores',
                        Math.max(1, Number(e.target.value))
                      )
                    }
                    min="1"
                    max="32"
                    className={styles.cpuCoresInput}
                    disabled={calculationStatus === 'running'}
                  />
                  <div className={styles.spinnerArrows}>
                    <button
                      type="button"
                      className={`${styles.spinnerBtn} ${styles.up}`}
                      onClick={() =>
                        handleParamChange(
                          'cpu_cores',
                          Math.min(32, (params.cpu_cores || 1) + 1)
                        )
                      }
                      disabled={calculationStatus === 'running'}
                    >
                      ▲
                    </button>
                    <button
                      type="button"
                      className={`${styles.spinnerBtn} ${styles.down}`}
                      onClick={() =>
                        handleParamChange(
                          'cpu_cores',
                          Math.max(1, (params.cpu_cores || 1) - 1)
                        )
                      }
                      disabled={calculationStatus === 'running'}
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
                      handleParamChange(
                        'memory_mb',
                        Math.max(128, Number(e.target.value))
                      )
                    }
                    min="128"
                    className={styles.memoryValueInput}
                    disabled={calculationStatus === 'running'}
                  />
                  <span className={styles.memoryUnit}>MB</span>
                </div>
              </div>
              <button
                className={`${styles.startCalculationBtn} ${
                  calculationStatus === 'completed'
                    ? styles.completed
                    : calculationStatus === 'running'
                      ? styles.running
                      : calculationStatus === 'error'
                        ? styles.error
                        : styles.pending
                }`}
                onClick={handleStartCalculation}
                disabled={
                  !hasValidMolecule ||
                  calculationStatus === 'running' ||
                  calculationStatus === 'completed'
                }
              >
                {getCalculationButtonText()}
              </button>
              {calculationError && (
                <div
                  className={styles.calculationError}
                  style={{
                    marginTop: '10px',
                    color: '#e74c3c',
                    fontSize: '14px',
                  }}
                >
                  ❌ {calculationError}
                </div>
              )}
            </div>
          </div>
          <div className={styles.calculationColumn}>
            <section className={styles.calculationSettingsSection}>
              <div className={styles.settingRow}>
                <label>Calculation Method</label>
                <select
                  value={params.calculation_method || 'DFT'}
                  onChange={e =>
                    handleParamChange('calculation_method', e.target.value)
                  }
                  disabled={calculationStatus === 'running'}
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
                  value={params.basis_function || '6-31G(d)'}
                  onChange={e =>
                    handleParamChange('basis_function', e.target.value)
                  }
                  disabled={calculationStatus === 'running'}
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
                <label>Exchange-Correlation Functional</label>
                <select
                  value={params.exchange_correlation || 'B3LYP'}
                  onChange={e =>
                    handleParamChange('exchange_correlation', e.target.value)
                  }
                  disabled={
                    !(
                      params.calculation_method === 'DFT' ||
                      params.calculation_method === 'TDDFT'
                    ) || calculationStatus === 'running'
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
                  value={params.charges || 0}
                  onChange={e =>
                    handleParamChange('charges', Number(e.target.value))
                  }
                  className={`${styles.numberInput} ${styles.withSpinner}`}
                  disabled={calculationStatus === 'running'}
                />
              </div>
              <div className={styles.settingRow}>
                <label>Spin Multiplicity (2S)</label>
                <input
                  type="number"
                  value={params.spin_multiplicity || 0}
                  onChange={e =>
                    handleParamChange(
                      'spin_multiplicity',
                      Number(e.target.value)
                    )
                  }
                  min={0}
                  step={1}
                  className={`${styles.numberInput} ${styles.withSpinner}`}
                  disabled={calculationStatus === 'running'}
                />
              </div>
            </section>
            {params.calculation_method === 'TDDFT' && (
              <section className={styles.calculationSettingsSection}>
                <div className={styles.settingRow}>
                  <label>Number of Excited States</label>
                  <input
                    type="number"
                    value={(params as any).tddft_nstates || 10}
                    onChange={e =>
                      handleParamChange(
                        'tddft_nstates' as any,
                        Math.max(1, Math.min(50, Number(e.target.value)))
                      )
                    }
                    min={1}
                    max={50}
                    step={1}
                    className={`${styles.numberInput} ${styles.withSpinner}`}
                    disabled={calculationStatus === 'running'}
                  />
                </div>
                <div className={styles.settingRow}>
                  <label>TDDFT Method</label>
                  <select
                    value={(params as any).tddft_method || 'TDDFT'}
                    onChange={e =>
                      handleParamChange('tddft_method' as any, e.target.value)
                    }
                    disabled={calculationStatus === 'running'}
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
                        handleParamChange(
                          'tddft_analyze_nto' as any,
                          e.target.checked
                        )
                      }
                      disabled={calculationStatus === 'running'}
                    />
                    Perform Natural Transition Orbital Analysis
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
                      onChange={e =>
                        handleParamChange(
                          'frozen_core' as any,
                          e.target.checked
                        )
                      }
                      disabled={calculationStatus === 'running'}
                    />
                    Use Frozen Core Approximation
                  </label>
                  <div className={styles.frozenCoreHelp}>
                    Freeze core orbitals to reduce computational cost
                    (recommended)
                  </div>
                </div>
              </section>
            )}
            <section className={styles.calculationSettingsSection}>
              <div className={styles.settingRow}>
                <label>Solvent Effect Method</label>
                <select
                  value={params.solvent_method || 'none'}
                  onChange={e =>
                    handleParamChange('solvent_method', e.target.value)
                  }
                  disabled={calculationStatus === 'running'}
                >
                  {isLoadingParams ? (
                    <option value="">Loading...</option>
                  ) : paramsError ? (
                    <option value="">Error loading solvent methods</option>
                  ) : (
                    supportedParams?.solvent_methods
                      ?.map(method => {
                        // Organize methods by category
                        if (method === 'none') {
                          return (
                            <option key={method} value={method}>
                              None
                            </option>
                          );
                        } else if (
                          ['ief-pcm', 'c-pcm', 'cosmo', 'ssvpe'].includes(
                            method
                          )
                        ) {
                          return null; // Handle in PCM optgroup below
                        } else if (method === 'ddcosmo') {
                          return null; // Handle in ddCOSMO optgroup below
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
                              supportedParams.solvent_methods.includes(
                                method
                              ) && (
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
                        {supportedParams.solvent_methods.includes(
                          'ddcosmo'
                        ) && (
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
                  value={getSolventDisplayValue()}
                  onChange={e => handleParamChange('solvent', e.target.value)}
                  disabled={
                    params.solvent_method === 'none' ||
                    calculationStatus === 'running'
                  }
                >
                  {isLoadingParams ? (
                    <option value="">Loading...</option>
                  ) : paramsError ? (
                    <option value="">Error loading solvents</option>
                  ) : (
                    supportedParams?.solvents &&
                    Object.entries(supportedParams.solvents).map(
                      ([group, solvents]) => (
                        <optgroup key={group} label={group}>
                          {solvents.map(solvent => (
                            <option key={solvent.value} value={solvent.value}>
                              {solvent.display}
                            </option>
                          ))}
                        </optgroup>
                      )
                    )
                  )}
                  {!isLoadingParams && !paramsError && (
                    <option value="custom">
                      Custom (Enter dielectric constant below)
                    </option>
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
                      value={getCustomDielectricValue()}
                      onChange={e =>
                        handleParamChange('solvent', e.target.value || '78.36')
                      }
                      className={styles.numberInput}
                      disabled={calculationStatus === 'running'}
                    />
                  </div>
                )}
            </section>
          </div>
        </div>
        <MoleculeViewerSection
          moleculeViewerRef={moleculeViewerRef}
          hasValidMolecule={hasValidMolecule}
          onStyleChange={handleStyleChange}
          showAxes={showAxes}
          onShowAxesChange={setShowAxes}
          showCoordinates={showCoordinates}
          onShowCoordinatesChange={setShowCoordinates}
          useAtomicRadii={useAtomicRadii}
          onUseAtomicRadiiChange={setUseAtomicRadii}
        />
      </div>
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
                onChange={e => setInputMethod(e.target.value)}
                disabled={calculationStatus === 'running'}
              />
              <span className={styles.radioText}>
                Get from PubChem Name/CID
              </span>
            </label>
            <label className={styles.radioOption}>
              <input
                type="radio"
                name="inputMethod"
                value="smiles"
                checked={inputMethod === 'smiles'}
                onChange={e => setInputMethod(e.target.value)}
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
              placeholder={getInputPlaceholder()}
              onChange={e => {
                setPubchemInput(e.target.value);
                if (convertError) setConvertError(null);
              }}
              className={styles.pubchemInput}
              disabled={calculationStatus === 'running'}
            />
            <button
              onClick={handleXYZConvert}
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
          <XYZInput onXYZChange={handleXYZChange} value={xyzInputValue} />
        </div>
      </section>
    </div>
  );
};
