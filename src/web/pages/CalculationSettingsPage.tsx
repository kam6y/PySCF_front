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
          adjustedParams.basis_function = adjustedParams.basis_function || '6-31G(d)';
          adjustedParams.memory_mb = adjustedParams.memory_mb || 2000;
          // Initialize TDDFT parameters if not already set
          (adjustedParams as any).tddft_nstates = (adjustedParams as any).tddft_nstates || 10;
          (adjustedParams as any).tddft_method = (adjustedParams as any).tddft_method || 'TDDFT';
          (adjustedParams as any).tddft_analyze_nto = (adjustedParams as any).tddft_analyze_nto !== undefined 
            ? (adjustedParams as any).tddft_analyze_nto 
            : false;
        } else {
          // DFT/HF defaults
          adjustedParams.basis_function = adjustedParams.basis_function || '6-31G(d)';
          adjustedParams.memory_mb = adjustedParams.memory_mb || 2000;
        }
      }

      const safeParams: QuantumCalculationRequest & { frozen_core?: boolean } = {
        xyz: adjustedParams.xyz || '',
        calculation_method: adjustedParams.calculation_method || 'DFT',
        basis_function: adjustedParams.basis_function || '6-31G(d)',
        exchange_correlation: adjustedParams.exchange_correlation || 'B3LYP',
        charges: adjustedParams.charges || 0,
        spin_multiplicity: adjustedParams.spin_multiplicity || 1,
        solvent_method: adjustedParams.solvent_method || 'none',
        solvent: adjustedParams.solvent || '-',
        name:
          (adjustedParams as any).name ||
          (adjustedParams as any).molecule_name ||
          'Unnamed Calculation',
        cpu_cores: adjustedParams.cpu_cores || undefined,
        memory_mb: adjustedParams.memory_mb || undefined,
        tddft_nstates: (adjustedParams as any).tddft_nstates !== undefined 
          ? (adjustedParams as any).tddft_nstates 
          : 10,
        tddft_method: (adjustedParams as any).tddft_method || 'TDDFT',
        tddft_analyze_nto: (adjustedParams as any).tddft_analyze_nto !== undefined 
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

        const safeParams: QuantumCalculationRequest & { frozen_core?: boolean } = {
          xyz: xyzData,
          calculation_method: currentParams.calculation_method || 'DFT',
          basis_function: currentParams.basis_function || '6-31G(d)',
          exchange_correlation: currentParams.exchange_correlation || 'B3LYP',
          charges: currentParams.charges || 0,
          spin_multiplicity: currentParams.spin_multiplicity || 1,
          solvent_method: currentParams.solvent_method || 'none',
          solvent: currentParams.solvent || '-',
          name:
            (currentParams as any).name ||
            (currentParams as any).molecule_name ||
            'Unnamed Calculation',
          cpu_cores: currentParams.cpu_cores || undefined,
          memory_mb: currentParams.memory_mb || undefined,
          tddft_nstates: (currentParams as any).tddft_nstates !== undefined 
            ? (currentParams as any).tddft_nstates 
            : 10,
          tddft_method: (currentParams as any).tddft_method || 'TDDFT',
          tddft_analyze_nto: (currentParams as any).tddft_analyze_nto !== undefined 
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
      spin_multiplicity: currentParams.spin_multiplicity || 1,
      solvent_method: currentParams.solvent_method || 'none',
      solvent: currentParams.solvent || '-',
      name: moleculeName,
      cpu_cores: currentParams.cpu_cores || undefined,
      memory_mb: currentParams.memory_mb || undefined,
      tddft_nstates: (currentParams as any).tddft_nstates !== undefined 
        ? (currentParams as any).tddft_nstates 
        : 10,
      tddft_method: (currentParams as any).tddft_method || 'TDDFT',
      tddft_analyze_nto: (currentParams as any).tddft_analyze_nto !== undefined 
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
          <h2>Loading Calculation...</h2>
          <p>Please wait or select a calculation from the sidebar.</p>
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

      const safeParams: QuantumCalculationRequest & { frozen_core?: boolean } = {
        xyz: data.xyz,
        calculation_method: params.calculation_method || 'DFT',
        basis_function: params.basis_function || '6-31G(d)',
        exchange_correlation: params.exchange_correlation || 'B3LYP',
        charges: params.charges || 0,
        spin_multiplicity: params.spin_multiplicity || 1,
        solvent_method: params.solvent_method || 'none',
        solvent: params.solvent || '-',
        name: moleculeName,
        cpu_cores: params.cpu_cores || undefined,
        memory_mb: params.memory_mb || undefined,
        tddft_nstates: (params as any).tddft_nstates !== undefined 
          ? (params as any).tddft_nstates 
          : 10,
        tddft_method: (params as any).tddft_method || 'TDDFT',
        tddft_analyze_nto: (params as any).tddft_analyze_nto !== undefined 
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
        return 'Completed';
      case 'error':
        return 'Error';
      default:
        return '+ Start';
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
                className={`${styles.startCalculationBtn} ${calculationStatus || 'pending'}`}
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
                  <option value="DFT">DFT</option>
                  <option value="HF">HF</option>
                  <option value="MP2">MP2</option>
                  <option value="CCSD">CCSD</option>
                  <option value="CCSD_T">CCSD(T)</option>
                  <option value="TDDFT">TDDFT</option>
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
                  <optgroup label="Minimal">
                    <option value="STO-3G">STO-3G</option>
                    <option value="3-21G">3-21G</option>
                  </optgroup>
                  <optgroup label="Pople Style">
                    <option value="6-31G">6-31G</option>
                    <option value="6-31G(d)">6-31G(d)</option>
                    <option value="6-31+G(d,p)">6-31+G(d,p)</option>
                    <option value="6-311G(d,p)">6-311G(d,p)</option>
                    <option value="6-311++G(d,p)">6-311++G(d,p)</option>
                  </optgroup>
                  <optgroup label="Correlation Consistent">
                    <option value="cc-pVDZ">cc-pVDZ</option>
                    <option value="cc-pVTZ">cc-pVTZ</option>
                    <option value="cc-pVQZ">cc-pVQZ</option>
                    <option value="aug-cc-pVDZ">aug-cc-pVDZ</option>
                    <option value="aug-cc-pVTZ">aug-cc-pVTZ</option>
                  </optgroup>
                  <optgroup label="def2">
                    <option value="def2-SVP">def2-SVP</option>
                    <option value="def2-TZVP">def2-TZVP</option>
                  </optgroup>
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
                  <optgroup label="Hybrid">
                    <option value="B3LYP">B3LYP</option>
                    <option value="PBE0">PBE0</option>
                    <option value="M06-2X">M06-2X</option>
                    <option value="CAM-B3LYP">CAM-B3LYP</option>
                    <option value="wB97XD">ωB97X-D</option>
                  </optgroup>
                  <optgroup label="GGA">
                    <option value="PBE">PBE</option>
                    <option value="BLYP">BLYP</option>
                    <option value="BP86">BP86</option>
                    <option value="PW91">PW91</option>
                  </optgroup>
                  <optgroup label="Meta-GGA">
                    <option value="M06">M06</option>
                    <option value="M06-L">M06-L</option>
                    <option value="TPSS">TPSS</option>
                  </optgroup>
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
                <label>Spin Multiplicity (2S+1)</label>
                <input
                  type="number"
                  value={params.spin_multiplicity || 1}
                  onChange={e =>
                    handleParamChange(
                      'spin_multiplicity',
                      Number(e.target.value)
                    )
                  }
                  min={1}
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
                    <option value="TDDFT">Full TDDFT</option>
                    <option value="TDA">
                      Tamm-Dancoff Approximation (TDA)
                    </option>
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
            {(params.calculation_method === 'CCSD' || params.calculation_method === 'CCSD_T') && (
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
                    Freeze core orbitals to reduce computational cost (recommended)
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
                  <option value="none">None</option>
                  <optgroup label="PCM Methods">
                    <option value="ief-pcm">IEF-PCM</option>
                    <option value="c-pcm">C-PCM</option>
                    <option value="cosmo">COSMO</option>
                    <option value="ssvpe">SS(V)PE</option>
                  </optgroup>
                  <optgroup label="ddCOSMO Method">
                    <option value="ddcosmo">ddCOSMO</option>
                  </optgroup>
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
                  <optgroup label="Highly Polar">
                    <option value="water">Water (78.36)</option>
                    <option value="dimethylsulfoxide">
                      Dimethylsulfoxide (46.83)
                    </option>
                    <option value="n,n-dimethylformamide">
                      N,N-Dimethylformamide (37.22)
                    </option>
                    <option value="nitromethane">Nitromethane (36.56)</option>
                  </optgroup>
                  <optgroup label="Protic Solvents">
                    <option value="methanol">Methanol (32.61)</option>
                    <option value="ethanol">Ethanol (24.85)</option>
                  </optgroup>
                  <optgroup label="Polar Aprotic">
                    <option value="acetone">Acetone (20.49)</option>
                    <option value="dichloroethane">
                      Dichloroethane (10.13)
                    </option>
                    <option value="dichloromethane">
                      Dichloromethane (8.93)
                    </option>
                    <option value="tetrahydrofuran">
                      Tetrahydrofuran (7.43)
                    </option>
                    <option value="chlorobenzene">Chlorobenzene (5.70)</option>
                  </optgroup>
                  <optgroup label="Moderately Polar">
                    <option value="chloroform">Chloroform (4.71)</option>
                    <option value="diethylether">Diethylether (4.24)</option>
                  </optgroup>
                  <optgroup label="Nonpolar">
                    <option value="toluene">Toluene (2.37)</option>
                    <option value="benzene">Benzene (2.27)</option>
                    <option value="1,4-dioxane">1,4-Dioxane (2.21)</option>
                    <option value="cyclohexane">Cyclohexane (2.02)</option>
                  </optgroup>
                  <option value="custom">
                    Custom (Enter dielectric constant below)
                  </option>
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
              <span className={styles.radioText}>Get from PubChem Name/CID</span>
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
