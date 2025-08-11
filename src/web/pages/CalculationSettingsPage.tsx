// src/web/pages/CalculationSettingsPage.tsx

import { useRef, useState, useEffect, useCallback } from "react";
import { MoleculeViewer, MoleculeViewerRef } from "../components/MoleculeViewer";
import { XYZInput } from "../components/XYZInput";
import { StyleControls } from "../components/StyleControls";
import { StyleSpec } from "../../types/3dmol";
import {
  QuantumCalculationRequest,
  CalculationInstance,
  PubChemSearchResponseData,
  SMILESConvertResponseData
} from "../types/api-types";
import { searchPubChem, convertSmilesToXyz } from "../apiClient";
import { useCalculationPolling } from "../hooks/useCalculationPolling";

interface CalculationSettingsPageProps {
  activeCalculation?: CalculationInstance;
  onCalculationUpdate: (updatedCalculation: CalculationInstance) => void;
  onStartCalculation: (params: QuantumCalculationRequest) => Promise<CalculationInstance>;
  onCalculationRename: (id: string, newName: string) => Promise<void>;
  createNewCalculationFromExisting: (originalCalc: CalculationInstance, newParams: QuantumCalculationRequest) => void;
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
  const [calculationError, setCalculationError] = useState<string | null>(null);
  const [inputMethod, setInputMethod] = useState("pubchem");
  const [pubchemInput, setPubchemInput] = useState("");
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);
  const [localName, setLocalName] = useState("");
  const [isEditingName, setIsEditingName] = useState(false);

  const { startPolling, stopPolling } = useCalculationPolling({
    calculationId: activeCalculation?.id || null,
    onUpdate: onCalculationUpdate,
    onError: (error: string) => setCalculationError(error)
  });


  useEffect(() => {
    const currentCalculationId = activeCalculation?.id || null;
    const previousCalculationId = previousCalculationIdRef.current;
    const isNewCalculation = currentCalculationId !== previousCalculationId;

    if (activeCalculation) {
      // Update localName only if:
      // 1. It's a new calculation (ID changed), OR
      // 2. User is not currently editing the name
      if (isNewCalculation || !isEditingName) {
        setLocalName(activeCalculation.name || activeCalculation.parameters?.molecule_name || "");
      }

      const xyz = activeCalculation.parameters?.xyz;
      if (xyz && xyz.trim() !== "") {
        moleculeViewerRef.current?.loadXYZ(xyz);
      } else {
        moleculeViewerRef.current?.clearModels();
      }
      // If calculation is running, start polling for its status
      if (activeCalculation.status === 'running') {
        startPolling();
      } else {
        stopPolling();
      }
    } else {
      setLocalName("");
      moleculeViewerRef.current?.clearModels();
      stopPolling();
      setIsEditingName(false); // Reset editing state when no calculation
    }

    // Update the previous calculation ID reference
    previousCalculationIdRef.current = currentCalculationId;
  }, [activeCalculation, isEditingName, startPolling, stopPolling]);

  const handleParamChange = useCallback((field: keyof QuantumCalculationRequest, value: string | number) => {
    if (!activeCalculation || !onCalculationUpdate) return;

    // For name, we handle it differently via the dedicated name input field
    // This prevents conflicts between the main name input and parameter updates
    if (field === 'name') {
      // For new calculations, allow direct parameter updates
      if (activeCalculation.id.startsWith('new-calculation-')) {
        const stringValue = String(value);
        const updatedParams = { ...activeCalculation.parameters, [field]: stringValue };
        onCalculationUpdate({
          ...activeCalculation,
          name: stringValue,
          parameters: updatedParams,
        });
      }
      return;
    }

    const isCompleted = activeCalculation.status === 'completed' || activeCalculation.status === 'error';
    const isParamChange = field !== 'xyz';
    const currentParams = activeCalculation.parameters;

    // Handle special case for solvent field when "custom" is selected
    let processedValue = value;
    if (field === 'solvent' && value === 'custom') {
      // When "custom" is selected from dropdown, set a default dielectric constant
      processedValue = '78.36';
    }

    // Ensure required fields have default values for QuantumCalculationRequest
    const safeParams: QuantumCalculationRequest = {
      xyz: currentParams.xyz || '',
      calculation_method: currentParams.calculation_method || 'DFT',
      basis_function: currentParams.basis_function || '6-31G(d)',
      exchange_correlation: currentParams.exchange_correlation || 'B3LYP',
      charges: currentParams.charges || 0,
      spin_multiplicity: currentParams.spin_multiplicity || 1,
      solvent_method: currentParams.solvent_method || 'none',
      solvent: currentParams.solvent || '-',
      name: (currentParams as any).name || (currentParams as any).molecule_name || 'Unnamed Calculation',
      cpu_cores: currentParams.cpu_cores || undefined,
      memory_mb: currentParams.memory_mb || undefined
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
  }, [activeCalculation, onCalculationUpdate, createNewCalculationFromExisting]);

  const handleXYZChange = useCallback((xyzData: string, isValid: boolean) => {
    if (isValid && activeCalculation) {
      const updatedParams = { ...activeCalculation.parameters, xyz: xyzData };
      onCalculationUpdate({ ...activeCalculation, parameters: updatedParams });
    }
  }, [activeCalculation, onCalculationUpdate]);

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

    // Handle new calculations (temporary IDs)
    if (activeCalculation.id.startsWith('new-calculation-')) {
      try {
        onCalculationUpdate({
          ...activeCalculation,
          name: newName,
          parameters: { ...activeCalculation.parameters, molecule_name: newName },
        });
      } catch (error) {
        console.error('Error updating new calculation name:', error);
        setLocalName(activeCalculation.name || ''); // Revert on error
      }
      setIsEditingName(false);
      return;
    }

    // Handle existing calculations (require confirmation and API call)
    if (window.confirm(`Are you sure you want to rename this calculation to "${newName}"?`)) {
      try {
        await onCalculationRename(activeCalculation.id, newName);
        // onCalculationRename should handle the UI updates via the callback chain
      } catch (error) {
        console.error('Error renaming calculation:', error);
        setLocalName(activeCalculation.name || ''); // Revert on error
        alert(`Failed to rename calculation: ${error instanceof Error ? error.message : 'Unknown error'}`);
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
    if (!activeCalculation || !activeCalculation.parameters?.xyz || !activeCalculation.parameters.xyz.trim()) {
      setCalculationError("A valid molecular structure is required.");
      return;
    }

    // Ensure molecule name is provided
    const moleculeName = localName.trim();
    if (!moleculeName) {
      setCalculationError("A molecule name is required.");
      return;
    }

    setCalculationError(null);
    
    // Ensure all parameters are properly set with current UI values and defaults
    const currentParams = activeCalculation.parameters;
    const finalParams: QuantumCalculationRequest = {
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
      memory_mb: currentParams.memory_mb || undefined
    };

    try {
      const runningCalculation = await onStartCalculation(finalParams);
      onCalculationUpdate(runningCalculation);
    } catch (error) {
      setCalculationError(error instanceof Error ? error.message : 'An unknown error occurred.');
      if (activeCalculation) {
        onCalculationUpdate({ ...activeCalculation, status: 'error' });
      }
    }
  };

  // Add a defensive check for activeCalculation and its parameters
  if (!activeCalculation || !activeCalculation.parameters) {
    return (
      <div className="calculation-settings-containers">
        <div style={{ textAlign: 'center', padding: '40px', color: '#666', width: '100%' }}>
          <h2>Loading Calculation...</h2>
          <p>Please wait or select a calculation from the sidebar.</p>
        </div>
      </div>
    );
  }

  const { parameters: params, status: calculationStatus } = activeCalculation;
  const xyzInputValue = params.xyz || "";
  const hasValidMolecule = !!(params.xyz && params.xyz.trim() !== "");

  const handleStyleChange = (style: StyleSpec) => {
    if (hasValidMolecule) {
      moleculeViewerRef.current?.setStyle(style);
    }
  };

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

      // Update local state and calculation parameters
      setLocalName(moleculeName);
      const updatedParams = { ...params, xyz: data.xyz, molecule_name: moleculeName };
      onCalculationUpdate({
        ...activeCalculation,
        name: moleculeName,
        parameters: updatedParams,
      });

    } catch (error) {
      setConvertError(error instanceof Error ? error.message : 'An unknown error occurred during conversion.');
    } finally {
      setIsConverting(false);
    }
  };

  const getInputPlaceholder = () => {
    return inputMethod === 'smiles'
      ? "e.g., CCO for ethanol"
      : "e.g., aspirin, or 2244";
  };

  const getCalculationButtonText = () => {
    switch (calculationStatus) {
      case 'running':
        return '⚛️ Running...';
      case 'completed':
        return '✅ Completed';
      case 'error':
        return '❌ Error';
      default:
        return '+ Start Calculation';
    }
  };

  const getCalculationButtonClass = () => {
    return `start-calculation-btn ${calculationStatus || 'pending'}`;
  };

  // Helper function to check if solvent value is a custom numeric dielectric constant
  const isCustomDielectricConstant = (solventValue: string | undefined): boolean => {
    if (!solventValue || solventValue === '-') return false;
    
    // List of predefined solvent names
    const predefinedSolvents = [
      'water', 'dimethylsulfoxide', 'n,n-dimethylformamide', 'nitromethane',
      'methanol', 'ethanol', 'acetone', 'dichloroethane', 'dichloromethane',
      'tetrahydrofuran', 'chlorobenzene', 'chloroform', 'diethylether',
      'toluene', 'benzene', '1,4-dioxane', 'cyclohexane', 'custom'
    ];
    
    // If it's a predefined solvent name, it's not custom
    if (predefinedSolvents.includes(solventValue.toLowerCase())) return false;
    
    // Check if it's a valid numeric string (dielectric constant)
    const numValue = parseFloat(solventValue);
    return !isNaN(numValue) && numValue > 0;
  };

  // Get the display value for solvent dropdown
  const getSolventDisplayValue = (): string => {
    const solventValue = params.solvent || '-';
    return isCustomDielectricConstant(solventValue) ? 'custom' : solventValue;
  };

  // Get the numeric value for custom dielectric input
  const getCustomDielectricValue = (): string => {
    const solventValue = params.solvent || '';
    return isCustomDielectricConstant(solventValue) ? solventValue : '';
  };

  return (
    <div className="calculation-settings-containers">
      <div className="calculation-settings-container">
        <div className="calculation-header">
          <div className="header-title-bar">
            <div className="molecule-name-section">
              <input
                type="text"
                placeholder="Molecule name..."
                value={localName}
                onBlur={handleNameBlur}
                onKeyDown={handleNameKeyDown}
                onChange={handleNameChange}
                className="molecule-name-input"
                disabled={calculationStatus === 'running'}
              />
              {localName && (
                <button onClick={() => setLocalName('')} className="clear-molecule-name" aria-label="Clear molecule name"> × </button>
              )}
            </div>
            <div className="computation-settings">
              <div className="setting-group cpu-setting">
                <label>CPU Cores</label>
                <div className="cpu-input-container">
                  <input
                    type="number"
                    value={params.cpu_cores || 1}
                    onChange={(e) => handleParamChange('cpu_cores', Math.max(1, Number(e.target.value)))}
                    min="1" max="32"
                    className="cpu-cores-input"
                    disabled={calculationStatus === 'running'}
                  />
                  <div className="spinner-arrows">
                    <button type="button" className="spinner-btn up" onClick={() => handleParamChange('cpu_cores', Math.min(32, (params.cpu_cores || 1) + 1))} disabled={calculationStatus === 'running'}>▲</button>
                    <button type="button" className="spinner-btn down" onClick={() => handleParamChange('cpu_cores', Math.max(1, (params.cpu_cores || 1) - 1))} disabled={calculationStatus === 'running'}>▼</button>
                  </div>
                </div>
              </div>
              <div className="setting-group memory-setting">
                <label>Memory Usage</label>
                <div className="memory-input-container">
                  <input
                    type="number"
                    value={params.memory_mb || 2000}
                    onChange={(e) => handleParamChange('memory_mb', Math.max(128, Number(e.target.value)))}
                    min="128"
                    className="memory-value-input"
                    disabled={calculationStatus === 'running'}
                  />
                  <span className="memory-unit">MB</span>
                </div>
              </div>
              <button
                className={getCalculationButtonClass()}
                onClick={handleStartCalculation}
                disabled={!hasValidMolecule || calculationStatus === 'running' || calculationStatus === 'completed'}
              >
                {getCalculationButtonText()}
              </button>
              {calculationError && (
                <div className="calculation-error" style={{ marginTop: '10px', color: '#e74c3c', fontSize: '14px' }}>
                  ❌ {calculationError}
                </div>
              )}
            </div>
          </div>
          <div className="calculation-column">
            <section className="calculation-settings-section">
              <div className="setting-row">
                <label>Calculation Method</label>
                <select value={params.calculation_method || 'DFT'} onChange={(e) => handleParamChange('calculation_method', e.target.value)} disabled={calculationStatus === 'running'}>
                  <option value="DFT">DFT</option>
                  <option value="HF">HF</option>
                  <option value="MP2">MP2</option>
                </select>
              </div>
              <div className="setting-row">
                <label>Basis Function</label>
                <select value={params.basis_function || '6-31G(d)'} onChange={(e) => handleParamChange('basis_function', e.target.value)} disabled={calculationStatus === 'running'}>
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
              <div className="setting-row">
                <label>Exchange-Correlation Functional</label>
                <select value={params.exchange_correlation || 'B3LYP'} onChange={(e) => handleParamChange('exchange_correlation', e.target.value)} disabled={params.calculation_method !== 'DFT' || calculationStatus === 'running'}>
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
              <div className="setting-row">
                <label>Charge</label>
                <input
                  type="number"
                  value={params.charges || 0}
                  onChange={(e) => handleParamChange('charges', Number(e.target.value))}
                  className="number-input with-spinner"
                  disabled={calculationStatus === 'running'}
                />
              </div>
              <div className="setting-row">
                <label>Spin Multiplicity (2S+1)</label>
                <input
                  type="number"
                  value={params.spin_multiplicity || 1}
                  onChange={(e) => handleParamChange('spin_multiplicity', Number(e.target.value))}
                  min={1} step={1}
                  className="number-input with-spinner"
                  disabled={calculationStatus === 'running'}
                />
              </div>
            </section>
            <section className="solvent-settings-section">
              <div className="setting-row">
                <label>Solvent Effect Method</label>
                <select value={params.solvent_method || 'none'} onChange={(e) => handleParamChange('solvent_method', e.target.value)} disabled={calculationStatus === 'running'}>
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
              <div className="setting-row">
                <label>Solvent(dielectric constant)</label>
                <select value={getSolventDisplayValue()} onChange={(e) => handleParamChange('solvent', e.target.value)} disabled={params.solvent_method === "none" || calculationStatus === 'running'}>
                  <option value="-">-</option>
                  <optgroup label="Highly Polar">
                    <option value="water">Water (78.36)</option>
                    <option value="dimethylsulfoxide">Dimethylsulfoxide (46.83)</option>
                    <option value="n,n-dimethylformamide">N,N-Dimethylformamide (37.22)</option>
                    <option value="nitromethane">Nitromethane (36.56)</option>
                  </optgroup>
                  <optgroup label="Protic Solvents">
                    <option value="methanol">Methanol (32.61)</option>
                    <option value="ethanol">Ethanol (24.85)</option>
                  </optgroup>
                  <optgroup label="Polar Aprotic">
                    <option value="acetone">Acetone (20.49)</option>
                    <option value="dichloroethane">Dichloroethane (10.13)</option>
                    <option value="dichloromethane">Dichloromethane (8.93)</option>
                    <option value="tetrahydrofuran">Tetrahydrofuran (7.43)</option>
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
                  <option value="custom">Custom (Enter dielectric constant below)</option>
                </select>
              </div>
              {(params.solvent === 'custom' || isCustomDielectricConstant(params.solvent)) && params.solvent_method !== 'none' && (
                <div className="setting-row">
                  <label>Custom Dielectric Constant</label>
                  <input
                    type="number"
                    min="1"
                    placeholder="e.g., 15.5"
                    value={getCustomDielectricValue()}
                    onChange={(e) => handleParamChange('solvent', e.target.value || '78.36')}
                    className="number-input"
                    disabled={calculationStatus === 'running'}
                  />
                </div>
              )}
            </section>
          </div>
        </div>
        <div className="main-content">
          <div className="left-column">
            <div className="viewer-container">
              <MoleculeViewer ref={moleculeViewerRef} width={600} height={500} backgroundColor="white" className="main-viewer" />
              {!hasValidMolecule && (
                <div className="viewer-placeholder">
                  <div className="placeholder-content">
                    <h3>No Molecule Loaded</h3>
                    <p>Enter a molecular structure in the right panel to see the 3D visualization</p>
                  </div>
                </div>
              )}
            </div>
          </div>
          <div className="right-column">
            <section className="visualization-section">
              <StyleControls onStyleChange={handleStyleChange} />
            </section>
          </div>
        </div>
      </div>
      <section className="molecular-input-section">
        <h3>Molecular Structure Input</h3>
        <div className="input-method-selection">
          <div className="radio-options">
            <label className="radio-option">
              <input type="radio" name="inputMethod" value="pubchem" checked={inputMethod === "pubchem"} onChange={(e) => setInputMethod(e.target.value)} disabled={calculationStatus === 'running'} />
              <span className="radio-text">Get from PubChem Name/CID</span>
            </label>
            <label className="radio-option">
              <input type="radio" name="inputMethod" value="smiles" checked={inputMethod === "smiles"} onChange={(e) => setInputMethod(e.target.value)} disabled={calculationStatus === 'running'} />
              <span className="radio-text">Get from SMILES</span>
            </label>
          </div>
        </div>
        <div className="pubchem-input-section">
          <div className="input-with-button">
            <input
              type="text"
              value={pubchemInput}
              placeholder={getInputPlaceholder()}
              onChange={(e) => {
                setPubchemInput(e.target.value);
                if (convertError) setConvertError(null);
              }}
              className="pubchem-input"
              disabled={calculationStatus === 'running'}
            />
            <button onClick={handleXYZConvert} className="convert-button" disabled={isConverting || !pubchemInput.trim() || calculationStatus === 'running'}>
              {isConverting ? 'Converting...' : 'Convert to XYZ'}
            </button>
          </div>
          {isConverting && (
            <div className="validation-message validating" style={{ marginTop: '8px' }}>
              Converting...
            </div>
          )}
          {convertError && <div className="validation-message invalid" style={{ marginTop: '8px' }}>❌ {convertError}</div>}
        </div>
        <div className="xyz-direct-input">
          <h4 className="subsection-title">Direct XYZ Input/Edit</h4>
          <XYZInput onXYZChange={handleXYZChange} value={xyzInputValue} />
        </div>
      </section>
    </div>
  );
};