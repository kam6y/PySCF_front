import { useRef, useState, useEffect, useCallback } from "react";
import { MoleculeViewer, MoleculeViewerRef } from "../components/MoleculeViewer";
import { XYZInput } from "../components/XYZInput";
import { StyleControls } from "../components/StyleControls";
import { StyleSpec } from "../../types/3dmol";
import {
    CalculationParameters,
    CalculationStatus,
    CalculationInstance
} from "../types/calculation";
import { startCalculation, searchPubChem, convertSmilesToXyz } from "../apiClient";
import type { PubChemSearchResponse, SmilesConvertResponse } from "../apiClient";

interface CalculationSettingsPageProps {
    activeCalculation?: CalculationInstance;
    onCalculationUpdate: (updatedCalculation: CalculationInstance) => void;
    onCalculationSuccess: (completedCalculation: CalculationInstance) => void;
    onCalculationRename: (id: string, newName: string) => Promise<void>;
    createNewCalculationFromExisting: (originalCalc: CalculationInstance, newParams: CalculationParameters) => void;
}

export const CalculationSettingsPage = ({
    activeCalculation,
    onCalculationUpdate,
    onCalculationSuccess,
    onCalculationRename,
    createNewCalculationFromExisting,
}: CalculationSettingsPageProps) => {
    const moleculeViewerRef = useRef<MoleculeViewerRef>(null);
    const [calculationStatus, setCalculationStatus] = useState<CalculationStatus>('idle');
    const [calculationError, setCalculationError] = useState<string | null>(null);
    const [inputMethod, setInputMethod] = useState("pubchem");
    const [pubchemInput, setPubchemInput] = useState("");
    const [isConverting, setIsConverting] = useState(false);
    const [convertError, setConvertError] = useState<string | null>(null);

    useEffect(() => {
        const xyz = activeCalculation?.parameters?.xyz;
        if (xyz && xyz.trim() !== "") {
            moleculeViewerRef.current?.loadXYZ(xyz);
        } else {
            moleculeViewerRef.current?.clearModels();
        }
    }, [activeCalculation?.id, activeCalculation?.parameters?.xyz]); // Depend on ID to reload viewer for new instances

    const handleParamChange = useCallback((field: keyof CalculationParameters, value: string | number) => {
        if (!activeCalculation || !onCalculationUpdate) return;

        const isCompleted = activeCalculation.status === 'completed' || activeCalculation.status === 'error';
        const isParamChange = field !== 'molecule_name' && field !== 'xyz';
        const currentParams = activeCalculation.parameters;
        
        // If calculation is completed and a parameter (not name/xyz) is changed, create a new derived calculation
        if (isCompleted && isParamChange) {
            const newParams = { ...currentParams, [field]: value };
            createNewCalculationFromExisting(activeCalculation, newParams);
        } else {
            // Otherwise, just update the current active calculation's state
            const updatedParams = { ...currentParams, [field]: value };
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

    const handleNameBlur = (e: React.FocusEvent<HTMLInputElement>) => {
        if (activeCalculation && activeCalculation.name !== e.target.value.trim() && e.target.value.trim()) {
            onCalculationRename(activeCalculation.id, e.target.value.trim());
        }
    };
    
    if (!activeCalculation) {
        return (
            <div className="calculation-settings-containers">
                <div style={{ textAlign: 'center', padding: '40px', color: '#666', width: '100%' }}>
                    <h2>No Calculation Selected</h2>
                    <p>Please select a calculation from the sidebar, or create a new one using the '+' button.</p>
                </div>
            </div>
        );
    }

    const { parameters: params } = activeCalculation;
    const xyzInputValue = params.xyz || "";
    const hasValidMolecule = !!(params.xyz && params.xyz.trim() !== "");

    const handleStyleChange = (style: StyleSpec) => {
        if (hasValidMolecule) {
            moleculeViewerRef.current?.setStyle(style);
        }
    };

    const handleStartCalculation = async () => {
        if (!hasValidMolecule || !params.xyz) {
            setCalculationError("A valid molecular structure is required.");
            return;
        }

        setCalculationStatus('running');
        setCalculationError(null);

        try {
            const responseData = await startCalculation(params);
            setCalculationStatus('completed');

            if (onCalculationSuccess) {
                const completedInstance: CalculationInstance = {
                    ...activeCalculation,
                    id: responseData.calculation_id,
                    status: 'completed',
                    updatedAt: new Date().toISOString(),
                    parameters: params,
                    results: responseData.calculation_results,
                    workingDirectory: responseData.calculation_results.working_directory,
                };
                onCalculationSuccess(completedInstance);
            }
            alert(`Calculation successful! ID: ${responseData.calculation_id}`);
        } catch (error) {
            setCalculationStatus('error');
            setCalculationError(error instanceof Error ? error.message : 'An unknown error occurred.');
        }
    };

    const handleXYZConvert = async () => {
        if (!pubchemInput.trim()) return;
        setIsConverting(true);
        setConvertError(null);

        try {
            let data: PubChemSearchResponse | SmilesConvertResponse;

            if (inputMethod === 'smiles') {
                data = await convertSmilesToXyz(pubchemInput.trim());
            } else {
                const searchType = /^\d+$/.test(pubchemInput.trim()) ? 'cid' : 'name';
                data = await searchPubChem(pubchemInput.trim(), searchType);
            }

            let moleculeName = params.molecule_name;
            if ('compound_info' in data && data.compound_info?.iupac_name) {
                moleculeName = data.compound_info.iupac_name;
            } else {
                moleculeName = pubchemInput.trim();
            }
            
            const updatedParams = { ...params, xyz: data.xyz, molecule_name: moleculeName };
            onCalculationUpdate({
                ...activeCalculation,
                parameters: updatedParams,
            });

        } catch (error) {
            setConvertError(error instanceof Error ? error.message : 'Unknown error during conversion.');
        } finally {
            setIsConverting(false);
        }
    };
    
    const getInputPlaceholder = () => {
        return inputMethod === 'smiles'
            ? "e.g., CCO for ethanol"
            : "e.g., aspirin, or 2244";
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
                                value={params.molecule_name || ''}
                                onChange={(e) => handleParamChange('molecule_name', e.target.value)}
                                className="molecule-name-input"
                            />
                            {params.molecule_name && (
                                <button onClick={() => handleParamChange('molecule_name', '')} className="clear-molecule-name" aria-label="Clear molecule name"> × </button>
                            )}
                        </div>
                        <div className="computation-settings">
                            <div className="setting-group cpu-setting">
                                <label>CPU cores</label>
                                <div className="cpu-input-container">
                                    <input
                                        type="number"
                                        value={params.cpu_cores || 1}
                                        onChange={(e) => handleParamChange('cpu_cores', Math.max(1, Number(e.target.value)))}
                                        min="1" max="32"
                                        className="cpu-cores-input"
                                    />
                                    <div className="spinner-arrows">
                                        <button type="button" className="spinner-btn up" onClick={() => handleParamChange('cpu_cores', Math.min(32, (params.cpu_cores || 1) + 1))}>▲</button>
                                        <button type="button" className="spinner-btn down" onClick={() => handleParamChange('cpu_cores', Math.max(1, (params.cpu_cores || 1) - 1))}>▼</button>
                                    </div>
                                </div>
                            </div>
                            <div className="setting-group memory-setting">
                                <label>Memory usage</label>
                                <div className="memory-input-container">
                                    <input
                                        type="number"
                                        value={params.memory_mb || 2000}
                                        onChange={(e) => handleParamChange('memory_mb', Math.max(128, Number(e.target.value)))}
                                        min="128"
                                        className="memory-value-input"
                                    />
                                    <span className="memory-unit">MB</span>
                                </div>
                            </div>
                            <button
                                className={`start-calculation-btn ${calculationStatus === 'running' ? 'calculating' : ''}`}
                                onClick={handleStartCalculation}
                                disabled={!hasValidMolecule || calculationStatus === 'running'}
                            >
                                {calculationStatus === 'running' ? '⚛️ Calculating...' : '+ Start Calculation'}
                            </button>
                            {calculationError && (
                                <div className="calculation-error" style={{ marginTop: '10px', color: '#e74c3c', fontSize: '14px' }}>
                                    ❌ {calculationError}
                                </div>
                            )}
                            {calculationStatus === 'running' && (
                                <div className="calculation-status" style={{ marginTop: '10px', color: '#3498db', fontSize: '14px' }}>
                                    ⚛️ Running quantum chemistry calculation... This may take several minutes.
                                </div>
                            )}
                        </div>
                    </div>
                    <div className="calculation-column">
                        <section className="calculation-settings-section">
                            <div className="setting-row">
                                <label>Calculation method</label>
                                <select value={params.calculation_method || 'DFT'} onChange={(e) => handleParamChange('calculation_method', e.target.value)}>
                                    <option value="DFT">DFT</option>
                                    <option value="HF">HF</option>
                                    <option value="MP2">MP2</option>
                                </select>
                            </div>
                            <div className="setting-row">
                                <label>Basis functions</label>
                                <select value={params.basis_function || 'STO-3G'} onChange={(e) => handleParamChange('basis_function', e.target.value)}>
                                    <option value="STO-3G">STO-3G</option>
                                    <option value="3-21G">3-21G</option>
                                    <option value="6-31G">6-31G</option>
                                    <option value="cc-pVDZ">cc-pVDZ</option>
                                </select>
                            </div>
                            <div className="setting-row">
                                <label>Exchange-correlation function</label>
                                <select value={params.exchange_correlation || 'B3LYP'} onChange={(e) => handleParamChange('exchange_correlation', e.target.value)} disabled={params.calculation_method !== 'DFT'}>
                                    <option value="B3LYP">B3LYP</option>
                                    <option value="PBE">PBE</option>
                                    <option value="BP86">BP86</option>
                                    <option value="wB97XD">wB97XD</option>
                                </select>
                            </div>
                            <div className="setting-row">
                                <label>Charges</label>
                                <input
                                    type="number"
                                    value={params.charges || 0}
                                    onChange={(e) => handleParamChange('charges', Number(e.target.value))}
                                    className="number-input with-spinner"
                                />
                            </div>
                            <div className="setting-row">
                                <label>Spin multiplicity (2S+1)</label>
                                <input
                                    type="number"
                                    value={params.spin_multiplicity || 1}
                                    onChange={(e) => handleParamChange('spin_multiplicity', Number(e.target.value))}
                                    min={1} step={1}
                                    className="number-input with-spinner"
                                />
                            </div>
                        </section>
                        <section className="solvent-settings-section">
                            <div className="setting-row">
                                <label>Solvent effect method</label>
                                <select value={params.solvent_method || 'none'} onChange={(e) => handleParamChange('solvent_method', e.target.value)}>
                                    <option value="none">none</option>
                                    <option value="PCM">PCM</option>
                                    <option value="SMD">SMD</option>
                                </select>
                            </div>
                            <div className="setting-row">
                                <label>Solvent</label>
                                <select value={params.solvent || '-'} onChange={(e) => handleParamChange('solvent', e.target.value)} disabled={params.solvent_method === "none"}>
                                    <option value="-">-</option>
                                    <option value="water">Water</option>
                                    <option value="methanol">Methanol</option>
                                    <option value="acetone">Acetone</option>
                                </select>
                            </div>
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
                                        <h3>No molecule loaded</h3>
                                        <p>Enter molecular structure in the right panel to see the 3D visualization</p>
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
                            <input type="radio" name="inputMethod" value="pubchem" checked={inputMethod === "pubchem"} onChange={(e) => setInputMethod(e.target.value)} />
                            <span className="radio-text">Get from PubChem name/CID</span>
                        </label>
                        <label className="radio-option">
                            <input type="radio" name="inputMethod" value="smiles" checked={inputMethod === "smiles"} onChange={(e) => setInputMethod(e.target.value)} />
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
                        />
                        <button onClick={handleXYZConvert} className="convert-button" disabled={isConverting || !pubchemInput.trim()}>
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