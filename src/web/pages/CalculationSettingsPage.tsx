import { useRef, useState } from "react";
import { MoleculeViewer, MoleculeViewerRef } from "../components/MoleculeViewer";
import { XYZInput } from "../components/XYZInput";
import { StyleControls } from "../components/StyleControls";
import { StyleSpec } from "../../types/3dmol";

export const CalculationSettingsPage = () => {
  const moleculeViewerRef = useRef<MoleculeViewerRef>(null);
  const [hasValidMolecule, setHasValidMolecule] = useState(false);

  // Header state
  const [moleculeName, setMoleculeName] = useState("");
  const [cpuCores, setCpuCores] = useState(1);
  const [memoryMB, setMemoryMB] = useState(2000);

  // Calculation settings state
  const [calculationMethod, setCalculationMethod] = useState("DFT");
  const [basisFunction, setBasisFunction] = useState("STO-3G");
  const [exchangeCorrelation, setExchangeCorrelation] = useState("B3LYP");
  const [charges, setCharges] = useState(0);
  const [spinMultiplicity, setSpinMultiplicity] = useState(1);

  // Solvent settings state
  const [solventMethod, setSolventMethod] = useState("none");
  const [solvent, setSolvent] = useState("-");

  // Molecular input state
  const [inputMethod, setInputMethod] = useState("pubchem");
  const [pubchemInput, setPubchemInput] = useState("");

  const handleXYZChange = (xyzData: string, isValid: boolean) => {
    setHasValidMolecule(isValid);

    if (isValid && xyzData.trim()) {
      moleculeViewerRef.current?.loadXYZ(xyzData);
    } else {
      moleculeViewerRef.current?.clearModels();
    }
  };

  const handleStyleChange = (style: StyleSpec) => {
    if (hasValidMolecule) {
      moleculeViewerRef.current?.setStyle(style);
    }
  };

  const handleStartCalculation = () => {
    console.log("Starting calculation with settings:", {
      moleculeName,
      cpuCores,
      memoryMB,
      calculationMethod,
      basisFunction,
      exchangeCorrelation,
      charges,
      spinMultiplicity,
      solventMethod,
      solvent
    });
  };

  const handleXYZConvert = () => {
    console.log("Converting from PubChem/SMILES:", pubchemInput);
  };

  return (
    <div className="calculation-settings-containers">
      <div className="calculation-settings-container">
        {/* Header Section */}
        <div className="calculation-header">
          <div className="header-title-bar">
            <div className="molecule-name-section">
              <input
                type="text"
                placeholder="Molecule name..."
                value={moleculeName}
                onChange={(e) => setMoleculeName(e.target.value)}
                className="molecule-name-input"
              />
              <button className="clear-molecule-name">×</button>
            </div>
            <div className="computation-settings">
              <div className="setting-group cpu-setting">
                <label>CPU cores</label>
                <div className="cpu-input-container">
                  <input
                    type="number"
                    value={cpuCores}
                    onChange={(e) => setCpuCores(Math.max(1, Number(e.target.value)))}
                    min="1"
                    max="32"
                    className="cpu-cores-input"
                  />
                  <div className="spinner-arrows">
                    <button 
                      type="button" 
                      className="spinner-btn up" 
                      onClick={() => setCpuCores(prev => Math.min(32, prev + 1))}
                    >
                      ▲
                    </button>
                    <button 
                      type="button" 
                      className="spinner-btn down" 
                      onClick={() => setCpuCores(prev => Math.max(1, prev - 1))}
                    >
                      ▼
                    </button>
                  </div>
                </div>
              </div>
              <div className="setting-group memory-setting">
                <label>Memory usage</label>
                <div className="memory-input-container">
                  <input
                    type="number"
                    value={memoryMB}
                    onChange={(e) => setMemoryMB(Math.max(128, Number(e.target.value)))}
                    min="128"
                    className="memory-value-input"
                  />
                  <span className="memory-unit">MB</span>
                </div>
              </div>
              <button className="start-calculation-btn" onClick={handleStartCalculation}>
                + Start calculation
              </button>
            </div>
          </div>
          {/* Left Column - Calculation Settings */}
          <div className="calculation-column">
            {/* Calculation Method Settings */}
            <section className="calculation-settings-section">
              <div className="setting-row">
                <label>Calculation method</label>
                <select value={calculationMethod} onChange={(e) => setCalculationMethod(e.target.value)}>
                  <option value="DFT">DFT</option>
                  <option value="HF">HF</option>
                  <option value="MP2">MP2</option>
                </select>
              </div>
              <div className="setting-row">
                <label>Basis functions</label>
                <select value={basisFunction} onChange={(e) => setBasisFunction(e.target.value)}>
                  <option value="STO-3G">STO-3G</option>
                  <option value="3-21G">3-21G</option>
                  <option value="6-31G">6-31G</option>
                  <option value="cc-pVDZ">cc-pVDZ</option>
                </select>
              </div>
              <div className="setting-row">
                <label>Exchange-correlation function</label>
                <select value={exchangeCorrelation} onChange={(e) => setExchangeCorrelation(e.target.value)}>
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
                  value={charges}
                  onChange={(e) => setCharges(Number(e.target.value))}
                  className="number-input"
                />
              </div>
              <div className="setting-row">
                <label>Spin multiplicity (2S+1)</label>
                <input
                  type="number"
                  value={spinMultiplicity}
                  onChange={(e) => setSpinMultiplicity(Number(e.target.value))}
                  min={1}
                  className="number-input"
                />
              </div>
            </section>
            {/* Solvent Settings */}
            <section className="solvent-settings-section">
              <div className="setting-row">
                <label>Solvent effect method</label>
                <select value={solventMethod} onChange={(e) => setSolventMethod(e.target.value)}>
                  <option value="none">none</option>
                  <option value="PCM">PCM</option>
                  <option value="SMD">SMD</option>
                </select>
              </div>
              <div className="setting-row">
                <label>Solvent</label>
                <select value={solvent} onChange={(e) => setSolvent(e.target.value)} disabled={solventMethod === "none"}>
                  <option value="-">-</option>
                  <option value="water">Water</option>
                  <option value="methanol">Methanol</option>
                  <option value="acetone">Acetone</option>
                </select>
              </div>
            </section>
          </div>
        </div>

        {/* Main Content - Three Column Layout */}
        <div className="main-content">
          {/* Left Column - 3D Molecular Viewer */}
          <div className="left-column">
            <div className="viewer-container">
              <MoleculeViewer
                ref={moleculeViewerRef}
                width={600}
                height={500}
                backgroundColor="white"
                className="main-viewer"
              />

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
          {/* Right Column - Visualization Controls and Molecular Input */}
          <div className="right-column">
            {/* Visualization Style Controls */}
            <section className="visualization-section">
              <StyleControls onStyleChange={handleStyleChange} />
            </section>
          </div>
        </div>
      </div>
      {/* Molecular Structure Input */}
      <section className="molecular-input-section">
        <h3>Molecular Structure Input</h3>

        <div className="input-method-selection">
          <h4>Input Method</h4>
          <label className="radio-option">
            <input
              type="radio"
              name="inputMethod"
              value="pubchem"
              checked={inputMethod === "pubchem"}
              onChange={(e) => setInputMethod(e.target.value)}
            />
            Get from PubChem name/ID
          </label>

          <label className="radio-option">
            <input
              type="radio"
              name="inputMethod"
              value="smiles"
              checked={inputMethod === "smiles"}
              onChange={(e) => setInputMethod(e.target.value)}
            />
            Get from SMILES
          </label>
        </div>

        <div className="pubchem-input-section">
          <label>Compound name or CID</label>
          <div className="pubchem-input-container">
            <input
              type="text"
              value={pubchemInput}
              onChange={(e) => setPubchemInput(e.target.value)}
              className="pubchem-input"
            />
            <button onClick={handleXYZConvert} className="convert-button">
              Convert to XYZ
            </button>
          </div>
        </div>

        <div className="xyz-direct-input">
          <h4>Direct XYZ Input/Edit</h4>
          <XYZInput onXYZChange={handleXYZChange} />
        </div>
      </section>
    </div>
  );
};