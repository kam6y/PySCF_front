import { useRef, useState } from "react";
import { MoleculeViewer, MoleculeViewerRef } from "../components/MoleculeViewer";
import { XYZInput } from "../components/XYZInput";
import { StyleControls } from "../components/StyleControls";
import { StyleSpec } from "../../types/3dmol";

// --- API通信と定数 ---

// APIサーバーのURLを定数として定義
const API_BASE_URL = 'http://127.0.0.1:5000';

/**
 * PubChem APIを介して分子情報を非同期で取得する関数
 * @param query 化合物名またはCID
 * @param searchType 検索タイプ（'name'など）
 * @returns 成功した場合は分子データ、失敗した場合はエラーをスロー
 */
const fetchMoleculeFromPubChem = async (query: string, searchType: string) => {
  const response = await fetch(`${API_BASE_URL}/api/pubchem/search`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ query, search_type: searchType }),
  });

  const data = await response.json();

  if (!response.ok || !data.success) {
    throw new Error(data.error || 'Failed to fetch data from the server.');
  }
  return data.data;
};

// --- メインコンポーネント ---

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
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);
  const [xyzInputValue, setXyzInputValue] = useState<string>("");

  /**
   * XYZデータが変更されたときのハンドラ。分子ビューアを更新する。
   */
  const handleXYZChange = (xyzData: string, isValid: boolean) => {
    setHasValidMolecule(isValid);
    setXyzInputValue(xyzData);

    if (isValid && xyzData.trim()) {
      moleculeViewerRef.current?.loadXYZ(xyzData);
    } else {
      moleculeViewerRef.current?.clearModels();
    }
  };

  /**
   * 分子表示スタイルが変更されたときのハンドラ。
   */
  const handleStyleChange = (style: StyleSpec) => {
    if (hasValidMolecule) {
      moleculeViewerRef.current?.setStyle(style);
    }
  };

  /**
   * 計算開始ボタンが押されたときのハンドラ。
   */
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
      solvent,
      xyz: xyzInputValue,
    });
    // NOTE: ここで実際の計算プロセスを呼び出す処理を実装します
  };

  /**
   * PubChemから分子構造を取得し、XYZ形式に変換するハンドラ。
   */
  const handleXYZConvert = async () => {
    if (!pubchemInput.trim()) {
      setConvertError("Please enter a compound name or CID");
      return;
    }

    setIsConverting(true);
    setConvertError(null);

    try {
      const searchType = inputMethod === 'pubchem' ? 'name' : 'name'; // 将来的にSMILESなども考慮
      const data = await fetchMoleculeFromPubChem(pubchemInput.trim(), searchType);

      if (data.compound_info?.iupac_name) {
        setMoleculeName(data.compound_info.iupac_name);
      }

      // handleXYZChange を呼び出してUI全体を更新
      handleXYZChange(data.xyz, true);

    } catch (error) {
      console.error('Error converting PubChem data:', error);
      setConvertError(error instanceof Error ? error.message : 'An unknown error occurred.');
    } finally {
      setIsConverting(false);
    }
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
              {moleculeName && (
                <button
                  onClick={() => setMoleculeName('')}
                  className="clear-molecule-name"
                  aria-label="Clear molecule name"
                >
                  ×
                </button>
              )}
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
              <button
                className="start-calculation-btn"
                onClick={handleStartCalculation}
                disabled={!hasValidMolecule}
              >
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
                <select value={exchangeCorrelation} onChange={(e) => setExchangeCorrelation(e.target.value)} disabled={calculationMethod !== 'DFT'}>
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
                  className="number-input with-spinner"
                />
              </div>
              <div className="setting-row">
                <label>Spin multiplicity (2S+1)</label>
                <input
                  type="number"
                  value={spinMultiplicity}
                  onChange={(e) => setSpinMultiplicity(Number(e.target.value))}
                  min={1}
                  step={1}
                  className="number-input with-spinner"
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

        {/* Main Content - Two Column Layout */}
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
              <input
                type="radio"
                name="inputMethod"
                value="pubchem"
                checked={inputMethod === "pubchem"}
                onChange={(e) => setInputMethod(e.target.value)}
              />
              <span className="radio-text">Get from PubChem name/ID</span>
            </label>
            <label className="radio-option">
              <input type="radio" name="inputMethod" value="smiles" disabled />
              Get from SMILES
            </label>
          </div>
        </div>
        <div className="pubchem-input-section">
          <div className="input-with-button">
            <input
              type="text"
              value={pubchemInput}
              placeholder="e.g., water, caffeine, or CID"
              onChange={(e) => {
                setPubchemInput(e.target.value);
                if (convertError) setConvertError(null);
              }}
              className="pubchem-input"
            />
            <button
              onClick={handleXYZConvert}
              className="convert-button"
              disabled={isConverting || !pubchemInput.trim()}
            >
              {isConverting ? 'Converting...' : 'Convert to XYZ'}
            </button>
          </div>
          {isConverting && <div className="validation-message validating" style={{ marginTop: '8px' }}>Searching PubChem...</div>}
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