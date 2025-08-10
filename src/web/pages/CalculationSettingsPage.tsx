import { useRef, useState, useEffect } from "react";
import { MoleculeViewer, MoleculeViewerRef } from "../components/MoleculeViewer";
import { XYZInput } from "../components/XYZInput";
import { StyleControls } from "../components/StyleControls";
import { StyleSpec } from "../../types/3dmol";
import { 
  CalculationParameters, 
  CalculationStatus,
  QuantumCalculationResponse,
  CalculationInstance 
} from "../types/calculation";

interface CalculationSettingsPageProps {
  activeCalculation?: CalculationInstance;
  onCalculationUpdate?: (updatedCalculation: CalculationInstance) => void;
  onCalculationStatusUpdate?: (calculationId: string, status: 'pending' | 'running' | 'completed' | 'error') => Promise<void>;
}

// --- API通信と定数 ---

// APIサーバーのURLを定数として定義
const API_BASE_URL = 'http://127.0.0.1:5000';

/**
 * サーバーのエンドポイントに対して非同期でPOSTリクエストを送信する汎用関数
 * @param endpoint APIエンドポイント (例: '/api/pubchem/search')
 * @param body リクエストボディに含めるオブジェクト
 * @returns 成功した場合はレスポンスのdata部、失敗した場合はエラーをスロー
 */
const fetchAPI = async (endpoint: string, body: object) => {
  const response = await fetch(`${API_BASE_URL}${endpoint}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });

  const data = await response.json();

  if (!response.ok || !data.success) {
    throw new Error(data.error || 'Failed to fetch data from the server.');
  }
  return data.data;
};


// --- メインコンポーネント ---

export const CalculationSettingsPage = ({ 
  activeCalculation, 
  onCalculationUpdate,
  onCalculationStatusUpdate
}: CalculationSettingsPageProps) => {
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
  
  // Track if we're editing an existing calculation
  const isEditingExisting = Boolean(activeCalculation && activeCalculation.parameters);
  
  // Load active calculation parameters into form state
  useEffect(() => {
    if (activeCalculation && activeCalculation.parameters) {
      const params = activeCalculation.parameters;
      setMoleculeName(params.molecule_name || '');
      setCpuCores(params.cpu_cores || 1);
      setMemoryMB(params.memory_mb || 2000);
      setCalculationMethod(params.calculation_method || 'DFT');
      setBasisFunction(params.basis_function || 'STO-3G');
      setExchangeCorrelation(params.exchange_correlation || 'B3LYP');
      setCharges(params.charges || 0);
      setSpinMultiplicity(params.spin_multiplicity || 1);
      setSolventMethod(params.solvent_method || 'none');
      setSolvent(params.solvent || '-');
      if (params.xyz) {
        setXyzInputValue(params.xyz);
        handleXYZChange(params.xyz, true);
      }
    } else {
      // Reset form for new calculation
      setMoleculeName('');
      setCpuCores(1);
      setMemoryMB(2000);
      setCalculationMethod('DFT');
      setBasisFunction('STO-3G');
      setExchangeCorrelation('B3LYP');
      setCharges(0);
      setSpinMultiplicity(1);
      setSolventMethod('none');
      setSolvent('-');
      setXyzInputValue('');
      handleXYZChange('', false);
    }
  }, [activeCalculation]);

  // Molecular input state
  const [inputMethod, setInputMethod] = useState("pubchem");
  const [pubchemInput, setPubchemInput] = useState("");
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);
  const [currentSearchType, setCurrentSearchType] = useState<'name' | 'cid' | null>(null);
  const [xyzInputValue, setXyzInputValue] = useState<string>("");

  // Calculation status state
  const [calculationStatus, setCalculationStatus] = useState<CalculationStatus>('idle');
  const [calculationError, setCalculationError] = useState<string | null>(null);
  const [currentCalculationId, setCurrentCalculationId] = useState<string | null>(null);

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
  const handleStartCalculation = async () => {
    if (!hasValidMolecule || !xyzInputValue.trim()) {
      setCalculationError("Valid molecular structure is required for calculation.");
      return;
    }

    setCalculationStatus('running');
    setCalculationError(null);
    setCurrentCalculationId(null);

    const calculationParams: CalculationParameters = {
      calculation_method: calculationMethod as any,
      basis_function: basisFunction as any,
      exchange_correlation: exchangeCorrelation as any,
      charges,
      spin_multiplicity: spinMultiplicity,
      solvent_method: solventMethod as any,
      solvent,
      xyz: xyzInputValue,
      molecule_name: moleculeName,
      cpu_cores: cpuCores,
      memory_mb: memoryMB
    };

    try {
      console.log("Starting quantum calculation with parameters:", calculationParams);
      
      const response = await fetchAPI('/api/quantum/calculate', calculationParams);
      
      console.log("Calculation completed successfully:", response);
      setCalculationStatus('completed');
      setCurrentCalculationId(response.calculation_id);
      
      // Update parent component about the calculation status
      if (response.calculation_id && onCalculationStatusUpdate) {
        await onCalculationStatusUpdate(response.calculation_id, 'completed');
      }
      
      // Update calculation instance with results
      if (onCalculationUpdate && response.calculation_id) {
        const updatedCalculation: CalculationInstance = {
          id: response.calculation_id,
          name: moleculeName || 'Untitled Calculation',
          status: 'completed',
          createdAt: new Date().toISOString(),
          updatedAt: new Date().toISOString(),
          parameters: calculationParams,
          results: response.calculation_results,
          workingDirectory: response.calculation_results.working_directory
        };
        onCalculationUpdate(updatedCalculation);
      }
      
      // Show success message
      alert(`Calculation completed successfully!\nCalculation ID: ${response.calculation_id}\nHOMO Index: ${response.calculation_results.homo_index}\nLUMO Index: ${response.calculation_results.lumo_index}\nSCF Energy: ${response.calculation_results.scf_energy.toFixed(6)} hartree`);
      
    } catch (error) {
      console.error("Calculation failed:", error);
      setCalculationStatus('error');
      setCalculationError(error instanceof Error ? error.message : 'Unknown error occurred during calculation');
      
      // Update status if we have a calculation ID
      if (currentCalculationId && onCalculationStatusUpdate) {
        try {
          await onCalculationStatusUpdate(currentCalculationId, 'error');
        } catch (err) {
          console.error('Failed to update calculation status to error:', err);
        }
      }
    }
  };

  /**
   * PubChemまたはSMILESから分子構造を取得し、XYZ形式に変換するハンドラ。
   */
  const handleXYZConvert = async () => {
    if (!pubchemInput.trim()) {
      setConvertError("Please enter a value.");
      return;
    }

    setIsConverting(true);
    setConvertError(null);
    setCurrentSearchType(null);

    try {
      let data;
      if (inputMethod === 'smiles') {
        data = await fetchAPI('/api/smiles/convert', { smiles: pubchemInput.trim() });
        setMoleculeName(`Molecule from SMILES`);
      } else { // 'pubchem'
        const searchType = detectSearchType(pubchemInput.trim());
        setCurrentSearchType(searchType);
        data = await fetchAPI('/api/pubchem/search', { query: pubchemInput.trim(), search_type: searchType });
        if (data.compound_info?.iupac_name) {
          setMoleculeName(data.compound_info.iupac_name);
        }
      }

      handleXYZChange(data.xyz, true);

    } catch (error) {
      console.error(`Error converting from ${inputMethod}:`, error);
      setConvertError(error instanceof Error ? error.message : 'An unknown error occurred.');
    } finally {
      setIsConverting(false);
    }
  };

  /**
   * 入力されたクエリがCID（数値のみ）かどうかを判定する。
   */
  const detectSearchType = (query: string): 'name' | 'cid' => {
    const trimmedQuery = query.trim();
    // 純粋な数値（正の整数）の場合はCIDとして扱う
    return /^\d+$/.test(trimmedQuery) ? 'cid' : 'name';
  };

  /**
   * 選択された入力方法に応じて、入力欄のプレースホルダーテキストを返す。
   */
  const getInputPlaceholder = () => {
    if (inputMethod === 'smiles') {
      return "e.g., CCO for ethanol";
    }
    return "e.g., aspirin, or 2244";
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
                    <button type="button" className="spinner-btn up" onClick={() => setCpuCores(prev => Math.min(32, prev + 1))}>▲</button>
                    <button type="button" className="spinner-btn down" onClick={() => setCpuCores(prev => Math.max(1, prev - 1))}>▼</button>
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
                className={`start-calculation-btn ${calculationStatus === 'running' ? 'calculating' : ''}`}
                onClick={handleStartCalculation} 
                disabled={!hasValidMolecule || calculationStatus === 'running'}
              >
                {calculationStatus === 'running' 
                  ? '⚛️ Calculating...' 
                  : isEditingExisting 
                    ? '🔄 Re-run calculation'
                    : '+ Start calculation'
                }
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
              <span className="radio-text">Get from PubChem name/CID</span>
            </label>
            <label className="radio-option">
              <input
                type="radio"
                name="inputMethod"
                value="smiles"
                checked={inputMethod === "smiles"}
                onChange={(e) => setInputMethod(e.target.value)}
              />
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
                if (currentSearchType) setCurrentSearchType(null);
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
          {isConverting && (
            <div className="validation-message validating" style={{ marginTop: '8px' }}>
              {inputMethod === 'smiles' ? 'Converting SMILES...' : 
               currentSearchType === 'cid' ? `Searching by CID...` :
               currentSearchType === 'name' ? `Searching by name...` :
               'Converting...'}
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