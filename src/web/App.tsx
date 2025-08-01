import { useRef, useState } from "react";
import "./App.css";
import { MoleculeViewer, MoleculeViewerRef } from "./components/MoleculeViewer";
import { XYZInput } from "./components/XYZInput";
import { StyleControls } from "./components/StyleControls";
import { Header } from "./components/Header";
import { Sidebar } from "./components/Sidebar";
import { DropdownMenu, DropdownOption } from "./components/DropdownMenu";
import { StyleSpec } from "../types/3dmol";

export const App = () => {
  const moleculeViewerRef = useRef<MoleculeViewerRef>(null);
  const [hasValidMolecule, setHasValidMolecule] = useState(false);
  
  // New state for UI components
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const [selectedDropdownOption, setSelectedDropdownOption] = useState<DropdownOption>('calculation-settings');

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

  const handleZoomToFit = () => {
    if (hasValidMolecule) {
      moleculeViewerRef.current?.zoomToFit();
    }
  };

  const handleTakeScreenshot = () => {
    if (hasValidMolecule) {
      const screenshot = moleculeViewerRef.current?.takeScreenshot();
      if (screenshot) {
        // Create a download link for the screenshot
        const link = document.createElement('a');
        link.download = 'molecule-screenshot.png';
        link.href = screenshot;
        link.click();
      }
    }
  };

  // New event handlers for UI components
  const handleSidebarToggle = () => {
    setIsSidebarOpen(!isSidebarOpen);
  };

  const handleSidebarClose = () => {
    setIsSidebarOpen(false);
  };

  const handleDropdownToggle = () => {
    setIsDropdownOpen(!isDropdownOpen);
  };

  const handleDropdownClose = () => {
    setIsDropdownOpen(false);
  };

  const handleDropdownOptionSelect = (option: DropdownOption) => {
    setSelectedDropdownOption(option);
    // Handle navigation based on selected option
    switch (option) {
      case 'calculation-settings':
        // Currently the main view - no navigation needed
        break;
      case 'calculation-results':
        // TODO: Navigate to calculation results view
        console.log('Navigate to calculation results');
        break;
      case 'draw-molecule':
        // TODO: Navigate to draw molecule view
        console.log('Navigate to draw molecule');
        break;
    }
  };

  const handlePlusClick = () => {
    // Navigate to Calculation page
    console.log('Navigate to Calculation page');
    // TODO: Implement navigation to calculation page
  };

  return (
    <div className="app-container">
      {/* New Header */}
      <Header
        onSidebarToggle={handleSidebarToggle}
        onDropdownToggle={handleDropdownToggle}
        onPlusClick={handlePlusClick}
        isDropdownOpen={isDropdownOpen}
        isSidebarOpen={isSidebarOpen}
      />

      {/* Dropdown Menu */}
      <DropdownMenu
        isOpen={isDropdownOpen}
        selectedOption={selectedDropdownOption}
        onOptionSelect={handleDropdownOptionSelect}
        onClose={handleDropdownClose}
      />

      {/* Sidebar */}
      <Sidebar
        isOpen={isSidebarOpen}
        onClose={handleSidebarClose}
      />

      {/* Main Content */}
      <main className={`app-main ${isSidebarOpen ? 'sidebar-open' : ''}`}>
        {/* Left Panel - Input and Controls */}
        <div className="left-panel">
          <section className="input-section">
            <XYZInput onXYZChange={handleXYZChange} />
          </section>

          <section className="controls-section">
            <StyleControls 
              onStyleChange={handleStyleChange}
            />
          </section>

          {/* Action Buttons */}
          <section className="actions-section">
            <div className="action-buttons">
              <button 
                onClick={handleZoomToFit}
                disabled={!hasValidMolecule}
                className="action-button"
              >
                🔍 Zoom to Fit
              </button>
              <button 
                onClick={handleTakeScreenshot}
                disabled={!hasValidMolecule}
                className="action-button"
              >
                📷 Screenshot
              </button>
            </div>
          </section>
        </div>

        {/* Right Panel - 3D Viewer */}
        <div className="right-panel">
          <section className="viewer-section">
            <div className="viewer-header">
              <h2>3D Molecular Structure</h2>
              {hasValidMolecule && (
                <div className="viewer-status">
                  ✅ Molecule loaded successfully
                </div>
              )}
            </div>
            
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
                    <p>Enter valid XYZ coordinates in the left panel to see the 3D structure</p>
                    <p>Try loading a sample molecule to get started!</p>
                  </div>
                </div>
              )}
            </div>
          </section>
        </div>
      </main>

      {/* Footer */}
      <footer className="app-footer">
        <p>Powered by 3Dmol.js • Built with React and TypeScript</p>
      </footer>
    </div>
  );
};