import React, { useRef, useState } from "react";
import "./App.css";
import { MoleculeViewer, MoleculeViewerRef } from "./components/MoleculeViewer";
import { XYZInput } from "./components/XYZInput";
import { StyleControls } from "./components/StyleControls";
import { StyleSpec } from "../types/3dmol";

export const App = () => {
  const moleculeViewerRef = useRef<MoleculeViewerRef>(null);
  const [hasValidMolecule, setHasValidMolecule] = useState(false);
  const [currentXYZ, setCurrentXYZ] = useState("");

  const handleXYZChange = (xyzData: string, isValid: boolean) => {
    setCurrentXYZ(xyzData);
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

  return (
    <div className="app-container">
      {/* Header */}
      <header className="app-header">
        <h1>PySCF Molecular Visualizer</h1>
        <p>Enter XYZ coordinates to visualize molecular structures in 3D</p>
      </header>

      {/* Main Content */}
      <main className="app-main">
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