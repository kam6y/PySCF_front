import React, { useEffect, useRef } from 'react'

interface MoleculeViewerProps {
  molecule: any
  calculationResults: any
}

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ 
  molecule, 
  calculationResults 
}) => {
  const viewerRef = useRef<HTMLDivElement>(null)

  useEffect(() => {
    // Initialize 3Dpymol viewer here
    // For now, we'll just show a placeholder
    console.log('Molecule viewer initialized')
  }, [])

  useEffect(() => {
    if (molecule) {
      console.log('Loading molecule:', molecule)
      // Load molecule into 3Dpymol viewer
    }
  }, [molecule])

  useEffect(() => {
    if (calculationResults) {
      console.log('Displaying calculation results:', calculationResults)
      // Update visualization with calculation results
    }
  }, [calculationResults])

  return (
    <div className="w-full h-full relative">
      <div 
        ref={viewerRef} 
        className="w-full h-full bg-primary flex items-center justify-center"
        style={{ background: 'linear-gradient(135deg, #1e293b 0%, #0f172a 100%)' }}
      >
        {molecule ? (
          <div className="text-center">
            <div className="text-6xl mb-4">🧬</div>
            <h3 className="text-xl font-semibold mb-2">Molecule Loaded</h3>
            <div className="text-secondary">
              <p>Atoms: {molecule.natoms || 'Unknown'}</p>
              <p>Charge: {molecule.charge || 0}</p>
              <p>Spin: {molecule.spin || 0}</p>
            </div>
            {calculationResults && (
              <div className="mt-4 p-4 bg-secondary rounded-lg">
                <h4 className="font-medium mb-2">Latest Results</h4>
                <p className="text-sm">
                  Method: {calculationResults.method}
                </p>
                <p className="text-sm">
                  Energy: {calculationResults.energy?.toFixed(6)} Hartree
                </p>
              </div>
            )}
          </div>
        ) : (
          <div className="text-center text-muted">
            <div className="text-6xl mb-4">⚛️</div>
            <h3 className="text-xl font-semibold mb-2">No Molecule Loaded</h3>
            <p>Select or create a molecule to begin</p>
            <p className="text-sm mt-2">
              3Dpymol integration coming soon...
            </p>
          </div>
        )}
      </div>

      {/* Viewer Controls */}
      <div className="absolute top-4 right-4 flex gap-2">
        <button className="btn btn-secondary text-sm">
          Reset View
        </button>
        <button className="btn btn-secondary text-sm">
          Toggle Labels
        </button>
      </div>

      {/* Viewer Info */}
      <div className="absolute bottom-4 left-4 bg-secondary px-3 py-1 rounded text-sm">
        3D Molecular Viewer
      </div>
    </div>
  )
}

export default MoleculeViewer