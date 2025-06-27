import React, { useState, useEffect } from 'react'
import MoleculeViewer from './components/MoleculeViewer'
import CalculationPanel from './components/CalculationPanel'
import StatusBar from './components/StatusBar'
import ProjectPanel from './components/ProjectPanel'

const App: React.FC = () => {
  const [currentMolecule, setCurrentMolecule] = useState(null)
  const [calculationResults, setCalculationResults] = useState(null)
  const [isCalculating, setIsCalculating] = useState(false)
  const [logs, setLogs] = useState<Array<{ type: string; data: string }>>([])

  useEffect(() => {
    // Set up event listeners for Python backend
    window.electronAPI.onPythonLog((logData) => {
      setLogs(prev => [...prev, logData])
    })

    window.electronAPI.onPythonError((error) => {
      setLogs(prev => [...prev, { type: 'error', data: error }])
    })

    window.electronAPI.onCalculationProgress((progress) => {
      console.log('Calculation progress:', progress)
    })

    return () => {
      window.electronAPI.removeAllListeners('python-log')
      window.electronAPI.removeAllListeners('python-error')
      window.electronAPI.removeAllListeners('calculation-progress')
    }
  }, [])

  const handleCalculation = async (params: any) => {
    if (!currentMolecule) {
      alert('Please select a molecule first')
      return
    }

    setIsCalculating(true)
    setCalculationResults(null)

    try {
      const result = await window.electronAPI.calculate({
        molecule: currentMolecule,
        ...params
      })
      
      if (result.status === 'success') {
        setCalculationResults(result.results)
      } else {
        alert(`Calculation failed: ${result.message}`)
      }
    } catch (error) {
      console.error('Calculation error:', error)
      alert(`Calculation error: ${error}`)
    } finally {
      setIsCalculating(false)
    }
  }

  const handleMoleculeLoad = async (moleculeData: any) => {
    try {
      const result = await window.electronAPI.buildMolecule(moleculeData)
      if (result.status === 'success') {
        setCurrentMolecule(result.molecule)
      } else {
        alert(`Failed to build molecule: ${result.message}`)
      }
    } catch (error) {
      console.error('Molecule building error:', error)
      alert(`Molecule building error: ${error}`)
    }
  }

  return (
    <div className="flex h-screen bg-primary">
      {/* Left Panel - Project and Molecule List */}
      <div className="w-80 bg-secondary border-r border-border flex flex-col">
        <ProjectPanel 
          onMoleculeLoad={handleMoleculeLoad}
          currentMolecule={currentMolecule}
        />
      </div>

      {/* Center Panel - 3D Viewer */}
      <div className="flex-1 flex flex-col">
        <div className="flex-1 bg-tertiary">
          <MoleculeViewer 
            molecule={currentMolecule}
            calculationResults={calculationResults}
          />
        </div>
        
        {/* Status Bar */}
        <StatusBar 
          isCalculating={isCalculating}
          logs={logs}
          onClearLogs={() => setLogs([])}
        />
      </div>

      {/* Right Panel - Calculation Controls */}
      <div className="w-96 bg-secondary border-l border-border">
        <CalculationPanel 
          onCalculate={handleCalculation}
          isCalculating={isCalculating}
          results={calculationResults}
        />
      </div>
    </div>
  )
}

export default App