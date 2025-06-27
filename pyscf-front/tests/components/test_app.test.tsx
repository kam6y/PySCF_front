/**
 * Tests for App.tsx main React component
 * Tests application state management and component integration
 */

import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import '@testing-library/jest-dom'
import App from '../../src/renderer/App'

// Mock child components
jest.mock('../../src/renderer/components/MoleculeViewer', () => {
  return function MockMoleculeViewer({ molecule, calculationResults }: any) {
    return (
      <div data-testid="molecule-viewer">
        <div>Molecule: {molecule ? `${molecule.natoms} atoms` : 'None'}</div>
        <div>Results: {calculationResults ? calculationResults.method : 'None'}</div>
      </div>
    )
  }
})

jest.mock('../../src/renderer/components/CalculationPanel', () => {
  return function MockCalculationPanel({ onCalculate, isCalculating, results }: any) {
    return (
      <div data-testid="calculation-panel">
        <button 
          onClick={() => onCalculate({ method: 'HF', basis: 'STO-3G' })}
          disabled={isCalculating}
        >
          {isCalculating ? 'Calculating...' : 'Start Calculation'}
        </button>
        <div>Results: {results ? 'Available' : 'None'}</div>
      </div>
    )
  }
})

jest.mock('../../src/renderer/components/StatusBar', () => {
  return function MockStatusBar({ isCalculating, logs, onClearLogs }: any) {
    return (
      <div data-testid="status-bar">
        <div>Status: {isCalculating ? 'Calculating' : 'Ready'}</div>
        <div>Logs: {logs.length} entries</div>
        <button onClick={onClearLogs}>Clear Logs</button>
      </div>
    )
  }
})

jest.mock('../../src/renderer/components/ProjectPanel', () => {
  return function MockProjectPanel({ onMoleculeLoad, currentMolecule }: any) {
    const loadTestMolecule = () => {
      onMoleculeLoad({
        type: 'coordinates',
        coordinates: [['H', 0, 0, 0], ['H', 0.74, 0, 0]],
        charge: 0,
        spin: 0
      })
    }
    
    return (
      <div data-testid="project-panel">
        <button onClick={loadTestMolecule}>Load Test Molecule</button>
        <div>Current: {currentMolecule ? 'Loaded' : 'None'}</div>
      </div>
    )
  }
})

// Mock electron API
const mockElectronAPI = {
  calculate: jest.fn(),
  buildMolecule: jest.fn(),
  onPythonLog: jest.fn(),
  onPythonError: jest.fn(),
  onCalculationProgress: jest.fn(),
  removeAllListeners: jest.fn()
}

// Setup global window.electronAPI
Object.defineProperty(window, 'electronAPI', {
  value: mockElectronAPI
})

describe('App Component', () => {
  beforeEach(() => {
    jest.clearAllMocks()
  })

  describe('Initial Rendering', () => {
    test('should render all main components', () => {
      render(<App />)
      
      expect(screen.getByTestId('project-panel')).toBeInTheDocument()
      expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
      expect(screen.getByTestId('calculation-panel')).toBeInTheDocument()
      expect(screen.getByTestId('status-bar')).toBeInTheDocument()
    })

    test('should have correct initial state', () => {
      render(<App />)
      
      expect(screen.getByText('Molecule: None')).toBeInTheDocument()
      expect(screen.getByText('Results: None')).toBeInTheDocument()
      expect(screen.getByText('Status: Ready')).toBeInTheDocument()
      expect(screen.getByText('Current: None')).toBeInTheDocument()
    })

    test('should setup electron API event listeners', () => {
      render(<App />)
      
      expect(mockElectronAPI.onPythonLog).toHaveBeenCalled()
      expect(mockElectronAPI.onPythonError).toHaveBeenCalled()
      expect(mockElectronAPI.onCalculationProgress).toHaveBeenCalled()
    })
  })

  describe('Molecule Loading', () => {
    test('should load molecule successfully', async () => {
      mockElectronAPI.buildMolecule.mockResolvedValue({
        status: 'success',
        molecule: {
          natoms: 2,
          charge: 0,
          spin: 0,
          elements: ['H', 'H']
        }
      })

      render(<App />)
      
      fireEvent.click(screen.getByText('Load Test Molecule'))
      
      await waitFor(() => {
        expect(mockElectronAPI.buildMolecule).toHaveBeenCalledWith({
          type: 'coordinates',
          coordinates: [['H', 0, 0, 0], ['H', 0.74, 0, 0]],
          charge: 0,
          spin: 0
        })
      })

      await waitFor(() => {
        expect(screen.getByText('Molecule: 2 atoms')).toBeInTheDocument()
        expect(screen.getByText('Current: Loaded')).toBeInTheDocument()
      })
    })

    test('should handle molecule loading failure', async () => {
      mockElectronAPI.buildMolecule.mockResolvedValue({
        status: 'error',
        message: 'Invalid molecule data'
      })

      // Mock alert
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation()

      render(<App />)
      
      fireEvent.click(screen.getByText('Load Test Molecule'))
      
      await waitFor(() => {
        expect(alertSpy).toHaveBeenCalledWith('Failed to build molecule: Invalid molecule data')
      })

      expect(screen.getByText('Molecule: None')).toBeInTheDocument()
      
      alertSpy.mockRestore()
    })

    test('should handle molecule loading exception', async () => {
      mockElectronAPI.buildMolecule.mockRejectedValue(new Error('Network error'))

      const alertSpy = jest.spyOn(window, 'alert').mockImplementation()
      const consoleSpy = jest.spyOn(console, 'error').mockImplementation()

      render(<App />)
      
      fireEvent.click(screen.getByText('Load Test Molecule'))
      
      await waitFor(() => {
        expect(alertSpy).toHaveBeenCalledWith('Molecule building error: Error: Network error')
        expect(consoleSpy).toHaveBeenCalledWith('Molecule building error:', expect.any(Error))
      })

      alertSpy.mockRestore()
      consoleSpy.mockRestore()
    })
  })

  describe('Calculation Workflow', () => {
    beforeEach(async () => {
      // Setup molecule first
      mockElectronAPI.buildMolecule.mockResolvedValue({
        status: 'success',
        molecule: {
          natoms: 2,
          charge: 0,
          spin: 0,
          elements: ['H', 'H']
        }
      })

      render(<App />)
      fireEvent.click(screen.getByText('Load Test Molecule'))
      
      await waitFor(() => {
        expect(screen.getByText('Molecule: 2 atoms')).toBeInTheDocument()
      })
    })

    test('should run calculation successfully', async () => {
      mockElectronAPI.calculate.mockResolvedValue({
        status: 'success',
        results: {
          method: 'HF',
          energy: -1.123456,
          converged: true
        }
      })

      fireEvent.click(screen.getByText('Start Calculation'))
      
      await waitFor(() => {
        expect(mockElectronAPI.calculate).toHaveBeenCalledWith({
          molecule: {
            natoms: 2,
            charge: 0,
            spin: 0,
            elements: ['H', 'H']
          },
          method: 'HF',
          basis: 'STO-3G'
        })
      })

      await waitFor(() => {
        expect(screen.getByText('Results: HF')).toBeInTheDocument()
        expect(screen.getByText('Results: Available')).toBeInTheDocument()
      })
    })

    test('should handle calculation without molecule', async () => {
      // Reset to no molecule state
      const { rerender } = render(<App />)
      
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation()

      fireEvent.click(screen.getByText('Start Calculation'))
      
      expect(alertSpy).toHaveBeenCalledWith('Please select a molecule first')
      expect(mockElectronAPI.calculate).not.toHaveBeenCalled()
      
      alertSpy.mockRestore()
    })

    test('should handle calculation failure', async () => {
      mockElectronAPI.calculate.mockResolvedValue({
        status: 'error',
        message: 'Calculation failed to converge'
      })

      const alertSpy = jest.spyOn(window, 'alert').mockImplementation()

      fireEvent.click(screen.getByText('Start Calculation'))
      
      await waitFor(() => {
        expect(alertSpy).toHaveBeenCalledWith('Calculation failed: Calculation failed to converge')
      })

      expect(screen.getByText('Results: None')).toBeInTheDocument()
      
      alertSpy.mockRestore()
    })

    test('should handle calculation exception', async () => {
      mockElectronAPI.calculate.mockRejectedValue(new Error('IPC error'))

      const alertSpy = jest.spyOn(window, 'alert').mockImplementation()
      const consoleSpy = jest.spyOn(console, 'error').mockImplementation()

      fireEvent.click(screen.getByText('Start Calculation'))
      
      await waitFor(() => {
        expect(alertSpy).toHaveBeenCalledWith('Calculation error: Error: IPC error')
        expect(consoleSpy).toHaveBeenCalledWith('Calculation error:', expect.any(Error))
      })

      alertSpy.mockRestore()
      consoleSpy.mockRestore()
    })

    test('should manage calculation loading state', async () => {
      let resolveCalculation: (value: any) => void
      const calculationPromise = new Promise(resolve => {
        resolveCalculation = resolve
      })
      mockElectronAPI.calculate.mockReturnValue(calculationPromise)

      fireEvent.click(screen.getByText('Start Calculation'))
      
      // Should show calculating state
      await waitFor(() => {
        expect(screen.getByText('Calculating...')).toBeInTheDocument()
        expect(screen.getByText('Status: Calculating')).toBeInTheDocument()
      })

      // Resolve calculation
      resolveCalculation!({
        status: 'success',
        results: { method: 'HF', energy: -1.0 }
      })

      // Should return to ready state
      await waitFor(() => {
        expect(screen.getByText('Start Calculation')).toBeInTheDocument()
        expect(screen.getByText('Status: Ready')).toBeInTheDocument()
      })
    })
  })

  describe('Event Handling', () => {
    test('should handle Python log events', () => {
      render(<App />)
      
      const logCallback = mockElectronAPI.onPythonLog.mock.calls[0][0]
      logCallback({ type: 'stdout', data: 'Test log message' })
      
      expect(screen.getByText('Logs: 1 entries')).toBeInTheDocument()
    })

    test('should handle Python error events', () => {
      render(<App />)
      
      const errorCallback = mockElectronAPI.onPythonError.mock.calls[0][0]
      errorCallback('Test error message')
      
      expect(screen.getByText('Logs: 1 entries')).toBeInTheDocument()
    })

    test('should handle calculation progress events', () => {
      const consoleSpy = jest.spyOn(console, 'log').mockImplementation()
      
      render(<App />)
      
      const progressCallback = mockElectronAPI.onCalculationProgress.mock.calls[0][0]
      progressCallback({ progress: 0.5 })
      
      expect(consoleSpy).toHaveBeenCalledWith('Calculation progress:', { progress: 0.5 })
      
      consoleSpy.mockRestore()
    })

    test('should clear logs when requested', () => {
      render(<App />)
      
      // Add some logs
      const logCallback = mockElectronAPI.onPythonLog.mock.calls[0][0]
      logCallback({ type: 'stdout', data: 'Test log 1' })
      logCallback({ type: 'stdout', data: 'Test log 2' })
      
      expect(screen.getByText('Logs: 2 entries')).toBeInTheDocument()
      
      // Clear logs
      fireEvent.click(screen.getByText('Clear Logs'))
      
      expect(screen.getByText('Logs: 0 entries')).toBeInTheDocument()
    })
  })

  describe('Component Cleanup', () => {
    test('should cleanup event listeners on unmount', () => {
      const { unmount } = render(<App />)
      
      unmount()
      
      expect(mockElectronAPI.removeAllListeners).toHaveBeenCalledWith('python-log')
      expect(mockElectronAPI.removeAllListeners).toHaveBeenCalledWith('python-error')
      expect(mockElectronAPI.removeAllListeners).toHaveBeenCalledWith('calculation-progress')
    })
  })

  describe('State Persistence', () => {
    test('should maintain molecule state across calculations', async () => {
      // Load molecule
      mockElectronAPI.buildMolecule.mockResolvedValue({
        status: 'success',
        molecule: { natoms: 2, charge: 0, spin: 0 }
      })

      render(<App />)
      fireEvent.click(screen.getByText('Load Test Molecule'))
      
      await waitFor(() => {
        expect(screen.getByText('Molecule: 2 atoms')).toBeInTheDocument()
      })

      // Run first calculation
      mockElectronAPI.calculate.mockResolvedValue({
        status: 'success',
        results: { method: 'HF', energy: -1.0 }
      })

      fireEvent.click(screen.getByText('Start Calculation'))
      
      await waitFor(() => {
        expect(screen.getByText('Results: HF')).toBeInTheDocument()
      })

      // Molecule should still be loaded
      expect(screen.getByText('Molecule: 2 atoms')).toBeInTheDocument()
      expect(screen.getByText('Current: Loaded')).toBeInTheDocument()
    })

    test('should replace calculation results with new ones', async () => {
      // Setup molecule
      mockElectronAPI.buildMolecule.mockResolvedValue({
        status: 'success',
        molecule: { natoms: 2, charge: 0, spin: 0 }
      })

      render(<App />)
      fireEvent.click(screen.getByText('Load Test Molecule'))
      
      await waitFor(() => {
        expect(screen.getByText('Molecule: 2 atoms')).toBeInTheDocument()
      })

      // First calculation
      mockElectronAPI.calculate.mockResolvedValue({
        status: 'success',
        results: { method: 'HF', energy: -1.0 }
      })

      fireEvent.click(screen.getByText('Start Calculation'))
      
      await waitFor(() => {
        expect(screen.getByText('Results: HF')).toBeInTheDocument()
      })

      // Second calculation with different method
      mockElectronAPI.calculate.mockResolvedValue({
        status: 'success',
        results: { method: 'DFT', energy: -1.2 }
      })

      fireEvent.click(screen.getByText('Start Calculation'))
      
      await waitFor(() => {
        expect(screen.getByText('Results: DFT')).toBeInTheDocument()
      })
    })
  })
})