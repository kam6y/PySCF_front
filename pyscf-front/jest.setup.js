// Import React Testing Library setup
import '@testing-library/jest-dom'

// Mock Electron API for testing
const mockElectronAPI = {
  calculate: jest.fn().mockImplementation(() => 
    Promise.resolve({
      status: 'success',
      results: {
        method: 'HF',
        energy: -1.123456,
        converged: true,
        num_electrons: 2
      }
    })
  ),
  buildMolecule: jest.fn().mockImplementation(() =>
    Promise.resolve({
      status: 'success',
      molecule: {
        natoms: 2,
        charge: 0,
        spin: 0,
        elements: ['H', 'H']
      }
    })
  ),
  getCalculationStatus: jest.fn().mockImplementation(() =>
    Promise.resolve({
      status: 'success',
      backend_status: 'running',
      jobs: []
    })
  ),
  stopCalculation: jest.fn().mockImplementation(() =>
    Promise.resolve({ status: 'success' })
  ),
  onPythonLog: jest.fn(),
  onPythonError: jest.fn(),
  onCalculationProgress: jest.fn(),
  removeAllListeners: jest.fn(),
}

global.electronAPI = mockElectronAPI

// Mock window.electronAPI
Object.defineProperty(window, 'electronAPI', {
  value: mockElectronAPI,
  writable: true,
  configurable: true
})

// Reset all mock functions before each test
beforeEach(() => {
  jest.clearAllMocks()
  
  // Reset to default implementations
  mockElectronAPI.calculate.mockImplementation(() =>
    Promise.resolve({
      status: 'success',
      results: {
        method: 'HF',
        energy: -1.123456,
        converged: true,
        num_electrons: 2
      }
    })
  )
  
  mockElectronAPI.buildMolecule.mockImplementation(() =>
    Promise.resolve({
      status: 'success',
      molecule: {
        natoms: 2,
        charge: 0,
        spin: 0,
        elements: ['H', 'H']
      }
    })
  )
})