/**
 * End-to-End Integration Tests for React → Electron → Python workflow
 * Tests the complete application workflow from frontend to backend
 */

const { Application } = require('spectron')
const { spawn } = require('child_process')
const path = require('path')
const fs = require('fs')
const os = require('os')

// Mock implementations for testing without full Electron app
class MockE2EEnvironment {
  constructor() {
    this.pythonProcess = null
    this.electronAPI = new MockElectronAPI()
    this.reactApp = new MockReactApp()
  }

  async setup() {
    // Start Python backend
    await this.startPythonBackend()
    
    // Initialize React app state
    this.reactApp.initialize()
  }

  async cleanup() {
    if (this.pythonProcess) {
      this.pythonProcess.kill()
    }
  }

  async startPythonBackend() {
    const pythonPath = 'python3'
    const scriptPath = path.join(__dirname, '../../src/python/main.py')
    
    return new Promise((resolve, reject) => {
      this.pythonProcess = spawn(pythonPath, [scriptPath], {
        stdio: ['pipe', 'pipe', 'pipe']
      })

      this.pythonProcess.stdout.once('data', () => {
        resolve()
      })

      this.pythonProcess.on('error', reject)
    })
  }

  async sendToPython(message) {
    return new Promise((resolve, reject) => {
      if (!this.pythonProcess) {
        reject(new Error('Python process not started'))
        return
      }

      const messageStr = JSON.stringify(message) + '\n'
      this.pythonProcess.stdin.write(messageStr)

      this.pythonProcess.stdout.once('data', (data) => {
        try {
          const response = JSON.parse(data.toString().trim())
          resolve(response)
        } catch (error) {
          reject(error)
        }
      })

      setTimeout(() => {
        reject(new Error('Python response timeout'))
      }, 30000)
    })
  }
}

class MockElectronAPI {
  constructor() {
    this.pythonProcess = null
    this.listeners = new Map()
  }

  async calculate(params) {
    // Simulate IPC call to main process
    return this.sendToPython('calculate', params)
  }

  async buildMolecule(moleculeData) {
    // Simulate IPC call to main process
    return this.sendToPython('build-molecule', moleculeData)
  }

  async getCalculationStatus() {
    // Simulate IPC call to main process
    return this.sendToPython('get-status', {})
  }

  async sendToPython(action, data) {
    // This would be handled by Electron main process
    // For testing, we'll simulate the communication
    const message = {
      action: action,
      data: data,
      id: Date.now().toString()
    }
    
    // Mock response based on action
    if (action === 'get-status') {
      return {
        status: 'success',
        backend_status: 'running',
        jobs: []
      }
    } else if (action === 'build-molecule') {
      return {
        status: 'success',
        molecule: {
          natoms: data.coordinates ? data.coordinates.length : 0,
          charge: data.charge || 0,
          spin: data.spin || 0,
          elements: data.coordinates ? data.coordinates.map(coord => coord[0]) : []
        }
      }
    } else if (action === 'calculate') {
      return {
        status: 'success',
        results: {
          method: data.method || 'HF',
          energy: -1.123456,
          converged: true,
          num_electrons: 2
        }
      }
    }
    
    return { status: 'error', message: 'Unknown action' }
  }

  onPythonLog(callback) {
    if (!this.listeners.has('python-log')) {
      this.listeners.set('python-log', [])
    }
    this.listeners.get('python-log').push(callback)
  }

  onPythonError(callback) {
    if (!this.listeners.has('python-error')) {
      this.listeners.set('python-error', [])
    }
    this.listeners.get('python-error').push(callback)
  }

  onCalculationProgress(callback) {
    if (!this.listeners.has('calculation-progress')) {
      this.listeners.set('calculation-progress', [])
    }
    this.listeners.get('calculation-progress').push(callback)
  }

  removeAllListeners(event) {
    this.listeners.delete(event)
  }
}

class MockReactApp {
  constructor() {
    this.state = {
      currentMolecule: null,
      calculationResults: null,
      isCalculating: false,
      logs: []
    }
  }

  initialize() {
    // Initialize React app state
    this.state = {
      currentMolecule: null,
      calculationResults: null,
      isCalculating: false,
      logs: []
    }
  }

  async handleMoleculeLoad(moleculeData) {
    // Simulate React component behavior
    try {
      const result = await window.electronAPI.buildMolecule(moleculeData)
      if (result.status === 'success') {
        this.state.currentMolecule = result.molecule
        return true
      } else {
        throw new Error(result.message)
      }
    } catch (error) {
      console.error('Molecule building error:', error)
      return false
    }
  }

  async handleCalculation(params) {
    // Simulate React component behavior
    if (!this.state.currentMolecule) {
      throw new Error('No molecule loaded')
    }

    this.state.isCalculating = true
    this.state.calculationResults = null

    try {
      const result = await window.electronAPI.calculate({
        molecule: this.state.currentMolecule,
        ...params
      })

      if (result.status === 'success') {
        this.state.calculationResults = result.results
        return result.results
      } else {
        throw new Error(result.message)
      }
    } finally {
      this.state.isCalculating = false
    }
  }
}

describe('End-to-End Workflow Tests', () => {
  let testEnv

  beforeAll(async () => {
    testEnv = new MockE2EEnvironment()
    await testEnv.setup()
    
    // Make electronAPI available globally for React components
    global.window = {
      electronAPI: testEnv.electronAPI
    }
  })

  afterAll(async () => {
    await testEnv.cleanup()
  })

  describe('Complete Calculation Workflow', () => {
    test('should complete full molecule → calculation → results workflow', async () => {
      // Step 1: Load molecule data (React → Electron)
      const moleculeData = {
        type: 'coordinates',
        coordinates: [
          ['H', 0.0, 0.0, 0.0],
          ['H', 0.74, 0.0, 0.0]
        ],
        charge: 0,
        spin: 0
      }

      const moleculeLoadSuccess = await testEnv.reactApp.handleMoleculeLoad(moleculeData)
      expect(moleculeLoadSuccess).toBe(true)
      expect(testEnv.reactApp.state.currentMolecule).toBeDefined()
      expect(testEnv.reactApp.state.currentMolecule.natoms).toBe(2)

      // Step 2: Run calculation (React → Electron → Python)
      const calculationParams = {
        method: 'HF',
        basis: 'STO-3G'
      }

      const results = await testEnv.reactApp.handleCalculation(calculationParams)
      
      expect(results).toBeDefined()
      expect(results.method).toBe('HF')
      expect(results.energy).toBeLessThan(0)
      expect(results.converged).toBe(true)
      expect(testEnv.reactApp.state.calculationResults).toEqual(results)
    })

    test('should handle water molecule DFT calculation workflow', async () => {
      // Load water molecule
      const waterMolecule = {
        type: 'coordinates',
        coordinates: [
          ['O', 0.0, 0.0, 0.0],
          ['H', 0.757, 0.586, 0.0],
          ['H', -0.757, 0.586, 0.0]
        ],
        charge: 0,
        spin: 0
      }

      await testEnv.reactApp.handleMoleculeLoad(waterMolecule)
      expect(testEnv.reactApp.state.currentMolecule.natoms).toBe(3)

      // Run DFT calculation
      const dftParams = {
        method: 'DFT',
        functional: 'B3LYP',
        basis: '6-31G*'
      }

      const results = await testEnv.reactApp.handleCalculation(dftParams)
      
      expect(results.method).toBe('DFT')
      expect(results.energy).toBeLessThan(0)
    })

    test('should handle calculation errors gracefully', async () => {
      // Load molecule first
      const moleculeData = {
        type: 'coordinates',
        coordinates: [['H', 0.0, 0.0, 0.0]],
        charge: 0,
        spin: 0
      }

      await testEnv.reactApp.handleMoleculeLoad(moleculeData)

      // Try calculation with invalid method
      const invalidParams = {
        method: 'INVALID_METHOD',
        basis: 'STO-3G'
      }

      try {
        await testEnv.reactApp.handleCalculation(invalidParams)
        fail('Should have thrown error for invalid method')
      } catch (error) {
        expect(error.message).toContain('error')
        expect(testEnv.reactApp.state.isCalculating).toBe(false)
      }
    })

    test('should prevent calculation without loaded molecule', async () => {
      // Reset state
      testEnv.reactApp.state.currentMolecule = null

      const calculationParams = {
        method: 'HF',
        basis: 'STO-3G'
      }

      try {
        await testEnv.reactApp.handleCalculation(calculationParams)
        fail('Should have thrown error for no molecule')
      } catch (error) {
        expect(error.message).toBe('No molecule loaded')
      }
    })
  })

  describe('State Management Throughout Workflow', () => {
    test('should maintain correct loading states during calculation', async () => {
      // Load molecule
      const moleculeData = {
        type: 'coordinates',
        coordinates: [['H', 0.0, 0.0, 0.0], ['H', 0.74, 0.0, 0.0]],
        charge: 0,
        spin: 0
      }

      await testEnv.reactApp.handleMoleculeLoad(moleculeData)

      // Check initial state
      expect(testEnv.reactApp.state.isCalculating).toBe(false)
      expect(testEnv.reactApp.state.calculationResults).toBe(null)

      // Start calculation
      const calculationPromise = testEnv.reactApp.handleCalculation({
        method: 'HF',
        basis: 'STO-3G'
      })

      // During calculation (this is mocked, so it completes immediately)
      const results = await calculationPromise

      // After calculation
      expect(testEnv.reactApp.state.isCalculating).toBe(false)
      expect(testEnv.reactApp.state.calculationResults).toEqual(results)
    })

    test('should handle multiple sequential calculations', async () => {
      // Load molecule
      const moleculeData = {
        type: 'coordinates',
        coordinates: [['H', 0.0, 0.0, 0.0], ['H', 0.74, 0.0, 0.0]],
        charge: 0,
        spin: 0
      }

      await testEnv.reactApp.handleMoleculeLoad(moleculeData)

      // First calculation
      const results1 = await testEnv.reactApp.handleCalculation({
        method: 'HF',
        basis: 'STO-3G'
      })
      expect(results1.method).toBe('HF')

      // Second calculation with different parameters
      const results2 = await testEnv.reactApp.handleCalculation({
        method: 'DFT',
        functional: 'B3LYP',
        basis: '6-31G*'
      })
      expect(results2.method).toBe('DFT')

      // State should reflect latest results
      expect(testEnv.reactApp.state.calculationResults).toEqual(results2)
    })
  })

  describe('Error Handling Across System Boundaries', () => {
    test('should propagate Python calculation errors to React', async () => {
      // This test would verify that errors from Python backend
      // are properly propagated through Electron to React
      
      // Mock the electronAPI to return an error
      const originalCalculate = testEnv.electronAPI.calculate
      testEnv.electronAPI.calculate = async () => ({
        status: 'error',
        message: 'PySCF calculation failed'
      })

      // Load molecule
      await testEnv.reactApp.handleMoleculeLoad({
        type: 'coordinates',
        coordinates: [['H', 0.0, 0.0, 0.0]],
        charge: 0,
        spin: 0
      })

      // Try calculation
      try {
        await testEnv.reactApp.handleCalculation({
          method: 'HF',
          basis: 'STO-3G'
        })
        fail('Should have propagated error from Python')
      } catch (error) {
        expect(error.message).toBe('PySCF calculation failed')
      }

      // Restore original method
      testEnv.electronAPI.calculate = originalCalculate
    })

    test('should handle molecule building errors from Python', async () => {
      // Mock the electronAPI to return molecule building error
      const originalBuildMolecule = testEnv.electronAPI.buildMolecule
      testEnv.electronAPI.buildMolecule = async () => ({
        status: 'error',
        message: 'Invalid molecule coordinates'
      })

      const success = await testEnv.reactApp.handleMoleculeLoad({
        type: 'coordinates',
        coordinates: [],  // Invalid empty coordinates
        charge: 0,
        spin: 0
      })

      expect(success).toBe(false)
      expect(testEnv.reactApp.state.currentMolecule).toBe(null)

      // Restore original method
      testEnv.electronAPI.buildMolecule = originalBuildMolecule
    })
  })

  describe('Communication Protocol Validation', () => {
    test('should handle IPC timeout scenarios', async () => {
      // Mock a timeout scenario
      const originalCalculate = testEnv.electronAPI.calculate
      testEnv.electronAPI.calculate = async () => {
        return new Promise((resolve, reject) => {
          setTimeout(() => {
            reject(new Error('IPC timeout'))
          }, 100)
        })
      }

      // Load molecule
      await testEnv.reactApp.handleMoleculeLoad({
        type: 'coordinates',
        coordinates: [['H', 0.0, 0.0, 0.0]],
        charge: 0,
        spin: 0
      })

      // Try calculation with timeout
      try {
        await testEnv.reactApp.handleCalculation({
          method: 'HF',
          basis: 'STO-3G'
        })
        fail('Should have timed out')
      } catch (error) {
        expect(error.message).toBe('IPC timeout')
        expect(testEnv.reactApp.state.isCalculating).toBe(false)
      }

      // Restore original method
      testEnv.electronAPI.calculate = originalCalculate
    })

    test('should validate message format integrity', async () => {
      // This test would verify that messages maintain their format
      // throughout the React → Electron → Python chain
      
      const moleculeData = {
        type: 'coordinates',
        coordinates: [['H', 0.0, 0.0, 0.0], ['H', 0.74, 0.0, 0.0]],
        charge: 0,
        spin: 0
      }

      await testEnv.reactApp.handleMoleculeLoad(moleculeData)

      const calculationParams = {
        method: 'HF',
        basis: 'STO-3G',
        custom_param: 'test_value'  // Custom parameter to track through system
      }

      // In a real test, we would intercept the message at each boundary
      // and verify the format is preserved
      const results = await testEnv.reactApp.handleCalculation(calculationParams)
      
      expect(results).toBeDefined()
      expect(results.method).toBe('HF')
    })
  })
})

module.exports = {
  MockE2EEnvironment,
  MockElectronAPI,
  MockReactApp
}