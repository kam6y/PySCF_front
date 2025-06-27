/**
 * Integration tests for React components
 * Tests component interactions and data flow
 */

import React from 'react'
import { render, screen, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import '@testing-library/jest-dom'
import App from '../../src/renderer/App'

// Enhanced mock for electronAPI
const mockElectronAPI = {
  calculate: jest.fn(),
  buildMolecule: jest.fn(),
  getCalculationStatus: jest.fn(),
  stopCalculation: jest.fn(),
  onPythonLog: jest.fn(),
  onPythonError: jest.fn(),
  onCalculationProgress: jest.fn(),
  removeAllListeners: jest.fn(),
}

// Override the global mock
Object.defineProperty(window, 'electronAPI', {
  value: mockElectronAPI,
  writable: true,
})

describe('App Integration Tests', () => {
  beforeEach(() => {
    jest.clearAllMocks()
    
    // Set up default mock returns
    mockElectronAPI.buildMolecule.mockResolvedValue({
      status: 'success',
      molecule: {
        natoms: 3,
        charge: 0,
        spin: 0,
        coordinates: [[0, 0, 0], [0.757, 0.586, 0], [-0.757, 0.586, 0]],
        elements: ['O', 'H', 'H']
      }
    })
    
    mockElectronAPI.calculate.mockResolvedValue({
      status: 'success',
      results: {
        method: 'HF',
        energy: -74.962947,
        converged: true,
        homo_energy: -0.5,
        lumo_energy: 0.2,
        homo_lumo_gap: 0.7,
        dipole_moment: [0.0, 1.7, 0.0],
        dipole_magnitude: 1.7
      }
    })
  })

  it('renders main application layout', () => {
    render(<App />)
    
    // Check main layout elements
    expect(screen.getByText('PySCF Front')).toBeInTheDocument()
    expect(screen.getByText('Quantum Chemistry Calculation')).toBeInTheDocument()
    expect(screen.getByText('3D Molecular Viewer')).toBeInTheDocument()
  })

  it('loads molecule and performs calculation workflow', async () => {
    const user = userEvent.setup()
    render(<App />)
    
    // Step 1: Load a molecule
    const waterMolecule = screen.getByText('Water (H₂O)')
    await user.click(waterMolecule)
    
    // Verify molecule building was called
    expect(mockElectronAPI.buildMolecule).toHaveBeenCalledWith({
      type: 'coordinates',
      coordinates: [
        ['O', 0.0, 0.0, 0.0],
        ['H', 0.757, 0.586, 0.0],
        ['H', -0.757, 0.586, 0.0]
      ],
      charge: 0,
      spin: 0
    })
    
    // Wait for molecule to be loaded
    await waitFor(() => {
      expect(screen.getByText('Molecule Loaded')).toBeInTheDocument()
    })
    
    // Step 2: Start calculation
    const calculateButton = screen.getByRole('button', { name: 'Start Calculation' })
    await user.click(calculateButton)
    
    // Verify calculation was called
    expect(mockElectronAPI.calculate).toHaveBeenCalledWith({
      molecule: {
        natoms: 3,
        charge: 0,
        spin: 0,
        coordinates: [[0, 0, 0], [0.757, 0.586, 0], [-0.757, 0.586, 0]],
        elements: ['O', 'H', 'H']
      },
      method: 'DFT',
      functional: 'B3LYP',
      basis: '6-31G*'
    })
    
    // Step 3: Verify results are displayed
    await waitFor(() => {
      expect(screen.getByText('Calculation Results')).toBeInTheDocument()
      expect(screen.getByText('-74.96294700 Ha')).toBeInTheDocument()
    })
  })

  it('handles calculation errors gracefully', async () => {
    const user = userEvent.setup()
    
    // Mock calculation failure
    mockElectronAPI.calculate.mockResolvedValue({
      status: 'error',
      message: 'Calculation failed'
    })
    
    // Mock window.alert
    const alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {})
    
    render(<App />)
    
    // Load molecule first
    const waterMolecule = screen.getByText('Water (H₂O)')
    await user.click(waterMolecule)
    
    await waitFor(() => {
      expect(screen.getByText('Molecule Loaded')).toBeInTheDocument()
    })
    
    // Try calculation
    const calculateButton = screen.getByRole('button', { name: 'Start Calculation' })
    await user.click(calculateButton)
    
    // Verify error handling
    await waitFor(() => {
      expect(alertSpy).toHaveBeenCalledWith('Calculation failed: Calculation failed')
    })
    
    alertSpy.mockRestore()
  })

  it('prevents calculation without molecule', async () => {
    const user = userEvent.setup()
    
    // Mock window.alert
    const alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {})
    
    render(<App />)
    
    // Try to calculate without loading molecule
    const calculateButton = screen.getByRole('button', { name: 'Start Calculation' })
    await user.click(calculateButton)
    
    // Verify prevention
    expect(alertSpy).toHaveBeenCalledWith('Please select a molecule first')
    expect(mockElectronAPI.calculate).not.toHaveBeenCalled()
    
    alertSpy.mockRestore()
  })

  it('switches between molecules correctly', async () => {
    const user = userEvent.setup()
    
    // Mock different molecule response
    mockElectronAPI.buildMolecule.mockResolvedValueOnce({
      status: 'success',
      molecule: {
        natoms: 2,
        charge: 0,
        spin: 0,
        coordinates: [[0, 0, 0], [0.74, 0, 0]],
        elements: ['H', 'H']
      }
    })
    
    render(<App />)
    
    // Load hydrogen molecule
    const hydrogenMolecule = screen.getByText('Hydrogen (H₂)')
    await user.click(hydrogenMolecule)
    
    // Verify correct molecule data was sent
    expect(mockElectronAPI.buildMolecule).toHaveBeenCalledWith({
      type: 'coordinates',
      coordinates: [
        ['H', 0.0, 0.0, 0.0],
        ['H', 0.74, 0.0, 0.0]
      ],
      charge: 0,
      spin: 0
    })
    
    await waitFor(() => {
      expect(screen.getByText('Atoms: 2')).toBeInTheDocument()
    })
  })

  it('handles molecule building errors', async () => {
    const user = userEvent.setup()
    
    // Mock molecule building failure
    mockElectronAPI.buildMolecule.mockResolvedValue({
      status: 'error',
      message: 'Failed to build molecule'
    })
    
    // Mock window.alert
    const alertSpy = jest.spyOn(window, 'alert').mockImplementation(() => {})
    
    render(<App />)
    
    // Try to load molecule
    const waterMolecule = screen.getByText('Water (H₂O)')
    await user.click(waterMolecule)
    
    // Verify error handling
    await waitFor(() => {
      expect(alertSpy).toHaveBeenCalledWith('Failed to build molecule: Failed to build molecule')
    })
    
    alertSpy.mockRestore()
  })

  it('updates calculation method and calls with correct parameters', async () => {
    const user = userEvent.setup()
    render(<App />)
    
    // Load molecule
    const waterMolecule = screen.getByText('Water (H₂O)')
    await user.click(waterMolecule)
    
    await waitFor(() => {
      expect(screen.getByText('Molecule Loaded')).toBeInTheDocument()
    })
    
    // Change method to HF
    const methodSelect = screen.getByLabelText('Calculation Method')
    await user.selectOptions(methodSelect, 'HF')
    
    // Change basis set
    const basisSelect = screen.getByLabelText('Basis Set')
    await user.selectOptions(basisSelect, 'STO-3G')
    
    // Start calculation
    const calculateButton = screen.getByRole('button', { name: 'Start Calculation' })
    await user.click(calculateButton)
    
    // Verify calculation called with correct parameters
    expect(mockElectronAPI.calculate).toHaveBeenCalledWith({
      molecule: expect.any(Object),
      method: 'HF',
      functional: undefined,
      basis: 'STO-3G'
    })
  })

  it('displays different calculation results correctly', async () => {
    const user = userEvent.setup()
    
    // Mock DFT result
    mockElectronAPI.calculate.mockResolvedValue({
      status: 'success',
      results: {
        method: 'DFT',
        functional: 'PBE',
        energy: -75.123456,
        converged: true,
        dipole_moment: [0.1, 1.8, 0.0],
        dipole_magnitude: 1.8
      }
    })
    
    render(<App />)
    
    // Load molecule and calculate
    const waterMolecule = screen.getByText('Water (H₂O)')
    await user.click(waterMolecule)
    
    await waitFor(() => {
      expect(screen.getByText('Molecule Loaded')).toBeInTheDocument()
    })
    
    const calculateButton = screen.getByRole('button', { name: 'Start Calculation' })
    await user.click(calculateButton)
    
    // Verify DFT-specific results are shown
    await waitFor(() => {
      expect(screen.getByText('DFT')).toBeInTheDocument()
      expect(screen.getByText('PBE')).toBeInTheDocument()
      expect(screen.getByText('-75.12345600 Ha')).toBeInTheDocument()
      expect(screen.getByText('1.8000 Debye')).toBeInTheDocument()
    })
  })

  it('sets up event listeners on mount', () => {
    render(<App />)
    
    // Verify event listeners were set up
    expect(mockElectronAPI.onPythonLog).toHaveBeenCalled()
    expect(mockElectronAPI.onPythonError).toHaveBeenCalled()
    expect(mockElectronAPI.onCalculationProgress).toHaveBeenCalled()
  })

  it('cleans up event listeners on unmount', () => {
    const { unmount } = render(<App />)
    
    unmount()
    
    // Verify cleanup
    expect(mockElectronAPI.removeAllListeners).toHaveBeenCalledWith('python-log')
    expect(mockElectronAPI.removeAllListeners).toHaveBeenCalledWith('python-error')
    expect(mockElectronAPI.removeAllListeners).toHaveBeenCalledWith('calculation-progress')
  })
})