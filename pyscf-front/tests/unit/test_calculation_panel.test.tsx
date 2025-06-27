/**
 * Unit tests for CalculationPanel component
 */

import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import '@testing-library/jest-dom'
import CalculationPanel from '../../src/renderer/components/CalculationPanel'

// Mock props
const mockProps = {
  onCalculate: jest.fn(),
  isCalculating: false,
  results: null,
}

describe('CalculationPanel', () => {
  beforeEach(() => {
    jest.clearAllMocks()
  })

  it('renders calculation panel with basic elements', () => {
    render(<CalculationPanel {...mockProps} />)
    
    expect(screen.getByText('Quantum Chemistry Calculation')).toBeInTheDocument()
    expect(screen.getByText('Configure and run calculations')).toBeInTheDocument()
    expect(screen.getByLabelText('Calculation Method')).toBeInTheDocument()
    expect(screen.getByLabelText('Basis Set')).toBeInTheDocument()
    expect(screen.getByRole('button', { name: 'Start Calculation' })).toBeInTheDocument()
  })

  it('shows functional selector for DFT method', async () => {
    const user = userEvent.setup()
    render(<CalculationPanel {...mockProps} />)
    
    // DFT should be available and functional selector should be visible
    const methodSelect = screen.getByLabelText('Calculation Method')
    await user.selectOptions(methodSelect, 'DFT')
    
    expect(screen.getByLabelText('Exchange-Correlation Functional')).toBeInTheDocument()
  })

  it('hides functional selector for non-DFT methods', async () => {
    const user = userEvent.setup()
    render(<CalculationPanel {...mockProps} />)
    
    const methodSelect = screen.getByLabelText('Calculation Method')
    await user.selectOptions(methodSelect, 'HF')
    
    expect(screen.queryByLabelText('Exchange-Correlation Functional')).not.toBeInTheDocument()
  })

  it('calls onCalculate with correct parameters when submitted', async () => {
    const user = userEvent.setup()
    render(<CalculationPanel {...mockProps} />)
    
    const methodSelect = screen.getByLabelText('Calculation Method')
    const basisSelect = screen.getByLabelText('Basis Set')
    const submitButton = screen.getByRole('button', { name: 'Start Calculation' })
    
    await user.selectOptions(methodSelect, 'HF')
    await user.selectOptions(basisSelect, '6-31G*')
    await user.click(submitButton)
    
    expect(mockProps.onCalculate).toHaveBeenCalledWith({
      method: 'HF',
      functional: undefined,
      basis: '6-31G*'
    })
  })

  it('calls onCalculate with functional for DFT calculations', async () => {
    const user = userEvent.setup()
    render(<CalculationPanel {...mockProps} />)
    
    const methodSelect = screen.getByLabelText('Calculation Method')
    await user.selectOptions(methodSelect, 'DFT')
    
    const functionalSelect = screen.getByLabelText('Exchange-Correlation Functional')
    const submitButton = screen.getByRole('button', { name: 'Start Calculation' })
    
    await user.selectOptions(functionalSelect, 'PBE')
    await user.click(submitButton)
    
    expect(mockProps.onCalculate).toHaveBeenCalledWith({
      method: 'DFT',
      functional: 'PBE',
      basis: '6-31G*'
    })
  })

  it('disables form elements when calculating', () => {
    const calculatingProps = { ...mockProps, isCalculating: true }
    render(<CalculationPanel {...calculatingProps} />)
    
    expect(screen.getByLabelText('Calculation Method')).toBeDisabled()
    expect(screen.getByLabelText('Basis Set')).toBeDisabled()
    expect(screen.getByRole('button', { name: /calculating/i })).toBeDisabled()
  })

  it('shows calculating state in button', () => {
    const calculatingProps = { ...mockProps, isCalculating: true }
    render(<CalculationPanel {...calculatingProps} />)
    
    expect(screen.getByRole('button', { name: /calculating/i })).toBeInTheDocument()
    expect(screen.getByText('Calculating...')).toBeInTheDocument()
  })

  it('displays calculation results when provided', () => {
    const resultsProps = {
      ...mockProps,
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
    }
    
    render(<CalculationPanel {...resultsProps} />)
    
    expect(screen.getByText('Calculation Results')).toBeInTheDocument()
    expect(screen.getByText('Energy Information')).toBeInTheDocument()
    expect(screen.getByText('HF')).toBeInTheDocument()
    expect(screen.getByText('-74.96294700 Ha')).toBeInTheDocument()
    expect(screen.getByText('Yes')).toBeInTheDocument()
    
    expect(screen.getByText('Orbital Information')).toBeInTheDocument()
    expect(screen.getByText('-0.5000 Ha')).toBeInTheDocument()
    expect(screen.getByText('0.2000 Ha')).toBeInTheDocument()
    expect(screen.getByText('0.7000 Ha')).toBeInTheDocument()
    
    expect(screen.getByText('Molecular Properties')).toBeInTheDocument()
    expect(screen.getByText('1.7000 Debye')).toBeInTheDocument()
  })

  it('shows DFT functional in results', () => {
    const dftResultsProps = {
      ...mockProps,
      results: {
        method: 'DFT',
        functional: 'B3LYP',
        energy: -75.123456,
        converged: true
      }
    }
    
    render(<CalculationPanel {...dftResultsProps} />)
    
    expect(screen.getByText('B3LYP')).toBeInTheDocument()
  })

  it('shows convergence failure status', () => {
    const failedResultsProps = {
      ...mockProps,
      results: {
        method: 'HF',
        energy: -74.962947,
        converged: false
      }
    }
    
    render(<CalculationPanel {...failedResultsProps} />)
    
    expect(screen.getByText('No')).toBeInTheDocument()
  })

  it('has all required method options', () => {
    render(<CalculationPanel {...mockProps} />)
    
    const methodSelect = screen.getByLabelText('Calculation Method')
    const options = Array.from(methodSelect.querySelectorAll('option')).map(
      option => option.textContent
    )
    
    expect(options).toContain('HF')
    expect(options).toContain('DFT')
    expect(options).toContain('MP2')
  })

  it('has all required functional options', async () => {
    const user = userEvent.setup()
    render(<CalculationPanel {...mockProps} />)
    
    const methodSelect = screen.getByLabelText('Calculation Method')
    await user.selectOptions(methodSelect, 'DFT')
    
    const functionalSelect = screen.getByLabelText('Exchange-Correlation Functional')
    const options = Array.from(functionalSelect.querySelectorAll('option')).map(
      option => option.textContent
    )
    
    expect(options).toContain('B3LYP')
    expect(options).toContain('PBE')
    expect(options).toContain('PBE0')
    expect(options).toContain('M06-2X')
    expect(options).toContain('wB97X-D')
  })

  it('has all required basis set options', () => {
    render(<CalculationPanel {...mockProps} />)
    
    const basisSelect = screen.getByLabelText('Basis Set')
    const options = Array.from(basisSelect.querySelectorAll('option')).map(
      option => option.textContent
    )
    
    expect(options).toContain('STO-3G')
    expect(options).toContain('6-31G')
    expect(options).toContain('6-31G*')
    expect(options).toContain('6-31G**')
    expect(options).toContain('6-311G*')
    expect(options).toContain('cc-pVDZ')
    expect(options).toContain('cc-pVTZ')
  })
})