/**
 * Tests for MoleculeViewer component
 * Tests 3D visualization placeholder functionality and molecule display
 */

import React from 'react'
import { render, screen } from '@testing-library/react'
import '@testing-library/jest-dom'
import MoleculeViewer from '../../src/renderer/components/MoleculeViewer'

describe('MoleculeViewer Component', () => {
  const mockMolecule = {
    natoms: 3,
    charge: 0,
    spin: 0,
    elements: ['O', 'H', 'H'],
    coordinates: [
      [0.0, 0.0, 0.0],
      [0.757, 0.586, 0.0],
      [-0.757, 0.586, 0.0]
    ]
  }

  const mockCalculationResults = {
    method: 'DFT',
    functional: 'B3LYP',
    energy: -76.123456,
    converged: true,
    num_electrons: 10
  }

  describe('Molecule Display', () => {
    test('should display placeholder when no molecule is loaded', () => {
      render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      expect(screen.getByText('No Molecule Loaded')).toBeInTheDocument()
      expect(screen.getByText('Select or create a molecule to begin')).toBeInTheDocument()
      expect(screen.getByText('3Dpymol integration coming soon...')).toBeInTheDocument()
      expect(screen.getByText('⚛️')).toBeInTheDocument()
    })

    test('should display molecule information when molecule is loaded', () => {
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Molecule Loaded')).toBeInTheDocument()
      expect(screen.getByText('Atoms: 3')).toBeInTheDocument()
      expect(screen.getByText('Charge: 0')).toBeInTheDocument()
      expect(screen.getByText('Spin: 0')).toBeInTheDocument()
      expect(screen.getByText('🧬')).toBeInTheDocument()
    })

    test('should handle molecule with unknown atom count', () => {
      const incompleteeMolecule = {
        charge: 0,
        spin: 0
      }
      
      render(<MoleculeViewer molecule={incompleteeMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Atoms: Unknown')).toBeInTheDocument()
    })

    test('should display default values for missing charge and spin', () => {
      const incompleteMolecule = {
        natoms: 2
      }
      
      render(<MoleculeViewer molecule={incompleteMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Charge: 0')).toBeInTheDocument()
      expect(screen.getByText('Spin: 0')).toBeInTheDocument()
    })
  })

  describe('Calculation Results Display', () => {
    test('should not display results section when no results available', () => {
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.queryByText('Latest Results')).not.toBeInTheDocument()
    })

    test('should display calculation results when available', () => {
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={mockCalculationResults} />)
      
      expect(screen.getByText('Latest Results')).toBeInTheDocument()
      expect(screen.getByText('Method: DFT')).toBeInTheDocument()
      expect(screen.getByText('Energy: -76.123456 Hartree')).toBeInTheDocument()
    })

    test('should handle calculation results without energy', () => {
      const resultsWithoutEnergy = {
        method: 'HF',
        converged: true
      }
      
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={resultsWithoutEnergy} />)
      
      expect(screen.getByText('Method: HF')).toBeInTheDocument()
      // Should handle undefined energy gracefully
      expect(screen.getByText(/Energy:/)).toBeInTheDocument()
    })

    test('should format energy with 6 decimal places', () => {
      const resultsWithPreciseEnergy = {
        method: 'MP2',
        energy: -1.23456789
      }
      
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={resultsWithPreciseEnergy} />)
      
      expect(screen.getByText('Energy: -1.234568 Hartree')).toBeInTheDocument()
    })
  })

  describe('Viewer Controls', () => {
    test('should display viewer control buttons', () => {
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Reset View')).toBeInTheDocument()
      expect(screen.getByText('Toggle Labels')).toBeInTheDocument()
    })

    test('should display viewer info label', () => {
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.getByText('3D Molecular Viewer')).toBeInTheDocument()
    })

    test('control buttons should be present even without molecule', () => {
      render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      expect(screen.getByText('Reset View')).toBeInTheDocument()
      expect(screen.getByText('Toggle Labels')).toBeInTheDocument()
    })
  })

  describe('Component Structure', () => {
    test('should have correct CSS classes for styling', () => {
      const { container } = render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      const mainDiv = container.firstChild
      expect(mainDiv).toHaveClass('w-full', 'h-full', 'relative')
    })

    test('should apply gradient background', () => {
      const { container } = render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      const viewerDiv = container.querySelector('[style*="background"]')
      expect(viewerDiv).toHaveStyle({
        background: 'linear-gradient(135deg, #1e293b 0%, #0f172a 100%)'
      })
    })

    test('should position controls in top-right corner', () => {
      const { container } = render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      const controlsDiv = container.querySelector('.absolute.top-4.right-4')
      expect(controlsDiv).toBeInTheDocument()
    })

    test('should position info label in bottom-left corner', () => {
      const { container } = render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      const infoDiv = container.querySelector('.absolute.bottom-4.left-4')
      expect(infoDiv).toBeInTheDocument()
    })
  })

  describe('Component Updates', () => {
    test('should update display when molecule changes', () => {
      const { rerender } = render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      expect(screen.getByText('No Molecule Loaded')).toBeInTheDocument()
      
      rerender(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Molecule Loaded')).toBeInTheDocument()
      expect(screen.getByText('Atoms: 3')).toBeInTheDocument()
    })

    test('should update results display when calculation results change', () => {
      const { rerender } = render(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.queryByText('Latest Results')).not.toBeInTheDocument()
      
      rerender(<MoleculeViewer molecule={mockMolecule} calculationResults={mockCalculationResults} />)
      
      expect(screen.getByText('Latest Results')).toBeInTheDocument()
      expect(screen.getByText('Method: DFT')).toBeInTheDocument()
    })

    test('should clear results when set to null', () => {
      const { rerender } = render(<MoleculeViewer molecule={mockMolecule} calculationResults={mockCalculationResults} />)
      
      expect(screen.getByText('Latest Results')).toBeInTheDocument()
      
      rerender(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(screen.queryByText('Latest Results')).not.toBeInTheDocument()
    })
  })

  describe('Console Logging (Development Features)', () => {
    test('should log molecule loading', () => {
      const consoleSpy = jest.spyOn(console, 'log').mockImplementation()
      
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={null} />)
      
      expect(consoleSpy).toHaveBeenCalledWith('Loading molecule:', mockMolecule)
      
      consoleSpy.mockRestore()
    })

    test('should log calculation results display', () => {
      const consoleSpy = jest.spyOn(console, 'log').mockImplementation()
      
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={mockCalculationResults} />)
      
      expect(consoleSpy).toHaveBeenCalledWith('Displaying calculation results:', mockCalculationResults)
      
      consoleSpy.mockRestore()
    })

    test('should log viewer initialization', () => {
      const consoleSpy = jest.spyOn(console, 'log').mockImplementation()
      
      render(<MoleculeViewer molecule={null} calculationResults={null} />)
      
      expect(consoleSpy).toHaveBeenCalledWith('Molecule viewer initialized')
      
      consoleSpy.mockRestore()
    })
  })

  describe('Edge Cases', () => {
    test('should handle molecule with zero atoms', () => {
      const emptyMolecule = {
        natoms: 0,
        charge: 0,
        spin: 0
      }
      
      render(<MoleculeViewer molecule={emptyMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Atoms: 0')).toBeInTheDocument()
    })

    test('should handle negative charge and spin values', () => {
      const chargedMolecule = {
        natoms: 2,
        charge: -1,
        spin: 1
      }
      
      render(<MoleculeViewer molecule={chargedMolecule} calculationResults={null} />)
      
      expect(screen.getByText('Charge: -1')).toBeInTheDocument()
      expect(screen.getByText('Spin: 1')).toBeInTheDocument()
    })

    test('should handle very large energy values', () => {
      const largeEnergyResults = {
        method: 'HF',
        energy: -12345.6789012345
      }
      
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={largeEnergyResults} />)
      
      expect(screen.getByText('Energy: -12345.678901 Hartree')).toBeInTheDocument()
    })

    test('should handle zero energy', () => {
      const zeroEnergyResults = {
        method: 'HF',
        energy: 0.0
      }
      
      render(<MoleculeViewer molecule={mockMolecule} calculationResults={zeroEnergyResults} />)
      
      expect(screen.getByText('Energy: 0.000000 Hartree')).toBeInTheDocument()
    })
  })
})