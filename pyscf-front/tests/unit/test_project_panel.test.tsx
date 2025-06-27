/**
 * Unit tests for ProjectPanel component
 */

import React from 'react'
import { render, screen, fireEvent } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import '@testing-library/jest-dom'
import ProjectPanel from '../../src/renderer/components/ProjectPanel'

// Mock props
const mockProps = {
  onMoleculeLoad: jest.fn(),
  currentMolecule: null,
}

describe('ProjectPanel', () => {
  beforeEach(() => {
    jest.clearAllMocks()
  })

  it('renders project panel with basic elements', () => {
    render(<ProjectPanel {...mockProps} />)
    
    expect(screen.getByText('PySCF Front')).toBeInTheDocument()
    expect(screen.getByText('Quantum Chemistry GUI')).toBeInTheDocument()
    expect(screen.getByRole('button', { name: 'Molecules' })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: 'Projects' })).toBeInTheDocument()
  })

  it('shows molecules tab by default', () => {
    render(<ProjectPanel {...mockProps} />)
    
    expect(screen.getByText('Sample Molecules')).toBeInTheDocument()
    expect(screen.getByText('Select a molecule to load into the viewer')).toBeInTheDocument()
  })

  it('displays sample molecules', () => {
    render(<ProjectPanel {...mockProps} />)
    
    expect(screen.getByText('Water (H₂O)')).toBeInTheDocument()
    expect(screen.getByText('H₂O')).toBeInTheDocument()
    expect(screen.getByText('Methane (CH₄)')).toBeInTheDocument()
    expect(screen.getByText('CH₄')).toBeInTheDocument()
    expect(screen.getByText('Hydrogen (H₂)')).toBeInTheDocument()
    expect(screen.getByText('H₂')).toBeInTheDocument()
  })

  it('calls onMoleculeLoad when molecule is clicked', async () => {
    const user = userEvent.setup()
    render(<ProjectPanel {...mockProps} />)
    
    const waterMolecule = screen.getByText('Water (H₂O)')
    await user.click(waterMolecule)
    
    expect(mockProps.onMoleculeLoad).toHaveBeenCalledWith({
      type: 'coordinates',
      coordinates: [
        ['O', 0.0, 0.0, 0.0],
        ['H', 0.757, 0.586, 0.0],
        ['H', -0.757, 0.586, 0.0]
      ],
      charge: 0,
      spin: 0
    })
  })

  it('highlights selected molecule', () => {
    const propsWithMolecule = {
      ...mockProps,
      currentMolecule: {
        coordinates: [
          ['O', 0.0, 0.0, 0.0],
          ['H', 0.757, 0.586, 0.0],
          ['H', -0.757, 0.586, 0.0]
        ]
      }
    }
    
    render(<ProjectPanel {...propsWithMolecule} />)
    
    const waterMoleculeContainer = screen.getByText('Water (H₂O)').closest('div')
    expect(waterMoleculeContainer).toHaveClass('border-primary')
  })

  it('switches to projects tab', async () => {
    const user = userEvent.setup()
    render(<ProjectPanel {...mockProps} />)
    
    const projectsTab = screen.getByRole('button', { name: 'Projects' })
    await user.click(projectsTab)
    
    expect(screen.getByText('Recent Projects')).toBeInTheDocument()
    expect(screen.getByText('Manage your calculation projects')).toBeInTheDocument()
    expect(screen.getByText('No projects yet')).toBeInTheDocument()
  })

  it('shows input method buttons', () => {
    render(<ProjectPanel {...mockProps} />)
    
    expect(screen.getByText('Input Methods')).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /load xyz file/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /draw molecule/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /smiles input/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /pubchem search/i })).toBeInTheDocument()
  })

  it('shows project management buttons in projects tab', async () => {
    const user = userEvent.setup()
    render(<ProjectPanel {...mockProps} />)
    
    const projectsTab = screen.getByRole('button', { name: 'Projects' })
    await user.click(projectsTab)
    
    expect(screen.getByRole('button', { name: /new project/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /open project/i })).toBeInTheDocument()
  })

  it('shows active tab styling', () => {
    render(<ProjectPanel {...mockProps} />)
    
    const moleculesTab = screen.getByRole('button', { name: 'Molecules' })
    const projectsTab = screen.getByRole('button', { name: 'Projects' })
    
    expect(moleculesTab).toHaveClass('border-primary', 'text-primary')
    expect(projectsTab).not.toHaveClass('border-primary', 'text-primary')
  })

  it('changes active tab styling when clicked', async () => {
    const user = userEvent.setup()
    render(<ProjectPanel {...mockProps} />)
    
    const moleculesTab = screen.getByRole('button', { name: 'Molecules' })
    const projectsTab = screen.getByRole('button', { name: 'Projects' })
    
    await user.click(projectsTab)
    
    expect(projectsTab).toHaveClass('border-primary', 'text-primary')
    expect(moleculesTab).not.toHaveClass('border-primary', 'text-primary')
  })

  it('calls onMoleculeLoad for methane molecule', async () => {
    const user = userEvent.setup()
    render(<ProjectPanel {...mockProps} />)
    
    const methaneMolecule = screen.getByText('Methane (CH₄)')
    await user.click(methaneMolecule)
    
    expect(mockProps.onMoleculeLoad).toHaveBeenCalledWith({
      type: 'coordinates',
      coordinates: [
        ['C', 0.0, 0.0, 0.0],
        ['H', 1.089, 0.0, 0.0],
        ['H', -0.363, 1.027, 0.0],
        ['H', -0.363, -0.513, 0.889],
        ['H', -0.363, -0.513, -0.889]
      ],
      charge: 0,
      spin: 0
    })
  })

  it('calls onMoleculeLoad for hydrogen molecule', async () => {
    const user = userEvent.setup()
    render(<ProjectPanel {...mockProps} />)
    
    const hydrogenMolecule = screen.getByText('Hydrogen (H₂)')
    await user.click(hydrogenMolecule)
    
    expect(mockProps.onMoleculeLoad).toHaveBeenCalledWith({
      type: 'coordinates',
      coordinates: [
        ['H', 0.0, 0.0, 0.0],
        ['H', 0.74, 0.0, 0.0]
      ],
      charge: 0,
      spin: 0
    })
  })

  it('renders molecule atoms emoji', () => {
    render(<ProjectPanel {...mockProps} />)
    
    // Each molecule should have the atom emoji
    const atomEmojis = screen.getAllByText('⚛️')
    expect(atomEmojis).toHaveLength(3) // One for each sample molecule
  })
})