/**
 * Unit tests for StatusBar component
 */

import React from 'react'
import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import '@testing-library/jest-dom'
import StatusBar from '../../src/renderer/components/StatusBar'

// Mock props
const mockProps = {
  isCalculating: false,
  logs: [],
  onClearLogs: jest.fn(),
}

describe('StatusBar', () => {
  beforeEach(() => {
    jest.clearAllMocks()
  })

  it('renders status bar with basic elements', () => {
    render(<StatusBar {...mockProps} />)
    
    expect(screen.getByText('Ready')).toBeInTheDocument()
    expect(screen.getByText('PySCF Backend Connected')).toBeInTheDocument()
    expect(screen.getByText(/logs/i)).toBeInTheDocument()
  })

  it('shows calculating state', () => {
    const calculatingProps = { ...mockProps, isCalculating: true }
    render(<StatusBar {...calculatingProps} />)
    
    expect(screen.getByText('Calculating...')).toBeInTheDocument()
  })

  it('shows error state when error logs present', () => {
    const errorProps = {
      ...mockProps,
      logs: [
        { type: 'error', data: 'Some error occurred' }
      ]
    }
    render(<StatusBar {...errorProps} />)
    
    expect(screen.getByText('Error')).toBeInTheDocument()
  })

  it('shows error state when stderr logs present', () => {
    const stderrProps = {
      ...mockProps,
      logs: [
        { type: 'stderr', data: 'Some stderr message' }
      ]
    }
    render(<StatusBar {...stderrProps} />)
    
    expect(screen.getByText('Error')).toBeInTheDocument()
  })

  it('shows log count when logs are present', () => {
    const logsProps = {
      ...mockProps,
      logs: [
        { type: 'stdout', data: 'Log message 1' },
        { type: 'stdout', data: 'Log message 2' },
        { type: 'stdout', data: 'Log message 3' }
      ]
    }
    render(<StatusBar {...logsProps} />)
    
    expect(screen.getByText('(3)')).toBeInTheDocument()
  })

  it('toggles log panel when logs button clicked', async () => {
    const user = userEvent.setup()
    const logsProps = {
      ...mockProps,
      logs: [{ type: 'stdout', data: 'Test log' }]
    }
    render(<StatusBar {...logsProps} />)
    
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    expect(screen.getByText('System Logs')).toBeInTheDocument()
    expect(screen.getByText('Test log')).toBeInTheDocument()
  })

  it('shows log panel with correct log content', async () => {
    const user = userEvent.setup()
    const logsProps = {
      ...mockProps,
      logs: [
        { type: 'stdout', data: 'Standard output message' },
        { type: 'stderr', data: 'Error message' },
        { type: 'error', data: 'Exception occurred' }
      ]
    }
    render(<StatusBar {...logsProps} />)
    
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    expect(screen.getByText('Standard output message')).toBeInTheDocument()
    expect(screen.getByText('Error message')).toBeInTheDocument()
    expect(screen.getByText('Exception occurred')).toBeInTheDocument()
  })

  it('calls onClearLogs when clear button clicked', async () => {
    const user = userEvent.setup()
    const logsProps = {
      ...mockProps,
      logs: [{ type: 'stdout', data: 'Test log' }]
    }
    render(<StatusBar {...logsProps} />)
    
    // Open log panel
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    // Click clear button
    const clearButton = screen.getByText('Clear')
    await user.click(clearButton)
    
    expect(mockProps.onClearLogs).toHaveBeenCalled()
  })

  it('closes log panel when close button clicked', async () => {
    const user = userEvent.setup()
    const logsProps = {
      ...mockProps,
      logs: [{ type: 'stdout', data: 'Test log' }]
    }
    render(<StatusBar {...logsProps} />)
    
    // Open log panel
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    expect(screen.getByText('System Logs')).toBeInTheDocument()
    
    // Close log panel
    const closeButton = screen.getByText('✕')
    await user.click(closeButton)
    
    expect(screen.queryByText('System Logs')).not.toBeInTheDocument()
  })

  it('shows empty logs message when no logs', async () => {
    const user = userEvent.setup()
    render(<StatusBar {...mockProps} />)
    
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    expect(screen.getByText('No logs yet...')).toBeInTheDocument()
  })

  it('applies correct text colors for different log types', async () => {
    const user = userEvent.setup()
    const logsProps = {
      ...mockProps,
      logs: [
        { type: 'stdout', data: 'Info message' },
        { type: 'stderr', data: 'Error message' },
        { type: 'error', data: 'Exception' }
      ]
    }
    render(<StatusBar {...logsProps} />)
    
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    const infoLog = screen.getByText('Info message')
    const errorLog = screen.getByText('Error message')
    const exceptionLog = screen.getByText('Exception')
    
    expect(infoLog).toHaveClass('text-info')
    expect(errorLog).toHaveClass('text-error')
    expect(exceptionLog).toHaveClass('text-error')
  })

  it('shows calculating indicator with pulsing animation', () => {
    const calculatingProps = { ...mockProps, isCalculating: true }
    render(<StatusBar {...calculatingProps} />)
    
    const indicator = screen.getByText('Calculating...').previousSibling
    expect(indicator).toHaveClass('bg-warning', 'animate-pulse')
  })

  it('shows success indicator when ready', () => {
    render(<StatusBar {...mockProps} />)
    
    const indicator = screen.getByText('Ready').previousSibling
    expect(indicator).toHaveClass('bg-success')
  })

  it('shows error indicator when errors present', () => {
    const errorProps = {
      ...mockProps,
      logs: [{ type: 'error', data: 'Error' }]
    }
    render(<StatusBar {...errorProps} />)
    
    const indicator = screen.getByText('Error').previousSibling
    expect(indicator).toHaveClass('bg-error')
  })

  it('updates logs button arrow direction', async () => {
    const user = userEvent.setup()
    render(<StatusBar {...mockProps} />)
    
    // Initially should show up arrow
    expect(screen.getByText(/logs ▲/i)).toBeInTheDocument()
    
    // Click to open
    const logsButton = screen.getByText(/logs/i)
    await user.click(logsButton)
    
    // Should now show down arrow
    expect(screen.getByText(/logs ▼/i)).toBeInTheDocument()
  })
})