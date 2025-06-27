import { contextBridge, ipcRenderer } from 'electron'

// Define the API that will be exposed to the renderer process
export interface ElectronAPI {
  // Calculation methods
  calculate: (calculationData: any) => Promise<any>
  buildMolecule: (moleculeData: any) => Promise<any>
  getCalculationStatus: () => Promise<any>
  stopCalculation: () => Promise<any>

  // Event listeners
  onPythonLog: (callback: (data: { type: string; data: string }) => void) => void
  onPythonError: (callback: (error: string) => void) => void
  onCalculationProgress: (callback: (progress: any) => void) => void

  // Utility methods
  removeAllListeners: (channel: string) => void
}

// Expose protected methods that allow the renderer process to use
// the ipcRenderer without exposing the entire object
const electronAPI: ElectronAPI = {
  // Calculation methods
  calculate: (calculationData: any) => ipcRenderer.invoke('calculate', calculationData),
  buildMolecule: (moleculeData: any) => ipcRenderer.invoke('build-molecule', moleculeData),
  getCalculationStatus: () => ipcRenderer.invoke('get-calculation-status'),
  stopCalculation: () => ipcRenderer.invoke('stop-calculation'),

  // Event listeners
  onPythonLog: (callback) => {
    ipcRenderer.on('python-log', (_event, data) => callback(data))
  },
  onPythonError: (callback) => {
    ipcRenderer.on('python-error', (_event, error) => callback(error))
  },
  onCalculationProgress: (callback) => {
    ipcRenderer.on('calculation-progress', (_event, progress) => callback(progress))
  },

  // Utility methods
  removeAllListeners: (channel) => {
    ipcRenderer.removeAllListeners(channel)
  },
}

contextBridge.exposeInMainWorld('electronAPI', electronAPI)

// TypeScript declaration for global window object
declare global {
  interface Window {
    electronAPI: ElectronAPI
  }
}