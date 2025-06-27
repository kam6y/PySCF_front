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

declare global {
  interface Window {
    electronAPI: ElectronAPI
  }
}