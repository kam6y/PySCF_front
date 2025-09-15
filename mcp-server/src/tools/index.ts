// Export all tools and their handlers
export * from './pubchem.js';
export * from './quantum.js';
export * from './system.js';
export * from './inspector.js';

import {
  searchPubChemTool,
  convertSmilesTool,
  validateXyzTool,
  handleSearchPubChem,
  handleConvertSmiles,
  handleValidateXYZ,
} from './pubchem.js';

import {
  startCalculationTool,
  listCalculationsTool,
  getCalculationDetailsTool,
  getOrbitalsTool,
  getOrbitalCubeTool,
  getIRSpectrumTool,
  getOptimizedGeometryTool,
  startStepwiseCalculationTool,
  handleStartCalculation,
  handleListCalculations,
  handleGetCalculationDetails,
  handleGetOrbitals,
  handleGetOrbitalCube,
  handleGetIRSpectrum,
  handleGetOptimizedGeometry,
  handleStartStepwiseCalculation,
} from './quantum.js';

import {
  getSupportedParametersTool,
  getSettingsTool,
  updateSettingsTool,
  getResourceStatusTool,
  testConnectionTool,
  diagnosticsServerTool,
  handleGetSupportedParameters,
  handleGetSettings,
  handleUpdateSettings,
  handleGetResourceStatus,
  handleTestConnection,
  handleDiagnosticsServer,
} from './system.js';

import {
  diagnosticsTool,
  testApiTool,
  getDebugLogsTool,
  validateConfigTool,
  handleDiagnostics,
  handleTestApi,
  handleGetDebugLogs,
  handleValidateConfig,
} from './inspector.js';

import { Tool } from '@modelcontextprotocol/sdk/types.js';
import { PySCFApiClient } from '../client.js';

// All available tools
export const tools: Tool[] = [
  // PubChem and SMILES tools
  searchPubChemTool,
  convertSmilesTool,
  validateXyzTool,
  
  // Quantum calculation tools
  startCalculationTool,
  listCalculationsTool,
  getCalculationDetailsTool,
  getOrbitalsTool,
  getOrbitalCubeTool,
  getIRSpectrumTool,
  getOptimizedGeometryTool,
  startStepwiseCalculationTool,
  
  // System management tools
  getSupportedParametersTool,
  getSettingsTool,
  updateSettingsTool,
  getResourceStatusTool,
  testConnectionTool,
  diagnosticsServerTool,
  
  // Inspector tools
  diagnosticsTool,
  testApiTool,
  getDebugLogsTool,
  validateConfigTool,
];

// Tool handlers map
export const toolHandlers: Record<string, (args: any, client: PySCFApiClient) => Promise<any>> = {
  // PubChem and SMILES handlers
  searchPubChem: handleSearchPubChem,
  convertSmiles: handleConvertSmiles,
  validateXYZ: handleValidateXYZ,
  
  // Quantum calculation handlers
  startCalculation: handleStartCalculation,
  listCalculations: handleListCalculations,
  getCalculationDetails: handleGetCalculationDetails,
  getOrbitals: handleGetOrbitals,
  getOrbitalCube: handleGetOrbitalCube,
  getIRSpectrum: handleGetIRSpectrum,
  getOptimizedGeometry: handleGetOptimizedGeometry,
  startStepwiseCalculation: handleStartStepwiseCalculation,
  
  // System management handlers
  getSupportedParameters: handleGetSupportedParameters,
  getSettings: handleGetSettings,
  updateSettings: handleUpdateSettings,
  getResourceStatus: handleGetResourceStatus,
  testConnection: handleTestConnection,
  diagnosticsServer: handleDiagnosticsServer,
  
  // Inspector handlers
  diagnostics: handleDiagnostics,
  testApi: handleTestApi,
  getDebugLogs: handleGetDebugLogs,
  validateConfig: handleValidateConfig,
};