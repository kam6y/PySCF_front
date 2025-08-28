/**
 * Convenience exports from generated OpenAPI types
 * This file provides easier access to the generated types for use throughout the application
 */

import { components, operations } from './generated-api';
import { ApiError } from '../apiClient';

// Export enum types as union types for convenience
export type SearchType = components['schemas']['SearchType'];
export type SolventMethod = components['schemas']['SolventMethod'];
export type CalculationMethod = components['schemas']['CalculationMethod'];
export type CalculationStatus = components['schemas']['CalculationStatus'];

// Request types
export type PubChemSearchRequest =
  components['schemas']['PubChemSearchRequest'];
export type SMILESConvertRequest =
  components['schemas']['SMILESConvertRequest'];
export type XYZValidateRequest = components['schemas']['XYZValidateRequest'];
export type QuantumCalculationRequest =
  components['schemas']['QuantumCalculationRequest'];
export type CalculationUpdateRequest =
  components['schemas']['CalculationUpdateRequest'];

// Response types
export type HealthResponse = components['schemas']['HealthResponse'];
export type PubChemCompoundInfo = components['schemas']['PubChemCompoundInfo'];
export type PubChemSearchResponse =
  components['schemas']['PubChemSearchResponse'];
export type SMILESConvertResponse =
  components['schemas']['SMILESConvertResponse'];
export type XYZValidateResponse = components['schemas']['XYZValidateResponse'];
export type StartCalculationResponse =
  components['schemas']['StartCalculationResponse'];
export type CalculationListResponse =
  components['schemas']['CalculationListResponse'];
export type CalculationDetailsResponse =
  components['schemas']['CalculationDetailsResponse'];
export type CalculationUpdateResponse =
  components['schemas']['CalculationUpdateResponse'];
export type CalculationDeleteResponse =
  components['schemas']['CalculationDeleteResponse'];
export type OrbitalsResponse = components['schemas']['OrbitalsResponse'];
export type OrbitalCubeResponse = components['schemas']['OrbitalCubeResponse'];
export type CubeFilesListResponse =
  components['schemas']['CubeFilesListResponse'];
export type CubeFilesDeleteResponse =
  components['schemas']['CubeFilesDeleteResponse'];
export type SupportedParametersResponse =
  components['schemas']['SupportedParametersResponse'];
export type ErrorResponse = components['schemas']['ErrorResponse'];

// Data models
export type CalculationParameters =
  components['schemas']['CalculationParameters'];
export type CalculationResults = components['schemas']['CalculationResults'];
export type CalculationInstance = components['schemas']['CalculationInstance'];
export type CalculationSummary = components['schemas']['CalculationSummary'];
export type OrbitalInfo = components['schemas']['OrbitalInfo'];

// Operation types for reference
export type Operations = {
  healthCheck: operations['healthCheck'];
  searchPubChem: operations['searchPubChem'];
  convertSmiles: operations['convertSmiles'];
  validateXYZ: operations['validateXYZ'];
  startCalculation: operations['startCalculation'];
  listCalculations: operations['listCalculations'];
  getCalculationDetails: operations['getCalculationDetails'];
  updateCalculation: operations['updateCalculation'];
  deleteCalculation: operations['deleteCalculation'];
  getOrbitals: operations['getOrbitals'];
  getOrbitalCube: operations['getOrbitalCube'];
  listCubeFiles: operations['listCubeFiles'];
  deleteCubeFiles: operations['deleteCubeFiles'];
  getSupportedParameters: operations['getSupportedParameters'];
};

// Response wrapper for API responses
export type ApiResponse<T> = {
  success: boolean;
  data: T;
  error?: string;
};

// Commonly used response data types
export type PubChemSearchResponseData = PubChemSearchResponse['data'];
export type SMILESConvertResponseData = SMILESConvertResponse['data'];
export type StartCalculationResponseData = StartCalculationResponse['data'];
export type CalculationListResponseData = CalculationListResponse['data'];
export type CalculationDetailsResponseData = CalculationDetailsResponse['data'];
export type CalculationUpdateResponseData = CalculationUpdateResponse['data'];
export type CalculationDeleteResponseData = CalculationDeleteResponse['data'];
export type OrbitalsResponseData = OrbitalsResponse['data'];
export type OrbitalCubeResponseData = OrbitalCubeResponse['data'];
export type CubeFilesListResponseData = CubeFilesListResponse['data'];
export type CubeFilesDeleteResponseData = CubeFilesDeleteResponse['data'];
export type SupportedParametersResponseData =
  SupportedParametersResponse['data'];

// Error handling
export { ApiError };
