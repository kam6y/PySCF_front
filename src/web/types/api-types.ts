/**
 * Convenience exports from generated OpenAPI types
 * This file provides easier access to the generated types for use throughout the application
 */

import { components } from './generated-api';
import { ApiError } from '../api/core';

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
export type PauseCalculationResponse =
  components['schemas']['PauseCalculationResponse'];
export type ResumeCalculationResponse =
  components['schemas']['ResumeCalculationResponse'];
export type OrbitalsResponse = components['schemas']['OrbitalsResponse'];
export type OrbitalCubeResponse = components['schemas']['OrbitalCubeResponse'];
export type SupportedParametersResponse =
  components['schemas']['SupportedParametersResponse'];
export type IRSpectrumResponse = components['schemas']['IRSpectrumResponse'];
export type ErrorResponse = components['schemas']['ErrorResponse'];
export type Gpu4PyscfStatus = components['schemas']['Gpu4PyscfStatus'];
export type Gpu4PyscfStatusResponse =
  components['schemas']['Gpu4PyscfStatusResponse'];
export type Gpu4PyscfInstallRequest =
  components['schemas']['Gpu4PyscfInstallRequest'];
export type Gpu4PyscfInstallResult =
  components['schemas']['Gpu4PyscfInstallResult'];
export type Gpu4PyscfInstallResponse =
  components['schemas']['Gpu4PyscfInstallResponse'];

// Data models
export type CalculationParameters =
  components['schemas']['CalculationParameters'];
export type CalculationResults = components['schemas']['CalculationResults'];
export type CalculationInstance = components['schemas']['CalculationInstance'];
export type CalculationSummary = components['schemas']['CalculationSummary'];
export type OrbitalInfo = components['schemas']['OrbitalInfo'];
export type IRSpectrumData = components['schemas']['IRSpectrumData'];
export type IRSpectrumDetails = components['schemas']['IRSpectrumDetails'];
export type IRSpectrumMetadata = components['schemas']['IRSpectrumMetadata'];

// Chat history data models
export type ChatSession = components['schemas']['ChatSession'];
export type ChatSessionSummary = components['schemas']['ChatSessionSummary'];
export type ChatSessionDetail = components['schemas']['ChatSessionDetail'];

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
export type PauseCalculationResponseData = PauseCalculationResponse['data'];
export type ResumeCalculationResponseData = ResumeCalculationResponse['data'];
export type OrbitalsResponseData = OrbitalsResponse['data'];
export type OrbitalCubeResponseData = OrbitalCubeResponse['data'];
export type SupportedParametersResponseData =
  SupportedParametersResponse['data'];
export type IRSpectrumResponseData = IRSpectrumResponse['data'];
export type Gpu4PyscfStatusResponseData = Gpu4PyscfStatusResponse['data'];
export type Gpu4PyscfInstallResponseData = Gpu4PyscfInstallResponse['data'];

// Error handling
export { ApiError };
