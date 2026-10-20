/**
 * Convenience exports from generated OpenAPI types
 * This file provides easier access to the generated types for use throughout the application
 */

import type { components } from './generated-api';

// Export enum types as union types for convenience
export type CalculationStatus = components['schemas']['CalculationStatus'];

// Request types
export type QuantumCalculationRequest =
  components['schemas']['QuantumCalculationRequest'];

// Data models
export type CalculationParameters =
  components['schemas']['CalculationParameters'];
export type CalculationResults = components['schemas']['CalculationResults'];
export type CalculationInstance = components['schemas']['CalculationInstance'];
export type CalculationSummary = components['schemas']['CalculationSummary'];
export type OrbitalInfo = components['schemas']['OrbitalInfo'];

// Chat history data models
export type ChatSessionSummary = components['schemas']['ChatSessionSummary'];

// Commonly used response data types
export type PubChemSearchResponseData =
  components['schemas']['PubChemSearchResponse']['data'];
export type SMILESConvertResponseData =
  components['schemas']['SMILESConvertResponse']['data'];
export type StartCalculationResponseData =
  components['schemas']['StartCalculationResponse']['data'];
export type CalculationListResponseData =
  components['schemas']['CalculationListResponse']['data'];
export type CalculationDetailsResponseData =
  components['schemas']['CalculationDetailsResponse']['data'];
export type CalculationUpdateResponseData =
  components['schemas']['CalculationUpdateResponse']['data'];
export type CalculationDeleteResponseData =
  components['schemas']['CalculationDeleteResponse']['data'];
export type PauseCalculationResponseData =
  components['schemas']['PauseCalculationResponse']['data'];
export type ResumeCalculationResponseData =
  components['schemas']['ResumeCalculationResponse']['data'];
export type OrbitalsResponseData =
  components['schemas']['OrbitalsResponse']['data'];
export type OrbitalCubeResponseData =
  components['schemas']['OrbitalCubeResponse']['data'];
export type SupportedParametersResponseData =
  components['schemas']['SupportedParametersResponse']['data'];
export type IRSpectrumResponseData =
  components['schemas']['IRSpectrumResponse']['data'];
