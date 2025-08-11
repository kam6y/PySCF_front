/**
 * Legacy type definitions for quantum chemistry calculations
 * 
 * This file now re-exports from the OpenAPI-generated types to maintain compatibility
 * while ensuring all types come from a single source of truth.
 * 
 * @deprecated Use types from './api-types' directly instead
 */

// Re-export the main types from the new OpenAPI-generated types
export {
  CalculationInstance,
  CalculationStatus,
  CalculationMethod,
  SolventMethod,
  CalculationResults,
  CalculationSummary,
  CalculationListResponseData as CalculationListResponse,
  CalculationDetailsResponseData as CalculationDetailsResponse,
  CalculationUpdateResponseData as RenameResponse,
  CalculationDeleteResponseData as DeleteResponse,
  // Legacy compatibility: map old CalculationParameters to new request type
  QuantumCalculationRequest as CalculationParameters,
} from './api-types';

// Legacy compatibility aliases
export type BasisFunction = string;
export type ExchangeCorrelation = string;