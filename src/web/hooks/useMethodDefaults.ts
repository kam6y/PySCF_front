/**
 * Custom hook for accessing calculation method defaults and parameter constraints.
 *
 * This hook provides a centralized interface for:
 * - Retrieving default parameter values for each calculation method
 * - Checking parameter applicability and constraints
 * - Applying method-specific defaults when switching calculation methods
 *
 * All business logic related to defaults and constraints is sourced from the backend
 * via the `/api/quantum/supported-parameters` endpoint, ensuring a single source of truth.
 */

import { useSupportedParameters } from './useCalculationQueries';
import { QuantumCalculationRequest, CalculationParameters } from '../types/api-types';

export const useMethodDefaults = () => {
  const { data: supportedParams } = useSupportedParameters();

  /**
   * Get default parameter values for a specific calculation method.
   *
   * @param method - Calculation method name (e.g., 'DFT', 'CCSD', 'TDDFT')
   * @returns Object containing default parameter values for the method
   *
   * @example
   * const defaults = getDefaultsForMethod('CCSD');
   * // Returns: { basis_function: 'cc-pVDZ', memory_mb: 4000, frozen_core: true, ... }
   */
  const getDefaultsForMethod = (
    method: string
  ): Partial<QuantumCalculationRequest> => {
    if (!supportedParams?.method_defaults) return {};
    return (supportedParams.method_defaults[method] ||
      {}) as Partial<QuantumCalculationRequest>;
  };

  /**
   * Check if a parameter is applicable for a specific calculation method.
   *
   * @param paramName - Name of the parameter (e.g., 'ncas', 'tddft_nstates')
   * @param method - Calculation method name
   * @returns True if the parameter is applicable for the method
   *
   * @example
   * isParameterApplicable('ncas', 'CASCI');  // Returns: true
   * isParameterApplicable('ncas', 'DFT');    // Returns: false
   */
  const isParameterApplicable = (paramName: string, method: string): boolean => {
    const constraint = supportedParams?.parameter_constraints?.[paramName];
    if (!constraint) return true;

    if (constraint.applicable_methods) {
      return constraint.applicable_methods.includes(method);
    }
    return true;
  };

  /**
   * Check if a parameter should be disabled for a specific calculation method.
   *
   * @param paramName - Name of the parameter (e.g., 'optimize_geometry')
   * @param method - Calculation method name
   * @returns True if the parameter should be disabled for the method
   *
   * @example
   * isParameterDisabled('optimize_geometry', 'TDDFT');  // Returns: true
   * isParameterDisabled('optimize_geometry', 'DFT');    // Returns: false
   */
  const isParameterDisabled = (paramName: string, method: string): boolean => {
    const constraint = supportedParams?.parameter_constraints?.[paramName];
    if (!constraint) return false;

    if (constraint.disabled_for) {
      return constraint.disabled_for.includes(method);
    }
    return false;
  };

  /**
   * Get constraint information for a parameter.
   *
   * @param paramName - Name of the parameter
   * @returns Constraint object containing min, max, applicable_methods, etc., or undefined
   *
   * @example
   * const constraint = getParameterConstraint('ncas');
   * // Returns: { min: 1, max: 20, applicable_methods: ['CASCI', 'CASSCF'], ... }
   */
  const getParameterConstraint = (paramName: string) => {
    return supportedParams?.parameter_constraints?.[paramName];
  };

  /**
   * Apply method-specific default values when switching calculation methods.
   * This completely replaces parameter values with the method's defaults.
   *
   * @param currentParams - Current calculation parameters
   * @param newMethod - New calculation method to switch to
   * @returns Updated parameters with method defaults applied
   *
   * @example
   * const updated = applyMethodDefaults(
   *   { calculation_method: 'DFT', basis_function: '6-31G(d)', ... },
   *   'CCSD'
   * );
   * // Returns: { calculation_method: 'CCSD', basis_function: 'cc-pVDZ', memory_mb: 4000, ... }
   *
   * @remarks
   * This function follows the user requirement to ALWAYS overwrite with defaults when
   * the calculation method changes, providing a consistent starting point for each method.
   */
  const applyMethodDefaults = (
    currentParams: Partial<QuantumCalculationRequest> | CalculationParameters,
    newMethod: string
  ): QuantumCalculationRequest => {
    const defaults = getDefaultsForMethod(newMethod);
    const params = currentParams as any;

    // Preserve core parameters that should not be overwritten
    const preservedParams = {
      xyz: params.xyz || '',
      name: params.name || params.molecule_name || '',
      charges: params.charges ?? 0,
      spin: params.spin ?? 0,
      solvent_method: params.solvent_method || 'none',
      solvent: params.solvent || '-',
      cpu_cores: params.cpu_cores,
      ketcher_data: params.ketcher_data,
    };

    // Apply defaults for the new method, overwriting everything except preserved params
    // Ensure all required fields are present
    const result: QuantumCalculationRequest = {
      ...params,
      ...defaults,
      ...preservedParams,
      calculation_method: newMethod as any,
      basis_function: (defaults.basis_function as string) || params.basis_function || '6-31G(d)',
      exchange_correlation: (defaults.exchange_correlation as string | null) ?? params.exchange_correlation ?? null,
      optimize_geometry: defaults.optimize_geometry ?? params.optimize_geometry ?? true,

      // Ensure required fields have defaults (handling null/undefined from CalculationParameters)
      tddft_nstates: (defaults.tddft_nstates ?? params.tddft_nstates ?? 10) as number,
      tddft_method: (defaults.tddft_method ?? params.tddft_method ?? 'TDDFT') as "TDDFT" | "TDA",
      tddft_analyze_nto: defaults.tddft_analyze_nto ?? params.tddft_analyze_nto ?? false,
      ncas: (defaults.ncas ?? params.ncas ?? 4) as number,
      nelecas: (defaults.nelecas ?? params.nelecas ?? 4) as number,
      max_cycle_macro: (defaults.max_cycle_macro ?? params.max_cycle_macro ?? 50) as number,
      max_cycle_micro: (defaults.max_cycle_micro ?? params.max_cycle_micro ?? 3) as number,
      natorb: defaults.natorb ?? params.natorb ?? true,
      conv_tol: (defaults.conv_tol ?? params.conv_tol ?? 0.000001) as number,
      conv_tol_grad: (defaults.conv_tol_grad ?? params.conv_tol_grad ?? 0.0001) as number,
    };

    return result;
  };

  /**
   * Check if supported parameters data is currently loading.
   */
  const isLoading = !supportedParams;

  return {
    getDefaultsForMethod,
    isParameterApplicable,
    isParameterDisabled,
    getParameterConstraint,
    applyMethodDefaults,
    isLoading,
  };
};
