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
   *
   * This function implements the Single Source of Truth principle by completely
   * trusting the backend's method_defaults. The backend already provides only
   * the applicable parameters for each method, so no frontend filtering is needed.
   *
   * @param currentParams - Current calculation parameters
   * @param newMethod - New calculation method to switch to
   * @returns Updated parameters with method defaults applied
   *
   * @example
   * const updated = applyMethodDefaults(
   *   { calculation_method: 'DFT', basis_function: '6-31G(d)', exchange_correlation: 'B3LYP', ... },
   *   'HF'
   * );
   * // Returns: { calculation_method: 'HF', basis_function: '6-31G(d)', ... }
   * // Note: exchange_correlation is automatically removed because it's not in HF's defaults
   *
   * @remarks
   * Design principles:
   * - Backend method_defaults is the single source of truth
   * - No hardcoded parameter lists in frontend
   * - Molecular and system settings are preserved across method changes
   */
  const applyMethodDefaults = (
    currentParams: Partial<QuantumCalculationRequest> | CalculationParameters,
    newMethod: string
  ): Partial<QuantumCalculationRequest> => {
    const defaults = getDefaultsForMethod(newMethod);
    const params = currentParams as any;

    // Preserve molecular structure and system settings that transcend calculation methods
    const preservedParams: Record<string, any> = {
      xyz: params.xyz || '',
      name: params.name || params.molecule_name || '',
      charges: params.charges ?? 0,
      spin: params.spin ?? 0,
      solvent_method: params.solvent_method || 'none',
      solvent: params.solvent || '-',
      ketcher_data: params.ketcher_data,
    };

    // Preserve cpu_cores if explicitly set by user
    if (params.cpu_cores !== undefined) {
      preservedParams.cpu_cores = params.cpu_cores;
    }

    // Trust backend defaults completely - they already contain only applicable parameters
    // Preserved params override defaults to maintain user's molecular structure and settings
    return {
      ...defaults,
      ...preservedParams,
      calculation_method: newMethod as QuantumCalculationRequest['calculation_method'],
    } as Partial<QuantumCalculationRequest>;
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
