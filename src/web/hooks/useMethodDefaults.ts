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
import {
  QuantumCalculationRequest,
  CalculationParameters,
} from '../types/api-types';

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
  const isParameterApplicable = (
    paramName: string,
    method: string
  ): boolean => {
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
   * This function uses `applicable_methods` constraint to determine if a parameter
   * should be disabled. If a method is NOT in the `applicable_methods` list,
   * the parameter is disabled for that method.
   *
   * @param paramName - Name of the parameter (e.g., 'optimize_geometry')
   * @param method - Calculation method name
   * @returns True if the parameter should be disabled for the method
   *
   * @example
   * isParameterDisabled('optimize_geometry', 'TDDFT');  // Returns: true (TDDFT not in ['DFT', 'HF', 'MP2'])
   * isParameterDisabled('optimize_geometry', 'DFT');    // Returns: false (DFT is in ['DFT', 'HF', 'MP2'])
   * isParameterDisabled('cpu_cores', 'DFT');            // Returns: false (no constraint, universal parameter)
   */
  const isParameterDisabled = (paramName: string, method: string): boolean => {
    const constraint = supportedParams?.parameter_constraints?.[paramName];
    if (!constraint) return false;

    // Check applicable_methods constraint (inverted logic)
    // If method is NOT in the applicable_methods list, the parameter is disabled
    if (constraint.applicable_methods) {
      return !constraint.applicable_methods.includes(method);
    }

    // No constraint means parameter is universally applicable
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
   * trusting the backend's method_defaults. It also removes any parameters that
   * are not applicable to the new method to prevent backend validation errors.
   *
   * @param currentParams - Current calculation parameters
   * @param newMethod - New calculation method to switch to
   * @returns Updated parameters with method defaults applied and inapplicable parameters removed
   *
   * @example
   * const updated = applyMethodDefaults(
   *   { calculation_method: 'DFT', basis_function: '6-31G(d)', optimize_geometry: true, ... },
   *   'TDDFT'
   * );
   * // Returns: { calculation_method: 'TDDFT', basis_function: '6-31G(d)', ... }
   * // Note: optimize_geometry is automatically removed because it's not applicable to TDDFT
   *
   * @remarks
   * Design principles:
   * - Backend method_defaults is the single source of truth
   * - Explicitly remove inapplicable parameters to prevent validation errors
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

    // Build final parameters from defaults and preserved params
    const finalParams = {
      ...defaults,
      ...preservedParams,
      calculation_method:
        newMethod as QuantumCalculationRequest['calculation_method'],
    };

    // IMPORTANT: Remove any parameters that are not applicable to the new method
    // This prevents backend validation errors when switching methods
    const cleanedParams: Record<string, any> = {};
    for (const [key, value] of Object.entries(finalParams)) {
      // Keep the parameter if:
      // 1. It has no constraint (universal parameter)
      // 2. It is applicable to the new method
      if (isParameterApplicable(key, newMethod)) {
        cleanedParams[key] = value;
      }
    }

    return cleanedParams as Partial<QuantumCalculationRequest>;
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
