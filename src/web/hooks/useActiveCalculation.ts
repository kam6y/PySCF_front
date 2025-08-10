import { useState, useEffect, useCallback } from 'react';
import { CalculationInstance, CalculationDetailsResponse } from '../types/calculation';

const ACTIVE_CALCULATION_KEY = 'pyscf-active-calculation-id';
const API_BASE_URL = 'http://127.0.0.1:5000';

export interface UseActiveCalculationReturn {
  activeCalculationId: string | null;
  activeCalculation: CalculationInstance | undefined;
  setActiveCalculation: (calculation: CalculationInstance | null) => void;
  setActiveCalculationById: (id: string | null) => void;
  loadCalculationDetails: (id: string) => Promise<CalculationInstance | null>;
  clearActiveCalculation: () => void;
}

export const useActiveCalculation = (
  calculations: CalculationInstance[]
): UseActiveCalculationReturn => {
  const [activeCalculationId, setActiveCalculationId] = useState<string | null>(() => {
    // Restore active calculation ID from localStorage
    try {
      return localStorage.getItem(ACTIVE_CALCULATION_KEY);
    } catch {
      return null;
    }
  });

  // Get active calculation from the calculations array
  const activeCalculation = activeCalculationId 
    ? calculations.find(calc => calc.id === activeCalculationId)
    : undefined;

  // Persist active calculation ID to localStorage
  useEffect(() => {
    try {
      if (activeCalculationId) {
        localStorage.setItem(ACTIVE_CALCULATION_KEY, activeCalculationId);
      } else {
        localStorage.removeItem(ACTIVE_CALCULATION_KEY);
      }
    } catch (err) {
      console.error('Failed to persist active calculation ID:', err);
    }
  }, [activeCalculationId]);

  // Load detailed calculation data from backend
  const loadCalculationDetails = useCallback(async (id: string): Promise<CalculationInstance | null> => {
    try {
      // For now, we'll use the existing calculation info and try to parse stored results
      const calculation = calculations.find(calc => calc.id === id);
      if (!calculation) {
        return null;
      }

      // Try to load detailed results if calculation is completed
      if (calculation.status === 'completed' && calculation.workingDirectory) {
        // TODO: Implement detailed calculation loading from backend when API is available
        // const response = await fetch(`${API_BASE_URL}/api/quantum/calculations/${id}`);
        // const data: CalculationDetailsResponse = await response.json();
        
        // For now, return the calculation as-is
        return calculation;
      }

      return calculation;
    } catch (err) {
      console.error('Failed to load calculation details:', err);
      return null;
    }
  }, [calculations]);

  // Set active calculation directly
  const setActiveCalculation = useCallback((calculation: CalculationInstance | null) => {
    setActiveCalculationId(calculation?.id || null);
  }, []);

  // Set active calculation by ID
  const setActiveCalculationById = useCallback((id: string | null) => {
    setActiveCalculationId(id);
  }, []);

  // Clear active calculation
  const clearActiveCalculation = useCallback(() => {
    setActiveCalculationId(null);
  }, []);

  // Auto-clear active calculation if it no longer exists in calculations list
  useEffect(() => {
    if (activeCalculationId && !calculations.find(calc => calc.id === activeCalculationId)) {
      setActiveCalculationId(null);
    }
  }, [activeCalculationId, calculations]);

  return {
    activeCalculationId,
    activeCalculation,
    setActiveCalculation,
    setActiveCalculationById,
    loadCalculationDetails,
    clearActiveCalculation
  };
};