// src/web/hooks/useActiveCalculation.ts

import { useState, useEffect, useCallback } from 'react';
import { CalculationInstance } from '../types/calculation';
import { getCalculationDetails } from '../apiClient'; // apiClientからインポート

const ACTIVE_CALCULATION_KEY = 'pyscf-active-calculation-id';

export interface UseActiveCalculationReturn {
  activeCalculationId: string | null;
  activeCalculation: CalculationInstance | undefined;
  isLoadingDetails: boolean;
  detailsError: string | null;
  setActiveCalculation: (calculation: CalculationInstance | null) => void;
  setActiveCalculationById: (id: string | null) => void;
  loadCalculationDetails: (id: string, forceRefresh?: boolean) => Promise<CalculationInstance | null>;
  clearActiveCalculation: () => void;
}

export const useActiveCalculation = (
  calculations: CalculationInstance[]
): UseActiveCalculationReturn => {
  const [activeCalculationId, setActiveCalculationId] = useState<string | null>(() => {
    try {
      const storedId = localStorage.getItem(ACTIVE_CALCULATION_KEY);
      // Guard against invalid strings from localStorage
      if (storedId && storedId !== "undefined" && storedId !== "null") {
        return storedId;
      }
      return null;
    } catch {
      return null;
    }
  });

  const [detailedCalculations, setDetailedCalculations] = useState<Map<string, CalculationInstance>>(new Map());
  const [isLoadingDetails, setIsLoadingDetails] = useState(false);
  const [detailsError, setDetailsError] = useState<string | null>(null);

  const activeCalculation = activeCalculationId
    ? detailedCalculations.get(activeCalculationId) || calculations.find(calc => calc.id === activeCalculationId)
    : undefined;

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

  const loadCalculationDetails = useCallback(async (id: string, forceRefresh = false): Promise<CalculationInstance | null> => {
    if (!forceRefresh && detailedCalculations.has(id)) {
      return detailedCalculations.get(id) || null;
    }

    setIsLoadingDetails(true);
    setDetailsError(null);

    try {
      const data = await getCalculationDetails(id);
      const detailedCalculation = data.data.calculation;
      
      setDetailedCalculations(prev => {
        const newMap = new Map(prev);
        newMap.set(id, detailedCalculation);
        return newMap;
      });
      
      return detailedCalculation;
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'An unknown error occurred.';
      setDetailsError(errorMessage);
      console.error('Failed to load calculation details:', err);
      return null;
    } finally {
      setIsLoadingDetails(false);
    }
  }, [detailedCalculations]);

  const setActiveCalculation = useCallback((calculation: CalculationInstance | null) => {
    setActiveCalculationId(calculation?.id || null);
  }, []);

  const setActiveCalculationById = useCallback((id: string | null) => {
    setActiveCalculationId(id);
  }, []);

  const clearActiveCalculation = useCallback(() => {
    setActiveCalculationId(null);
  }, []);

  // Sync active ID if it's no longer in the main list
  useEffect(() => {
    if (activeCalculationId && calculations.length > 0 && !calculations.some(calc => calc.id === activeCalculationId)) {
      setActiveCalculationId(calculations[0].id);
    } else if (activeCalculationId && calculations.length === 0) {
        setActiveCalculationId(null);
    }
  }, [activeCalculationId, calculations]);

  // Fetch details when active ID changes
  useEffect(() => {
    if (activeCalculationId && !detailedCalculations.has(activeCalculationId)) {
      loadCalculationDetails(activeCalculationId);
    }
  }, [activeCalculationId, detailedCalculations, loadCalculationDetails]);
  
  return {
    activeCalculationId,
    activeCalculation,
    isLoadingDetails,
    detailsError,
    setActiveCalculation,
    setActiveCalculationById,
    loadCalculationDetails,
    clearActiveCalculation
  };
};