// src/web/hooks/useActiveCalculation.ts

import { useState, useEffect, useCallback } from 'react';
import { CalculationInstance } from '../types/calculation';
import { getCalculationDetails } from '../apiClient';

const ACTIVE_CALCULATION_KEY = 'pyscf-active-calculation-id';

export interface UseActiveCalculationReturn {
  activeCalculationId: string | null;
  activeCalculation: CalculationInstance | undefined;
  isLoadingDetails: boolean;
  detailsError: string | null;
  setActiveCalculationById: (id: string | null) => void;
  loadCalculationDetails: (id: string, forceRefresh?: boolean) => Promise<CalculationInstance | null>;
  clearActiveCalculation: () => void;
  updateActiveCalculationInCache: (updatedCalculation: CalculationInstance) => void;
  renameIdInCache: (oldId: string, newId: string, newName: string) => void;
}

export const useActiveCalculation = (
  calculations: CalculationInstance[]
): UseActiveCalculationReturn => {
  const [activeCalculationId, setActiveCalculationId] = useState<string | null>(() => {
    try {
      const storedId = localStorage.getItem(ACTIVE_CALCULATION_KEY);
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
  
  const setActiveCalculationById = useCallback((id: string | null) => {
    setActiveCalculationId(id);
  }, []);

  const clearActiveCalculation = useCallback(() => {
    setActiveCalculationId(null);
  }, []);

  const updateActiveCalculationInCache = useCallback((updatedCalculation: CalculationInstance) => {
    setDetailedCalculations(prev => {
        const newMap = new Map(prev);
        newMap.set(updatedCalculation.id, updatedCalculation);
        return newMap;
    });
  }, []);

  const renameIdInCache = useCallback((oldId: string, newId: string, newName: string) => {
    setDetailedCalculations(prev => {
        const newMap = new Map(prev);
        const data = newMap.get(oldId);
        if (data) {
            newMap.delete(oldId);
            newMap.set(newId, { 
                ...data, 
                id: newId, 
                name: newName,
                parameters: { ...data.parameters, molecule_name: newName }
            });
        }
        return newMap;
    });
  }, []);

  useEffect(() => {
    if (activeCalculationId && calculations.length > 0 && !calculations.some(calc => calc.id === activeCalculationId)) {
      setActiveCalculationId(calculations[0].id);
    } else if (activeCalculationId && calculations.length === 0) {
        setActiveCalculationId(null);
    }
  }, [activeCalculationId, calculations]);

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
    setActiveCalculationById,
    loadCalculationDetails,
    clearActiveCalculation,
    updateActiveCalculationInCache,
    renameIdInCache
  };
};