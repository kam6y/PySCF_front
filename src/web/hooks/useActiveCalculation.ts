import { useState, useEffect, useCallback } from 'react';
import { CalculationInstance, CalculationDetailsResponse } from '../types/calculation';

const ACTIVE_CALCULATION_KEY = 'pyscf-active-calculation-id';
const API_BASE_URL = 'http://127.0.0.1:5000';

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
    // Restore active calculation ID from localStorage
    try {
      return localStorage.getItem(ACTIVE_CALCULATION_KEY);
    } catch {
      return null;
    }
  });
  
  const [detailedCalculations, setDetailedCalculations] = useState<Map<string, CalculationInstance>>(new Map());
  const [isLoadingDetails, setIsLoadingDetails] = useState(false);
  const [detailsError, setDetailsError] = useState<string | null>(null);

  // Get active calculation from the calculations array or detailed cache
  const activeCalculation = activeCalculationId 
    ? detailedCalculations.get(activeCalculationId) || calculations.find(calc => calc.id === activeCalculationId)
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
  const loadCalculationDetails = useCallback(async (id: string, forceRefresh = false): Promise<CalculationInstance | null> => {
    // Check cache first unless forced refresh
    if (!forceRefresh && detailedCalculations.has(id)) {
      return detailedCalculations.get(id) || null;
    }
    
    setIsLoadingDetails(true);
    setDetailsError(null);
    
    try {
      const response = await fetch(`${API_BASE_URL}/api/quantum/calculations/${id}`);
      
      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || 'Failed to load calculation details');
      }
      
      const data: CalculationDetailsResponse = await response.json();
      
      if (!data.success) {
        throw new Error(data.error || 'Failed to load calculation details');
      }
      
      const detailedCalculation = data.data.calculation;
      
      // Cache the detailed calculation
      setDetailedCalculations(prev => {
        const newMap = new Map(prev);
        newMap.set(id, detailedCalculation);
        return newMap;
      });
      
      return detailedCalculation;
      
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Unknown error occurred';
      setDetailsError(errorMessage);
      console.error('Failed to load calculation details:', err);
      return null;
    } finally {
      setIsLoadingDetails(false);
    }
  }, [detailedCalculations]);

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
  
  // Auto-load details when active calculation changes
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