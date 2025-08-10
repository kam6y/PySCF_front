import { useState, useEffect, useCallback } from 'react';
import { 
  CalculationInstance, 
  CalculationListResponse,
  CalculationDetailsResponse,
  CalculationParameters,
  CalculationResults 
} from '../types/calculation';

const API_BASE_URL = 'http://127.0.0.1:5000';

export interface UseCalculationsReturn {
  calculations: CalculationInstance[];
  isLoading: boolean;
  error: string | null;
  refreshCalculations: () => Promise<void>;
  getCalculationById: (id: string) => CalculationInstance | undefined;
  updateCalculationName: (id: string, newName: string) => Promise<void>;
  deleteCalculation: (id: string) => Promise<void>;
  createCalculation: (name: string, parameters: CalculationParameters) => CalculationInstance;
}

export const useCalculations = (): UseCalculationsReturn => {
  const [calculations, setCalculations] = useState<CalculationInstance[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Helper function to convert backend calculation data to CalculationInstance
  const convertToCalculationInstance = (backendCalc: any): CalculationInstance => {
    // Extract calculation ID from directory name (e.g., "water_test_20250810_123019")
    const pathParts = backendCalc.path.split('/');
    const dirName = pathParts[pathParts.length - 1];
    
    return {
      id: dirName,
      name: dirName.replace(/_\d{8}_\d{6}$/, ''), // Remove timestamp suffix for display
      status: backendCalc.has_checkpoint ? 'completed' : 'pending',
      createdAt: new Date(backendCalc.date).toISOString(),
      updatedAt: new Date(backendCalc.date).toISOString(),
      workingDirectory: backendCalc.path,
      parameters: {} as CalculationParameters, // Will be populated when needed
      results: undefined // Will be populated when needed
    };
  };

  // Fetch calculations from backend
  const fetchCalculations = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    
    try {
      const response = await fetch(`${API_BASE_URL}/api/quantum/calculations`);
      const data: CalculationListResponse = await response.json();
      
      if (!data.success) {
        throw new Error(data.error || 'Failed to fetch calculations');
      }
      
      const convertedCalculations = data.data.calculations.map(convertToCalculationInstance);
      setCalculations(convertedCalculations);
      
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Unknown error occurred';
      setError(errorMessage);
      console.error('Failed to fetch calculations:', err);
    } finally {
      setIsLoading(false);
    }
  }, []);

  // Refresh calculations
  const refreshCalculations = useCallback(async () => {
    await fetchCalculations();
  }, [fetchCalculations]);

  // Get calculation by ID
  const getCalculationById = useCallback((id: string): CalculationInstance | undefined => {
    return calculations.find(calc => calc.id === id);
  }, [calculations]);

  // Update calculation name
  const updateCalculationName = useCallback(async (id: string, newName: string) => {
    try {
      // For now, just update locally since backend doesn't support name updates yet
      setCalculations(prev => 
        prev.map(calc => 
          calc.id === id 
            ? { ...calc, name: newName, updatedAt: new Date().toISOString() }
            : calc
        )
      );
      
      // TODO: Implement backend API call when available
      // const response = await fetch(`${API_BASE_URL}/api/quantum/calculations/${id}`, {
      //   method: 'PUT',
      //   headers: { 'Content-Type': 'application/json' },
      //   body: JSON.stringify({ name: newName })
      // });
      
    } catch (err) {
      console.error('Failed to update calculation name:', err);
      throw err;
    }
  }, []);

  // Delete calculation
  const deleteCalculation = useCallback(async (id: string) => {
    try {
      // For now, just remove locally since backend doesn't support deletion yet
      setCalculations(prev => prev.filter(calc => calc.id !== id));
      
      // TODO: Implement backend API call when available
      // const response = await fetch(`${API_BASE_URL}/api/quantum/calculations/${id}`, {
      //   method: 'DELETE'
      // });
      
    } catch (err) {
      console.error('Failed to delete calculation:', err);
      throw err;
    }
  }, []);

  // Create new calculation instance (local only, not persisted until calculation starts)
  const createCalculation = useCallback((name: string, parameters: CalculationParameters): CalculationInstance => {
    const now = new Date().toISOString();
    const timestamp = new Date().toISOString().replace(/[:.]/g, '').slice(0, 15);
    const id = `${name}_${timestamp}`;
    
    const newCalculation: CalculationInstance = {
      id,
      name,
      status: 'pending',
      createdAt: now,
      updatedAt: now,
      parameters,
      results: undefined
    };
    
    setCalculations(prev => [newCalculation, ...prev]);
    return newCalculation;
  }, []);

  // Load calculations on component mount
  useEffect(() => {
    fetchCalculations();
  }, [fetchCalculations]);

  return {
    calculations,
    isLoading,
    error,
    refreshCalculations,
    getCalculationById,
    updateCalculationName,
    deleteCalculation,
    createCalculation
  };
};