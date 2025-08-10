import { useState, useEffect, useCallback } from 'react';
import { 
    CalculationInstance, 
    CalculationListResponse,
    CalculationParameters,
} from '../types/calculation';

const API_BASE_URL = 'http://127.0.0.1:5000';

export interface UseCalculationsReturn {
    calculations: CalculationInstance[];
    isLoading: boolean;
    error: string | null;
    refreshCalculations: () => Promise<void>;
    getCalculationById: (id: string) => CalculationInstance | undefined;
    updateCalculationName: (id: string, newName: string) => Promise<void>;
    updateCalculationStatus: (id: string, status: 'pending' | 'running' | 'completed' | 'error') => Promise<void>;
    deleteCalculation: (id: string) => Promise<void>;
    createCalculation: (name: string, parameters: CalculationParameters) => CalculationInstance;
    updateCalculation: (updatedCalculation: CalculationInstance) => void;
}

export const useCalculations = (): UseCalculationsReturn => {
    const [calculations, setCalculations] = useState<CalculationInstance[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    const convertToCalculationInstance = (backendCalc: any): CalculationInstance => {
        const pathParts = backendCalc.path.split(/[\\/]/);
        const dirName = pathParts[pathParts.length - 1];
        
        return {
            id: dirName,
            name: dirName.replace(/_\d{8}_\d{6}$/, ''),
            status: backendCalc.status || (backendCalc.has_checkpoint ? 'completed' : 'pending'),
            createdAt: new Date(backendCalc.date).toISOString(),
            updatedAt: new Date(backendCalc.date).toISOString(),
            parameters: {} as CalculationParameters,
            results: undefined
        };
    };

    const fetchCalculations = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const response = await fetch(`${API_BASE_URL}/api/quantum/calculations`);
            if (!response.ok) throw new Error('Network response was not ok');
            const data: CalculationListResponse = await response.json();
            if (!data.success) throw new Error(data.error || 'Failed to fetch calculations');
            
            const convertedCalculations = data.data.calculations.map(convertToCalculationInstance);
            setCalculations(convertedCalculations);
        } catch (err) {
            const errorMessage = err instanceof Error ? err.message : 'Unknown error occurred';
            setError(errorMessage);
        } finally {
            setIsLoading(false);
        }
    }, []);

    useEffect(() => {
        fetchCalculations();
    }, [fetchCalculations]);

    const refreshCalculations = useCallback(async () => {
        await fetchCalculations();
    }, [fetchCalculations]);

    const getCalculationById = useCallback((id: string): CalculationInstance | undefined => {
        return calculations.find(calc => calc.id === id);
    }, [calculations]);

    const updateCalculationName = useCallback(async (id: string, newName: string) => {
        // (Backend API call logic here)
        setCalculations(prev => 
            prev.map(calc => 
                calc.id === id ? { ...calc, name: newName, updatedAt: new Date().toISOString() } : calc
            )
        );
    }, []);

    const updateCalculationStatus = useCallback(async (id: string, status: 'pending' | 'running' | 'completed' | 'error') => {
        // (Backend API call logic here)
        setCalculations(prev => 
            prev.map(calc => 
                calc.id === id ? { ...calc, status, updatedAt: new Date().toISOString() } : calc
            )
        );
    }, []);

    const deleteCalculation = useCallback(async (id: string) => {
        // (Backend API call logic here)
        setCalculations(prev => prev.filter(calc => calc.id !== id));
    }, []);

    const createCalculation = useCallback((name: string, parameters: CalculationParameters): CalculationInstance => {
        const now = new Date().toISOString();
        const id = `new-calculation-${Date.now()}`;
        
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

    const updateCalculation = useCallback((updatedCalculation: CalculationInstance) => {
        setCalculations(prevCalculations => 
            prevCalculations.map(calc => 
                calc.id === updatedCalculation.id ? updatedCalculation : calc
            )
        );
    }, []);

    return {
        calculations,
        isLoading,
        error,
        refreshCalculations,
        getCalculationById,
        updateCalculationName,
        updateCalculationStatus,
        deleteCalculation,
        createCalculation,
        updateCalculation,
    };
};