// src/web/hooks/useCalculations.ts

import { useState, useEffect, useCallback } from 'react';
import {
    CalculationInstance,
    CalculationParameters,
    CalculationListResponse,
    RenameResponse,
} from '../types/calculation';
import * as apiClient from '../apiClient';

export interface UseCalculationsReturn {
    calculations: CalculationInstance[];
    isLoading: boolean;
    error: string | null;
    refreshCalculations: () => Promise<void>;
    getCalculationById: (id: string) => CalculationInstance | undefined;
    updateCalculationName: (id: string, newName: string) => Promise<RenameResponse>;
    deleteCalculation: (id: string) => Promise<void>;
    createCalculation: (name: string, parameters: CalculationParameters) => CalculationInstance;
    updateCalculation: (updatedCalculation: CalculationInstance) => void;
    createNewCalculationFromExisting: (originalCalc: CalculationInstance, newParams: CalculationParameters) => CalculationInstance;
}

export const useCalculations = (): UseCalculationsReturn => {
    const [calculations, setCalculations] = useState<CalculationInstance[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    const convertToCalculationInstance = (backendCalc: CalculationListResponse['data']['calculations'][0]): CalculationInstance => {
        return {
            id: backendCalc.id,
            name: backendCalc.name,
            status: backendCalc.status,
            createdAt: new Date(backendCalc.date).toISOString(),
            updatedAt: new Date(backendCalc.date).toISOString(),
            parameters: {} as CalculationParameters, // Details are loaded separately
            results: undefined
        };
    };

    const fetchCalculations = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const data = await apiClient.getCalculations();
            const convertedCalculations = data.calculations.map(convertToCalculationInstance);
            setCalculations(convertedCalculations);
        } catch (err) {
            const errorMessage = err instanceof Error ? err.message : 'An unknown error occurred.';
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
        try {
            const response = await apiClient.updateCalculationName(id, newName);
            await refreshCalculations();
            return response;
        } catch (error) {
            console.error("Failed to update calculation name:", error);
            throw error;
        }
    }, [refreshCalculations]);

    const deleteCalculation = useCallback(async (id: string) => {
        try {
            await apiClient.deleteCalculation(id);
            setCalculations(prev => prev.filter(calc => calc.id !== id));
        } catch (error) {
            console.error("Failed to delete calculation:", error);
            await refreshCalculations();
            throw error;
        }
    }, [refreshCalculations]);

    const createCalculation = useCallback((name: string, parameters: CalculationParameters): CalculationInstance => {
        const now = new Date().toISOString();
        const id = `new-calculation-${Date.now()}`;
        
        const newCalculation: CalculationInstance = {
            id,
            name: name || "", // no default name (must be set by user)
            status: 'pending',
            createdAt: now,
            updatedAt: now,
            parameters,
            results: undefined
        };
        
        setCalculations(prev => [newCalculation, ...prev]);
        return newCalculation;
    }, []);

    const createNewCalculationFromExisting = useCallback((originalCalc: CalculationInstance, newParams: CalculationParameters): CalculationInstance => {
        const now = new Date().toISOString();
        const id = `new-calculation-${Date.now()}`;

        const newCalculation: CalculationInstance = {
            ...originalCalc,
            id,
            status: 'pending',
            createdAt: now,
            updatedAt: now,
            parameters: newParams,
            results: undefined,
            workingDirectory: undefined,
            errorMessage: undefined,
        };

        setCalculations(prev => [newCalculation, ...prev]);
        return newCalculation;
    }, []);

    const updateCalculation = useCallback((updatedCalculation: CalculationInstance) => {
        setCalculations(prevCalculations => {
            // First, check if this is a server response meant to replace a temporary calculation.
            // This happens when a calculation is started.
            if (!updatedCalculation.id.startsWith('new-calculation-')) {
                const tempIndex = prevCalculations.findIndex(c => c.id.startsWith('new-calculation-'));
                if (tempIndex !== -1) {
                    const newCalculations = [...prevCalculations];
                    newCalculations[tempIndex] = updatedCalculation;
                    return newCalculations.sort((a, b) => new Date(b.updatedAt).getTime() - new Date(a.updatedAt).getTime());
                }
            }

            // If not replacing a temp calc, check if it's an update to an existing one.
            // This handles status updates during polling or parameter changes.
            const existingIndex = prevCalculations.findIndex(c => c.id === updatedCalculation.id);
            if (existingIndex !== -1) {
                const newCalculations = [...prevCalculations];
                newCalculations[existingIndex] = updatedCalculation;
                return newCalculations;
            }

            // If it's not an update and not a replacement, it must be a new (temporary) calculation.
            return [updatedCalculation, ...prevCalculations];
        });
    }, []);

    return {
        calculations,
        isLoading,
        error,
        refreshCalculations,
        getCalculationById,
        updateCalculationName,
        deleteCalculation,
        createCalculation,
        updateCalculation,
        createNewCalculationFromExisting,
    };
};