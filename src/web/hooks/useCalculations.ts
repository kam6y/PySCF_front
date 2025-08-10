// src/web/hooks/useCalculations.ts

import { useState, useEffect, useCallback } from 'react';
import {
    CalculationInstance,
    CalculationParameters,
} from '../types/calculation';
import { getCalculations } from '../apiClient'; // apiClientからインポート

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
        const dirName = backendCalc.name;
        
        return {
            id: dirName,
            name: dirName.replace(/_\d{8}_\d{6}$/, ''),
            status: backendCalc.status || (backendCalc.has_checkpoint ? 'completed' : 'pending'),
            createdAt: new Date(backendCalc.date).toISOString(),
            updatedAt: new Date(backendCalc.date).toISOString(),
            parameters: {} as CalculationParameters, // 詳細は別途読み込み
            results: undefined
        };
    };

    const fetchCalculations = useCallback(async () => {
        setIsLoading(true);
        setError(null);
        try {
            const data = await getCalculations();
            const convertedCalculations = data.data.calculations.map(convertToCalculationInstance);
            setCalculations(convertedCalculations);
        } catch (err) {
            const errorMessage = err instanceof Error ? err.message : '不明なエラーが発生しました。';
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
        // TODO: バックエンドAPI呼び出しを実装
        console.warn("updateCalculationNameはバックエンドと永続化されません。");
        setCalculations(prev =>
            prev.map(calc =>
                calc.id === id ? { ...calc, name: newName, updatedAt: new Date().toISOString() } : calc
            )
        );
    }, []);

    const updateCalculationStatus = useCallback(async (id: string, status: 'pending' | 'running' | 'completed' | 'error') => {
        // TODO: バックエンドAPI呼び出しを実装
        console.warn("updateCalculationStatusはバックエンドと永続化されません。");
        setCalculations(prev =>
            prev.map(calc =>
                calc.id === id ? { ...calc, status, updatedAt: new Date().toISOString() } : calc
            )
        );
    }, []);

    const deleteCalculation = useCallback(async (id: string) => {
        // TODO: バックエンドAPI呼び出しを実装
        console.warn("deleteCalculationはバックエンドと永続化されません。");
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
        setCalculations(prevCalculations => {
            const exists = prevCalculations.some(c => c.id === updatedCalculation.id);
            if (exists) {
                return prevCalculations.map(calc =>
                    calc.id === updatedCalculation.id ? updatedCalculation : calc
                );
            }
            // 新規計算（'pending'から'completed'になった場合など）の処理
            return [updatedCalculation, ...prevCalculations.filter(c => !c.id.startsWith('new-calculation-'))];
        });
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