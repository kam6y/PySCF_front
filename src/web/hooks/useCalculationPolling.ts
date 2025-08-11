// src/web/hooks/useCalculationPolling.ts

import { useRef, useCallback, useEffect } from 'react';
import { CalculationInstance } from '../types/api-types';
import { getCalculationDetails } from '../apiClient';

export interface UseCalculationPollingOptions {
  calculationId: string | null;
  onUpdate: (calculation: CalculationInstance) => void;
  onError?: (error: string) => void;
}

export interface UseCalculationPollingReturn {
  isPolling: boolean;
  startPolling: () => void;
  stopPolling: () => void;
}

export const useCalculationPolling = ({
  calculationId,
  onUpdate,
  onError
}: UseCalculationPollingOptions): UseCalculationPollingReturn => {
  const pollingIntervalRef = useRef<NodeJS.Timeout | null>(null);

  const stopPolling = useCallback(() => {
    if (pollingIntervalRef.current) {
      clearInterval(pollingIntervalRef.current);
      pollingIntervalRef.current = null;
    }
  }, []);

  const startPolling = useCallback(() => {
    if (!calculationId) return;
    
    stopPolling(); // To prevent multiple pollers from running

    pollingIntervalRef.current = setInterval(async () => {
      try {
        const response = await getCalculationDetails(calculationId);
        const updatedCalc = response.calculation;

        if (updatedCalc.status === 'completed' || updatedCalc.status === 'error') {
          stopPolling();
          onUpdate(updatedCalc);
        } else {
          onUpdate(updatedCalc);
        }
      } catch (error) {
        console.error(`Polling failed for calculation ${calculationId}:`, error);
        const errorMessage = error instanceof Error ? error.message : 'Unknown error occurred';
        if (onError) {
          onError(`Failed to get calculation status. ${errorMessage}`);
        }
        stopPolling();
      }
    }, 3000); // Poll every 3 seconds
  }, [calculationId, onUpdate, onError, stopPolling]);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      stopPolling();
    };
  }, [stopPolling]);

  const isPolling = pollingIntervalRef.current !== null;

  return {
    isPolling,
    startPolling,
    stopPolling,
  };
};