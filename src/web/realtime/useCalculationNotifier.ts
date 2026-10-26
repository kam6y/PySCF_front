import { useCallback, useRef } from 'react';
import { CalculationInstance } from '../types/api-types';
import {
  showSuccessNotification,
  showInfoNotification,
} from '../store/notificationStore';
import { handleError } from '../utils/errorHandler';
import { isResourceInsufficientError } from '../utils/errorClassifier';

export const useCalculationNotifier = (activeCalculationId: string | null) => {
  const lastErrorNotificationTimeRef = useRef<number>(0);

  const notifyCalculationUpdate = useCallback(
    (
      updatedCalculation: CalculationInstance,
      previousStatus: string | undefined
    ) => {
      const calculationId = updatedCalculation.id;
      const isActiveCalculation = calculationId === activeCalculationId;
      const molecularName = updatedCalculation.name || 'Unknown';
      const currentErrorMessage =
        updatedCalculation.error || updatedCalculation.results?.error;

      if (
        previousStatus === 'running' &&
        updatedCalculation.status === 'completed'
      ) {
        if (!isActiveCalculation) {
          showSuccessNotification(
            'Calculation completed',
            `Calculation for ${molecularName} has been completed.`,
            calculationId
          );
        }
      }

      if (
        (previousStatus === 'running' || previousStatus === 'pending') &&
        updatedCalculation.status === 'error'
      ) {
        const errorMessage =
          currentErrorMessage || 'Detailed error information is not available.';

        if (isResourceInsufficientError(errorMessage)) {
          handleError(new Error(errorMessage), 'Calculation failed');
        } else {
          handleError(
            new Error(errorMessage),
            `Calculation "${molecularName}" failed`
          );
        }
      }

      if (
        previousStatus === 'pending' &&
        updatedCalculation.status === 'running'
      ) {
        if (!isActiveCalculation) {
          showInfoNotification(
            'Calculation started',
            `Calculation for ${molecularName} has been started.`,
            calculationId
          );
        }
      }

      if (
        previousStatus === 'waiting' &&
        updatedCalculation.status === 'running'
      ) {
        if (!isActiveCalculation) {
          showInfoNotification(
            'Queued calculation started',
            `Calculation for ${molecularName} has started from queued status.`,
            calculationId
          );
        }
      }

      if (
        previousStatus === 'pending' &&
        updatedCalculation.status === 'waiting'
      ) {
        if (!isActiveCalculation) {
          const waitingReason =
            updatedCalculation.waitingReason ||
            'Please wait for available system resources or execution slots.';
          showInfoNotification(
            'Calculation is waiting',
            `Calculation for ${molecularName} is currently waiting. Reason: ${waitingReason}`,
            calculationId
          );
        }
      }
    },
    [activeCalculationId]
  );

  const handleStreamError = useCallback((error: string) => {
    console.error('[CalculationUpdates] Error:', error);

    const now = Date.now();
    const timeSinceLastNotification =
      now - lastErrorNotificationTimeRef.current;

    if (timeSinceLastNotification > 30000) {
      handleError(error, 'Real-time monitoring error');
      lastErrorNotificationTimeRef.current = now;
      return;
    }

    console.log(
      `[CalculationUpdates] Skipping duplicate error notification (${Math.floor(timeSinceLastNotification / 1000)}s since last)`
    );
  }, []);

  return { notifyCalculationUpdate, handleStreamError };
};
