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

      // 計算完了通知: running -> completed
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

      // 計算エラー通知: running -> error または pending -> error
      if (
        (previousStatus === 'running' || previousStatus === 'pending') &&
        updatedCalculation.status === 'error'
      ) {
        const errorMessage =
          currentErrorMessage || 'Detailed error information is not available.';

        const isResourceError = isResourceInsufficientError(errorMessage);

        if (isResourceError) {
          handleError(new Error(errorMessage), 'Calculation failed');
        } else {
          handleError(
            new Error(errorMessage),
            `Calculation "${molecularName}" failed`
          );
        }
      }

      // 計算開始通知: pending -> running
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

      // 待機からの開始通知: waiting -> running
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

      // 待機状態への移行通知: pending -> waiting
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

  // 30秒以内の重複通知を防止
  const handleWebSocketError = useCallback((error: string) => {
    console.error('[UnifiedWebSocket] Error:', error);

    const now = Date.now();
    const timeSinceLastNotification =
      now - lastErrorNotificationTimeRef.current;

    if (timeSinceLastNotification > 30000) {
      handleError(error, 'Real-time monitoring error');
      lastErrorNotificationTimeRef.current = now;
    } else {
      console.log(
        `[UnifiedWebSocket] Skipping duplicate error notification (${Math.floor(timeSinceLastNotification / 1000)}s since last)`
      );
    }
  }, []);

  return { notifyCalculationUpdate, handleWebSocketError };
};
