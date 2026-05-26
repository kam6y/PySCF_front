import { useCallback, useRef } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import { CalculationInstance } from '../types/api-types';
import { handleError } from '../utils/errorHandler';
import { calculationQueryKeys } from '../hooks/useCalculationQueries';
import { invalidateQueriesWithRetry } from './invalidateQueriesWithRetry';

export interface UseCalculationSyncOptions {
  activeCalculationId: string | null;
  onCalculationUpdate?: (
    updated: CalculationInstance,
    previousStatus: string | undefined
  ) => void;
  onStreamError?: (error: string) => void;
}

export interface UseCalculationSyncReturn {
  handleCalculationUpdate: (updated: CalculationInstance) => void;
  handleStreamError: (errorData: unknown) => void;
  handleReconnect: () => Promise<void>;
}

const resolveStreamErrorMessage = (errorData: unknown): string => {
  if (typeof errorData === 'string' && errorData.length > 0) {
    return errorData;
  }

  if (typeof errorData === 'object' && errorData !== null) {
    const record = errorData as Record<string, unknown>;
    const maybeError = record.error;
    if (typeof maybeError === 'string' && maybeError.length > 0) {
      return maybeError;
    }

    const maybeMessage = record.message;
    if (typeof maybeMessage === 'string' && maybeMessage.length > 0) {
      return maybeMessage;
    }

    const maybePayload = record.payload;
    if (typeof maybePayload === 'object' && maybePayload !== null) {
      const maybePayloadMessage = (maybePayload as Record<string, unknown>)
        .message;
      if (
        typeof maybePayloadMessage === 'string' &&
        maybePayloadMessage.length > 0
      ) {
        return maybePayloadMessage;
      }
    }
  }

  return 'A problem occurred with calculation monitoring.';
};

export const useCalculationSync = ({
  activeCalculationId,
  onCalculationUpdate,
  onStreamError,
}: UseCalculationSyncOptions): UseCalculationSyncReturn => {
  const queryClient = useQueryClient();
  const activeCalculationIdRef = useRef(activeCalculationId);
  activeCalculationIdRef.current = activeCalculationId;
  const lastNotifiedStatus = useRef<Map<string, string>>(new Map());
  const lastNotifiedUpdatedAt = useRef<Map<string, string>>(new Map());
  const syncFailureCountRef = useRef<number>(0);

  const onCalculationUpdateRef = useRef(onCalculationUpdate);
  onCalculationUpdateRef.current = onCalculationUpdate;
  const onStreamErrorRef = useRef(onStreamError);
  onStreamErrorRef.current = onStreamError;

  const updateCalculationCache = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;

      try {
        queryClient.invalidateQueries({
          queryKey: calculationQueryKeys.detail(calculationId),
        });

        queryClient.invalidateQueries({
          queryKey: calculationQueryKeys.list(),
        });

        console.log(
          `[CalculationUpdates] Invalidated queries for calculation ${calculationId}`
        );
      } catch (error) {
        console.error(
          `[CalculationUpdates] Failed to invalidate queries for calculation ${calculationId}:`,
          error
        );
      }
    },
    [queryClient]
  );

  const notifyStreamError = useCallback(
    (errorMessage: string) => {
      queryClient.invalidateQueries({
        queryKey: calculationQueryKeys.list(),
      });
      onStreamErrorRef.current?.(errorMessage);
    },
    [queryClient]
  );

  const handleCalculationUpdate = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;

      if (!updatedCalculation || !calculationId) {
        console.warn(
          '[CalculationUpdates] Received invalid calculation update data:',
          updatedCalculation
        );
        return;
      }

      console.log(
        `[CalculationUpdates] Processing update for calculation ${calculationId}: ${updatedCalculation.status}`
      );

      const previousStatus = lastNotifiedStatus.current.get(calculationId);
      const previousUpdatedAt =
        lastNotifiedUpdatedAt.current.get(calculationId);

      if (previousUpdatedAt === updatedCalculation.updatedAt) {
        console.log(
          `[CalculationUpdates] Skipping duplicate update for ${calculationId}: updatedAt unchanged (${updatedCalculation.updatedAt})`
        );
        return;
      }

      lastNotifiedStatus.current.set(calculationId, updatedCalculation.status);
      lastNotifiedUpdatedAt.current.set(
        calculationId,
        updatedCalculation.updatedAt
      );

      console.log(
        `[CalculationUpdates] Update detected for ${calculationId}: status=${previousStatus ?? 'none'} -> ${updatedCalculation.status}, updatedAt=${updatedCalculation.updatedAt}`
      );

      updateCalculationCache(updatedCalculation);
      onCalculationUpdateRef.current?.(updatedCalculation, previousStatus);
    },
    [updateCalculationCache]
  );

  const handleStreamError = useCallback(
    (errorData: unknown) => {
      console.error('[CalculationUpdates] Stream error:', errorData);

      const errorMessage = resolveStreamErrorMessage(errorData);

      if (errorMessage.includes('not found')) {
        console.log(
          '[CalculationUpdates] Calculation not found (likely directory changed), ignoring error'
        );
        return;
      }

      notifyStreamError(errorMessage);
    },
    [notifyStreamError]
  );

  const handleReconnect = useCallback(async () => {
    console.log('[CalculationUpdates] Reconnected calculation update stream');
    console.log('[CalculationUpdates] Syncing data after reconnection...');

    const activeId = activeCalculationIdRef.current;

    try {
      await Promise.all([
        invalidateQueriesWithRetry({
          queryClient,
          queryKey: [...calculationQueryKeys.list()],
        }),
        activeId && !activeId.startsWith('new-calculation-')
          ? invalidateQueriesWithRetry({
              queryClient,
              queryKey: [...calculationQueryKeys.detail(activeId)],
            })
          : Promise.resolve(),
      ]);

      syncFailureCountRef.current = 0;
      console.log('[CalculationUpdates] Data sync completed successfully');
    } catch (error) {
      syncFailureCountRef.current++;
      console.error(
        `[CalculationUpdates] Data sync failed (attempt ${syncFailureCountRef.current}/2):`,
        error
      );

      if (syncFailureCountRef.current >= 2) {
        handleError(
          new Error(
            'Unable to sync calculation data after reconnection. Please refresh the page if calculations appear outdated.'
          ),
          'Connection synchronization failed'
        );
        syncFailureCountRef.current = 0;
      }
    }
  }, [queryClient]);

  return {
    handleCalculationUpdate,
    handleStreamError,
    handleReconnect,
  };
};
