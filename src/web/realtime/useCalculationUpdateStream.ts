import { useCallback, useEffect, useRef, useState } from 'react';
import type { MutableRefObject } from 'react';
import { fetchEventSource } from '@microsoft/fetch-event-source';
import { getApiBaseUrl } from '../api/core';
import type { CalculationInstance } from '../types/api-types';

const STREAM_RETRY_INTERVAL_MS = 2000;
const GLOBAL_CALCULATION_UPDATES_PATH =
  '/api/quantum/calculations/updates/stream';

type CalculationUpdateStreamEventType =
  | 'calculation_update'
  | 'heartbeat'
  | 'error';

export interface ParsedCalculationUpdateStreamEvent {
  type: CalculationUpdateStreamEventType;
  calculation?: CalculationInstance;
  errorMessage?: string;
}

export interface UseCalculationUpdateStreamOptions {
  activeCalculationId: string | null;
  onCalculationUpdate?: (updated: CalculationInstance) => void;
  onStreamError?: (error: string) => void;
  onReconnect?: () => void | Promise<void>;
}

export interface UseCalculationUpdateStreamReturn {
  isConnected: boolean;
  reconnect: () => void;
  disconnect: () => void;
}

const isUsableCalculationId = (
  calculationId: string | null
): calculationId is string => {
  return !!calculationId && !calculationId.startsWith('new-calculation-');
};

const isAbortError = (error: unknown, signal: AbortSignal): boolean => {
  if (signal.aborted) {
    return true;
  }

  if (typeof error !== 'object' || error === null || !('name' in error)) {
    return false;
  }

  return error.name === 'AbortError';
};

const toStreamErrorMessage = (error: unknown): string => {
  if (error instanceof Error) {
    return error.message;
  }

  if (typeof error === 'string' && error.length > 0) {
    return error;
  }

  return 'Calculation update stream failed.';
};

export const parseCalculationUpdateStreamEvent = (
  data: string
): ParsedCalculationUpdateStreamEvent => {
  const parsed = JSON.parse(data) as {
    type?: unknown;
    payload?: Record<string, unknown>;
  };

  if (parsed.type === 'calculation_update') {
    return {
      type: 'calculation_update',
      calculation: parsed.payload?.calculation as
        | CalculationInstance
        | undefined,
    };
  }

  if (parsed.type === 'error') {
    const message = parsed.payload?.message;
    return {
      type: 'error',
      errorMessage:
        typeof message === 'string'
          ? message
          : 'Calculation update stream reported an error.',
    };
  }

  if (parsed.type === 'heartbeat' || parsed.type === undefined) {
    return { type: 'heartbeat' };
  }

  return { type: 'heartbeat' };
};

export const useCalculationUpdateStream = ({
  activeCalculationId,
  onCalculationUpdate,
  onStreamError,
  onReconnect,
}: UseCalculationUpdateStreamOptions): UseCalculationUpdateStreamReturn => {
  const globalAbortRef = useRef<AbortController | null>(null);
  const detailAbortRef = useRef<AbortController | null>(null);
  const hasGlobalConnectedRef = useRef<boolean>(false);
  const isMountedRef = useRef<boolean>(false);
  const activeCalculationIdRef = useRef(activeCalculationId);
  activeCalculationIdRef.current = activeCalculationId;
  const [isConnected, setIsConnected] = useState(false);

  const onCalculationUpdateRef = useRef(onCalculationUpdate);
  onCalculationUpdateRef.current = onCalculationUpdate;
  const onStreamErrorRef = useRef(onStreamError);
  onStreamErrorRef.current = onStreamError;
  const onReconnectRef = useRef(onReconnect);
  onReconnectRef.current = onReconnect;

  const setGlobalConnected = useCallback((connected: boolean) => {
    if (isMountedRef.current) {
      setIsConnected(connected);
    }
  }, []);

  const notifyStreamError = useCallback((error: unknown) => {
    onStreamErrorRef.current?.(toStreamErrorMessage(error));
  }, []);

  const startStream = useCallback(
    (
      path: string,
      abortRef: MutableRefObject<AbortController | null>,
      trackGlobalConnection = false
    ) => {
      abortRef.current?.abort();

      const controller = new AbortController();
      abortRef.current = controller;

      (async () => {
        try {
          const authToken = await window.electronAPI?.getAuthToken?.();

          if (controller.signal.aborted) {
            return;
          }

          const headers: HeadersInit = {
            Accept: 'text/event-stream',
          };

          if (authToken) {
            (headers as Record<string, string>)['X-Auth-Token'] = authToken;
          }

          await fetchEventSource(`${getApiBaseUrl()}${path}`, {
            method: 'GET',
            headers,
            signal: controller.signal,

            onopen: async response => {
              if (!response.ok) {
                const errorText = await response
                  .text()
                  .catch(() => response.statusText);
                const error = new Error(
                  `Failed to open calculation update stream: ${response.status} ${errorText}`
                );
                if (trackGlobalConnection) {
                  setGlobalConnected(false);
                }
                notifyStreamError(error);
                controller.abort();
                throw error;
              }

              if (!trackGlobalConnection) {
                return;
              }

              setGlobalConnected(true);

              if (!hasGlobalConnectedRef.current) {
                hasGlobalConnectedRef.current = true;
                return;
              }

              Promise.resolve(onReconnectRef.current?.()).catch(error => {
                console.error(
                  '[CalculationUpdates] Reconnect handler failed:',
                  error
                );
              });
            },

            onmessage: event => {
              try {
                const parsed = parseCalculationUpdateStreamEvent(event.data);

                if (
                  parsed.type === 'calculation_update' &&
                  parsed.calculation
                ) {
                  onCalculationUpdateRef.current?.(parsed.calculation);
                  return;
                }

                if (parsed.type === 'error') {
                  notifyStreamError(
                    parsed.errorMessage ||
                      'Calculation update stream reported an error.'
                  );
                  controller.abort();
                  return;
                }
              } catch (error) {
                notifyStreamError(
                  new Error('Failed to parse calculation update stream event.')
                );
              }
            },

            onclose: () => {
              if (trackGlobalConnection) {
                setGlobalConnected(false);
              }

              if (controller.signal.aborted) {
                return;
              }

              throw new Error('Calculation update stream closed unexpectedly.');
            },

            onerror: error => {
              if (isAbortError(error, controller.signal)) {
                return;
              }

              if (trackGlobalConnection) {
                setGlobalConnected(false);
              }
              notifyStreamError(error);
              return STREAM_RETRY_INTERVAL_MS;
            },
          });
        } catch (error) {
          if (isAbortError(error, controller.signal)) {
            return;
          }

          if (trackGlobalConnection) {
            setGlobalConnected(false);
          }
          notifyStreamError(error);
        }
      })().catch(error => {
        if (isAbortError(error, controller.signal)) {
          return;
        }
        notifyStreamError(error);
      });
    },
    [notifyStreamError, setGlobalConnected]
  );

  const disconnect = useCallback(() => {
    globalAbortRef.current?.abort();
    globalAbortRef.current = null;
    detailAbortRef.current?.abort();
    detailAbortRef.current = null;
    setGlobalConnected(false);
  }, [setGlobalConnected]);

  const reconnect = useCallback(() => {
    startStream(GLOBAL_CALCULATION_UPDATES_PATH, globalAbortRef, true);

    const activeId = activeCalculationIdRef.current;
    if (isUsableCalculationId(activeId)) {
      startStream(
        `/api/quantum/calculations/${encodeURIComponent(activeId)}/updates/stream`,
        detailAbortRef
      );
    }
  }, [startStream]);

  useEffect(() => {
    isMountedRef.current = true;
    startStream(GLOBAL_CALCULATION_UPDATES_PATH, globalAbortRef, true);

    return () => {
      isMountedRef.current = false;
      globalAbortRef.current?.abort();
      globalAbortRef.current = null;
      setIsConnected(false);
    };
  }, [startStream]);

  useEffect(() => {
    if (!isUsableCalculationId(activeCalculationId)) {
      detailAbortRef.current?.abort();
      detailAbortRef.current = null;
      return;
    }

    startStream(
      `/api/quantum/calculations/${encodeURIComponent(activeCalculationId)}/updates/stream`,
      detailAbortRef
    );

    return () => {
      detailAbortRef.current?.abort();
      detailAbortRef.current = null;
    };
  }, [activeCalculationId, startStream]);

  return {
    isConnected,
    reconnect,
    disconnect,
  };
};
