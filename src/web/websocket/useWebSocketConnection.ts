import { useEffect } from 'react';
import { CalculationInstance } from '../types/api-types';
import { useCalculationSync } from './useCalculationSync';
import { useSocketTransport } from './useSocketTransport';

export interface UseWebSocketConnectionOptions {
  activeCalculationId: string | null;
  onCalculationUpdate?: (
    updatedCalculation: CalculationInstance,
    previousStatus: string | undefined
  ) => void;
  onWebSocketError?: (error: string) => void;
}

export const useWebSocketConnection = ({
  activeCalculationId,
  onCalculationUpdate,
  onWebSocketError,
}: UseWebSocketConnectionOptions) => {
  const { transportHandlers, manageActiveCalculationRoom } = useCalculationSync({
    activeCalculationId,
    onCalculationUpdate,
    onWebSocketError,
  });

  const { isConnected, reconnect, disconnect, emit } = useSocketTransport({
    handlers: transportHandlers,
  });

  useEffect(() => {
    manageActiveCalculationRoom(activeCalculationId, emit);
  }, [activeCalculationId, manageActiveCalculationRoom, emit]);

  return {
    isConnected,
    reconnect,
    disconnect,
  };
};
