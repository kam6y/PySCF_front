import { useWebSocketConnection } from '../websocket/useWebSocketConnection';
import { useCalculationNotifier } from '../websocket/useCalculationNotifier';

export interface UseUnifiedWebSocketOptions {
  activeCalculationId: string | null;
}

export const useUnifiedWebSocket = ({
  activeCalculationId,
}: UseUnifiedWebSocketOptions) => {
  const { notifyCalculationUpdate, handleWebSocketError } =
    useCalculationNotifier(activeCalculationId);

  return useWebSocketConnection({
    activeCalculationId,
    onCalculationUpdate: notifyCalculationUpdate,
    onWebSocketError: handleWebSocketError,
  });
};
