import { useCalculationNotifier } from '../realtime/useCalculationNotifier';
import { useCalculationSync } from '../realtime/useCalculationSync';
import { useCalculationUpdateStream } from '../realtime/useCalculationUpdateStream';

export interface UseCalculationUpdatesOptions {
  activeCalculationId: string | null;
}

export const useCalculationUpdates = ({
  activeCalculationId,
}: UseCalculationUpdatesOptions) => {
  const { notifyCalculationUpdate, handleStreamError } =
    useCalculationNotifier(activeCalculationId);

  const {
    handleCalculationUpdate,
    handleStreamError: handleSyncStreamError,
    handleReconnect,
  } = useCalculationSync({
    activeCalculationId,
    onCalculationUpdate: notifyCalculationUpdate,
    onStreamError: handleStreamError,
  });

  return useCalculationUpdateStream({
    activeCalculationId,
    onCalculationUpdate: handleCalculationUpdate,
    onStreamError: handleSyncStreamError,
    onReconnect: handleReconnect,
  });
};
