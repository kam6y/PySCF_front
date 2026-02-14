import { useCallback, useState } from 'react';
import { useDeleteChatSession } from '../../hooks/useChatHistoryQueries';
import { useAgentStore } from '../../store/agentStore';
import { useChatHistoryStore } from '../../store/chatHistoryStore';
import { CalculationSummary } from '../../types/api-types';

interface UseDeleteModalsDeps {
  onCalculationDelete: (calculationId: string) => Promise<void>;
}

export interface DeleteModalsState {
  isDeleteModalOpen: boolean;
  calculationToDelete: {
    id: string;
    name: string;
  } | null;
  handleRequestDelete: (calculationId: string, calculationName: string) => void;
  handleConfirmDelete: () => Promise<void>;
  handleCancelDelete: () => void;
  isBulkDeleteModalOpen: boolean;
  calculationsToDelete: CalculationSummary[];
  handleBulkDeleteError: (errorCalculations: CalculationSummary[]) => void;
  handleConfirmBulkDelete: () => Promise<void>;
  handleCancelBulkDelete: () => void;
  isChatDeleteModalOpen: boolean;
  chatToDelete: {
    id: string;
    name: string;
  } | null;
  handleRequestChatDelete: (sessionId: string, sessionName: string) => void;
  handleConfirmChatDelete: () => Promise<void>;
  handleCancelChatDelete: () => void;
  deleteChatSessionPending: boolean;
}

export function useDeleteModals({
  onCalculationDelete,
}: UseDeleteModalsDeps): DeleteModalsState {
  const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
  const [calculationToDelete, setCalculationToDelete] = useState<{
    id: string;
    name: string;
  } | null>(null);

  const [isBulkDeleteModalOpen, setIsBulkDeleteModalOpen] = useState(false);
  const [calculationsToDelete, setCalculationsToDelete] = useState<
    CalculationSummary[]
  >([]);

  const [isChatDeleteModalOpen, setIsChatDeleteModalOpen] = useState(false);
  const [chatToDelete, setChatToDelete] = useState<{
    id: string;
    name: string;
  } | null>(null);

  const deleteChatSession = useDeleteChatSession();
  const activeSessionId = useChatHistoryStore(state => state.activeSessionId);
  const clearActiveSession = useChatHistoryStore(
    state => state.clearActiveSession
  );
  const clearHistory = useAgentStore(state => state.clearHistory);

  const handleRequestChatDelete = useCallback(
    (sessionId: string, sessionName: string) => {
      setChatToDelete({ id: sessionId, name: sessionName });
      setIsChatDeleteModalOpen(true);
    },
    []
  );

  const handleConfirmChatDelete = useCallback(async () => {
    if (!chatToDelete) {
      return;
    }

    if (chatToDelete.id === activeSessionId) {
      clearActiveSession();
      clearHistory();
    }

    await deleteChatSession.mutateAsync(chatToDelete.id);
    setIsChatDeleteModalOpen(false);
    setChatToDelete(null);
  }, [
    chatToDelete,
    activeSessionId,
    clearActiveSession,
    clearHistory,
    deleteChatSession,
  ]);

  const handleCancelChatDelete = useCallback(() => {
    setIsChatDeleteModalOpen(false);
    setChatToDelete(null);
  }, []);

  const handleBulkDeleteError = useCallback(
    (errorCalculations: CalculationSummary[]) => {
      if (errorCalculations.length === 0) {
        return;
      }

      setCalculationsToDelete(errorCalculations);
      setIsBulkDeleteModalOpen(true);
    },
    []
  );

  const handleConfirmBulkDelete = useCallback(async () => {
    try {
      await Promise.all(
        calculationsToDelete.map(calculation =>
          onCalculationDelete(calculation.id)
        )
      );
      setIsBulkDeleteModalOpen(false);
      setCalculationsToDelete([]);
    } catch (error) {
      console.error('一括削除中にエラーが発生しました:', error);
      alert('Failed to delete some calculations.');
    }
  }, [calculationsToDelete, onCalculationDelete]);

  const handleCancelBulkDelete = useCallback(() => {
    setIsBulkDeleteModalOpen(false);
    setCalculationsToDelete([]);
  }, []);

  const handleRequestDelete = useCallback(
    (calculationId: string, calculationName: string) => {
      setCalculationToDelete({ id: calculationId, name: calculationName });
      setIsDeleteModalOpen(true);
    },
    []
  );

  const handleConfirmDelete = useCallback(async () => {
    if (!calculationToDelete) {
      return;
    }

    await onCalculationDelete(calculationToDelete.id);
    setIsDeleteModalOpen(false);
    setCalculationToDelete(null);
  }, [calculationToDelete, onCalculationDelete]);

  const handleCancelDelete = useCallback(() => {
    setIsDeleteModalOpen(false);
    setCalculationToDelete(null);
  }, []);

  return {
    isDeleteModalOpen,
    calculationToDelete,
    handleRequestDelete,
    handleConfirmDelete,
    handleCancelDelete,
    isBulkDeleteModalOpen,
    calculationsToDelete,
    handleBulkDeleteError,
    handleConfirmBulkDelete,
    handleCancelBulkDelete,
    isChatDeleteModalOpen,
    chatToDelete,
    handleRequestChatDelete,
    handleConfirmChatDelete,
    handleCancelChatDelete,
    deleteChatSessionPending: deleteChatSession.isPending,
  };
}
