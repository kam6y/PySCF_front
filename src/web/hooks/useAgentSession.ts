import React, { useCallback, useEffect, useRef, useState } from 'react';
import { QueryClient, useQueryClient } from '@tanstack/react-query';
import { useNotificationStore } from '../store/notificationStore';
import { useAgentStore, ChatHistory } from '../store/agentStore';
import { useChatHistoryStore } from '../store/chatHistoryStore';
import {
  useCreateChatSession,
  useGetChatSessionDetail,
  useUpdateChatSession,
} from './useChatHistoryQueries';

type SessionDetailData = ReturnType<typeof useGetChatSessionDetail>['data'];

interface UseAgentSessionOptions {
  isLoading: boolean;
}

interface UseAgentSessionReturn {
  activeSessionId: string | null;
  setActiveSessionId: (id: string | null) => void;
  sessionDetailData: SessionDetailData;
  queryClient: QueryClient;
  createChatSession: ReturnType<typeof useCreateChatSession>;
  prevSessionIdRef: React.MutableRefObject<string | null>;
  // 新規チャット
  isNewChatConfirmationOpen: boolean;
  handleNewChatRequest: () => void;
  handleNewChatConfirm: () => Promise<void>;
  handleNewChatCancel: () => void;
  // タイトル編集
  isEditingTitle: boolean;
  editedTitle: string;
  setEditedTitle: (title: string) => void;
  handleStartEditTitle: () => void;
  handleSaveTitle: () => Promise<void>;
  handleCancelEditTitle: () => void;
  handleTitleKeyDown: (e: React.KeyboardEvent<HTMLInputElement>) => void;
}

export function useAgentSession({
  isLoading,
}: UseAgentSessionOptions): UseAgentSessionReturn {
  const setHistory = useAgentStore(state => state.setHistory);
  const clearHistory = useAgentStore(state => state.clearHistory);
  const historyLength = useAgentStore(state => state.history.length);

  // チャット履歴ストア（セッションID管理を一元化）
  const activeSessionId = useChatHistoryStore(state => state.activeSessionId);
  const setActiveSessionId = useChatHistoryStore(
    state => state.setActiveSessionId
  );
  const clearActiveSession = useChatHistoryStore(
    state => state.clearActiveSession
  );

  // TanStack Query
  const queryClient = useQueryClient();

  // チャット履歴のクエリ
  const createChatSession = useCreateChatSession();
  const { data: sessionDetailData } = useGetChatSessionDetail(activeSessionId);
  const updateChatSession = useUpdateChatSession();

  // タイトル編集用の状態
  const [isEditingTitle, setIsEditingTitle] = useState(false);
  const [editedTitle, setEditedTitle] = useState('');
  const addNotification = useNotificationStore(state => state.addNotification);

  // 前回のセッションIDを追跡（競合状態を回避）
  const prevSessionIdRef = useRef<string | null>(null);

  // New chat confirmation modal state
  const [isNewChatConfirmationOpen, setIsNewChatConfirmationOpen] =
    useState(false);

  // Load session detail when activeSessionId changes
  // This effect handles three scenarios:
  // 1. Switching to a new session → Load history from database
  // 2. Clearing session (null) → Clear history display
  // 3. During message streaming → Skip to preserve real-time updates
  useEffect(() => {
    // Guard: Don't modify history during message streaming to preserve real-time updates
    if (isLoading) {
      return;
    }

    // Guard: Only process if session ID actually changed
    const hasSessionChanged = activeSessionId !== prevSessionIdRef.current;
    if (!hasSessionChanged) {
      return;
    }

    // Case 1: Session was cleared (null) - immediately clear history
    if (!activeSessionId) {
      prevSessionIdRef.current = null;
      clearHistory();
      return;
    }

    // Case 2: Switched to a new session - load history from database
    // Wait for sessionDetailData to be available before loading
    if (sessionDetailData) {
      prevSessionIdRef.current = activeSessionId;
      const loadedHistory: ChatHistory[] = sessionDetailData.messages.map(
        msg => ({
          role: msg.role as 'user' | 'model',
          parts: [{ text: msg.content }],
        })
      );
      setHistory(loadedHistory);
    }
    // Note: If sessionDetailData is not ready yet, this effect will re-run
    // when it becomes available (dependency array includes sessionDetailData)
  }, [activeSessionId, sessionDetailData, setHistory, clearHistory, isLoading]);

  // Handle new chat confirm (defined first to avoid reference error)
  const handleNewChatConfirm = useCallback(async () => {
    // 現在の会話をクリア
    clearHistory();
    clearActiveSession();
    // prevSessionIdRefもリセット
    prevSessionIdRef.current = null;
    setIsNewChatConfirmationOpen(false);

    addNotification({
      type: 'success',
      title: 'New Chat Started',
      message: 'Started a new conversation.',
      autoClose: true,
      duration: 3000,
    });
  }, [
    clearHistory,
    clearActiveSession,
    addNotification,
    setIsNewChatConfirmationOpen,
  ]);

  // Handle new chat request
  const handleNewChatRequest = useCallback(() => {
    if (historyLength > 0) {
      // 会話がある場合は確認モーダルを表示
      setIsNewChatConfirmationOpen(true);
    } else {
      // 会話がない場合は直接新しいチャットを開始
      handleNewChatConfirm();
    }
  }, [historyLength, handleNewChatConfirm, setIsNewChatConfirmationOpen]);

  const handleNewChatCancel = useCallback(() => {
    setIsNewChatConfirmationOpen(false);
  }, [setIsNewChatConfirmationOpen]);

  // タイトル編集の開始
  const handleStartEditTitle = useCallback(() => {
    if (sessionDetailData?.session?.name) {
      setEditedTitle(sessionDetailData.session.name);
      setIsEditingTitle(true);
    }
  }, [sessionDetailData]);

  // タイトル編集のキャンセル
  const handleCancelEditTitle = useCallback(() => {
    setIsEditingTitle(false);
    setEditedTitle('');
  }, []);

  // タイトルの保存
  const handleSaveTitle = useCallback(async () => {
    const trimmedTitle = editedTitle.trim();

    // バリデーション
    if (!trimmedTitle) {
      addNotification({
        type: 'error',
        title: 'Invalid Title',
        message: 'Title cannot be empty.',
        autoClose: true,
        duration: 3000,
      });
      return;
    }

    if (trimmedTitle.length > 100) {
      addNotification({
        type: 'error',
        title: 'Invalid Title',
        message: 'Title must be 100 characters or less.',
        autoClose: true,
        duration: 3000,
      });
      return;
    }

    if (!activeSessionId) {
      return;
    }

    try {
      await updateChatSession.mutateAsync({
        sessionId: activeSessionId,
        name: trimmedTitle,
      });

      setIsEditingTitle(false);
      setEditedTitle('');

      addNotification({
        type: 'success',
        title: 'Title Updated',
        message: 'Conversation title has been updated.',
        autoClose: true,
        duration: 3000,
      });
    } catch (error) {
      console.error('Failed to update title:', error);
      addNotification({
        type: 'error',
        title: 'Update Failed',
        message: 'Failed to update conversation title. Please try again.',
        autoClose: true,
        duration: 5000,
      });
    }
  }, [editedTitle, activeSessionId, updateChatSession, addNotification]);

  // Enter キーで保存、Escape キーでキャンセル
  const handleTitleKeyDown = useCallback(
    (e: React.KeyboardEvent<HTMLInputElement>) => {
      if (e.key === 'Enter') {
        e.preventDefault();
        handleSaveTitle();
      } else if (e.key === 'Escape') {
        e.preventDefault();
        handleCancelEditTitle();
      }
    },
    [handleSaveTitle, handleCancelEditTitle]
  );

  return {
    activeSessionId,
    setActiveSessionId,
    sessionDetailData,
    queryClient,
    createChatSession,
    prevSessionIdRef,
    isNewChatConfirmationOpen,
    handleNewChatRequest,
    handleNewChatConfirm,
    handleNewChatCancel,
    isEditingTitle,
    editedTitle,
    setEditedTitle,
    handleStartEditTitle,
    handleSaveTitle,
    handleCancelEditTitle,
    handleTitleKeyDown,
  };
}
