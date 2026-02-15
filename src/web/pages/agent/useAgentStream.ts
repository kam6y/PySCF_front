import React, { useCallback, useRef, useState } from 'react';
import type { QueryClient } from '@tanstack/react-query';
import { streamChatWithAgent } from '../../api/agent';
import { useNotificationStore } from '../../store/notificationStore';
import { useAgentStore, AgentStatus, ChatHistory } from '../../store/agentStore';
import { chatHistoryKeys } from '../../hooks/useChatHistoryQueries';

type CreateChatSessionMutation = ReturnType<
  typeof import('../../hooks/useChatHistoryQueries').useCreateChatSession
>;

interface UseAgentStreamOptions {
  currentMessage: string;
  setCurrentMessage: (msg: string) => void;
  isLoading: boolean;
  setIsLoading: (loading: boolean) => void;
  activeSessionId: string | null;
  setActiveSessionId: (id: string | null) => void;
  prevSessionIdRef: React.MutableRefObject<string | null>;
  queryClient: QueryClient;
  createChatSession: CreateChatSessionMutation;
}

interface UseAgentStreamReturn {
  error: string | null;
  handleSendMessage: () => Promise<void>;
  handleCancelMessage: () => void;
}

export function useAgentStream({
  currentMessage,
  setCurrentMessage,
  isLoading,
  setIsLoading,
  activeSessionId,
  setActiveSessionId,
  prevSessionIdRef,
  queryClient,
  createChatSession,
}: UseAgentStreamOptions): UseAgentStreamReturn {
  const history = useAgentStore(state => state.history);
  const addMessages = useAgentStore(state => state.addMessages);
  const updateMessage = useAgentStore(state => state.updateMessage);
  const setHistory = useAgentStore(state => state.setHistory);
  const setAgentStatus = useAgentStore(state => state.setAgentStatus);
  const addNotification = useNotificationStore(state => state.addNotification);

  const [error, setError] = useState<string | null>(null);

  // 重複送信を防止するためのRef（setStateは非同期なので、useRefで同期的に管理）
  const isSendingRef = useRef(false);

  // ストリームのキャンセル関数を保存するRef
  const abortStreamRef = useRef<(() => void) | null>(null);

  // メッセージ送信のキャンセル処理
  const handleCancelMessage = useCallback(() => {
    // ストリームを中断
    if (abortStreamRef.current) {
      abortStreamRef.current();
      abortStreamRef.current = null;
    }

    // 状態をリセット
    setIsLoading(false);
    isSendingRef.current = false;
    setAgentStatus({ agent: null, status: 'idle' });

    // 最後のメッセージ（AIのプレースホルダー）を削除
    const currentHistory = useAgentStore.getState().history;
    if (
      currentHistory.length > 0 &&
      currentHistory[currentHistory.length - 1].role === 'model' &&
      currentHistory[currentHistory.length - 1].isStreaming
    ) {
      // 最後のストリーミング中のメッセージを削除
      setHistory(currentHistory.slice(0, -1));
    }

    // 通知を表示
    addNotification({
      type: 'info',
      title: 'Message Sending Cancelled',
      message: 'The message sending has been cancelled',
      autoClose: true,
      duration: 3000,
    });
  }, [setIsLoading, setAgentStatus, setHistory, addNotification]);

  // メッセージ送信処理
  const handleSendMessage = useCallback(async () => {
    const trimmedMessage = currentMessage.trim();

    // 重複送信を防止（同期的にチェック）
    if (!trimmedMessage || isLoading || isSendingRef.current) {
      return;
    }

    // 送信中フラグを即座に設定（同期的）
    isSendingRef.current = true;

    // チャットセッションの作成（最初のメッセージ送信時のみ）
    let sessionIdToUse = activeSessionId;
    if (!activeSessionId && history.length === 0) {
      try {
        // セッション名を生成（制御文字・改行を削除し、最初の50文字をタイトルとして使用）
        // 1. 制御文字を除去（\x00-\x1F, \x7F-\x9F）
        let cleanedMessage = trimmedMessage.replace(
          /[\x00-\x1F\x7F-\x9F]/g,
          ' '
        );
        // 2. 複数の連続した空白を1つにまとめる
        cleanedMessage = cleanedMessage.replace(/\s+/g, ' ').trim();
        // 3. 先頭と末尾の句読点を除去
        cleanedMessage = cleanedMessage.replace(/^[.,!?;:]+|[.,!?;:]+$/g, '');
        // 4. 最大長チェック（50文字）
        const sessionName =
          cleanedMessage.length > 50
            ? cleanedMessage.substring(0, 47).trim() + '...'
            : cleanedMessage || '新しいチャット'; // 空の場合はデフォルト名

        const result = await createChatSession.mutateAsync(sessionName);
        sessionIdToUse = result.session.id;

        // セッションIDを即座に設定（送信完了を待たずに設定して整合性を保つ）
        setActiveSessionId(sessionIdToUse);

        // 新しいセッションを作成した直後は、prevSessionIdRefも更新して
        // useEffectが空の履歴をロードしようとするのを防ぐ
        prevSessionIdRef.current = sessionIdToUse;

        // キャッシュに新規セッションの詳細を即座に追加（404エラーを防ぐ）
        queryClient.setQueryData(
          chatHistoryKeys.sessionDetail(sessionIdToUse),
          {
            session: result.session,
            messages: [], // 空のメッセージ配列で初期化
          }
        );
      } catch (error) {
        // セッション作成に失敗しても、会話は継続可能
        // ただし、会話履歴はデータベースに保存されない（一時的な会話として扱う）
        console.warn('Failed to create chat session:', error);
        addNotification({
          type: 'info',
          title: 'Session Not Saved',
          message:
            'Unable to save this conversation. You can continue chatting, but the history will not be saved.',
          autoClose: false,
          duration: 0,
        });
        // Continue without session_id (conversation won't be saved, but chat continues)
        sessionIdToUse = null;
      }
    }

    const userMessage: ChatHistory = {
      role: 'user',
      parts: [{ text: trimmedMessage }],
    };

    // ユニークIDを使用してメッセージを識別
    const tempMessageId = crypto.randomUUID();
    const aiPlaceholder: ChatHistory = {
      role: 'model',
      parts: [{ text: '' }],
      isStreaming: true,
      tempId: tempMessageId, // 一意識別子を追加
    };

    // ユーザーメッセージとAIのプレースホルダーを一括で履歴に追加
    // 一括追加により、状態更新の原子性を保証
    addMessages([userMessage, aiPlaceholder]);

    setCurrentMessage('');
    setIsLoading(true);
    setError(null);

    // ストリームを開始し、キャンセル関数を保存
    const abortFunction = streamChatWithAgent(
      trimmedMessage,
      history, // ストリーム開始前の履歴を渡す
      sessionIdToUse, // session_idを渡す
      {
        onMessage: chunk => {
          const currentHistory = useAgentStore.getState().history;
          const lastIndex = currentHistory.length - 1;
          const lastMessage = currentHistory[lastIndex];

          if (
            lastMessage &&
            lastMessage.role === 'model' &&
            lastMessage.tempId === tempMessageId &&
            lastMessage.isStreaming
          ) {
            // メッセージを更新（テキストを追加）
            updateMessage(lastIndex, {
              parts: [{ text: lastMessage.parts[0].text + chunk }],
            });
          }
        },
        onAgentStatus: (status, agent) => {
          // エージェントステータスを更新
          setAgentStatus({
            agent: agent as AgentStatus['agent'],
            status:
              status === 'running' || status === 'responding' ? status : 'idle',
          });
        },
        onClose: () => {
          setIsLoading(false);
          isSendingRef.current = false; // 送信完了フラグをリセット
          abortStreamRef.current = null; // キャンセル関数をクリア
          setAgentStatus({ agent: null, status: 'idle' }); // ステータスをリセット

          const currentHistory = useAgentStore.getState().history;
          const lastIndex = currentHistory.length - 1;
          const lastMessage = currentHistory[lastIndex];

          if (
            lastMessage &&
            lastMessage.role === 'model' &&
            lastMessage.tempId === tempMessageId
          ) {
            // ストリーミング終了時にtempIdを削除し、isStreamingをfalseに
            updateMessage(lastIndex, {
              isStreaming: false,
            });
          }

          // AI応答がデータベースに保存された後、セッション詳細とセッション一覧のキャッシュを無効化
          // これにより、次回の履歴確認時に最新データがフェッチされ、サイドバーのメッセージ数も即座に更新される
          if (sessionIdToUse) {
            queryClient.invalidateQueries({
              queryKey: chatHistoryKeys.sessionDetail(sessionIdToUse),
            });
            queryClient.invalidateQueries({
              queryKey: chatHistoryKeys.sessions(),
            });
          }
        },
        onError: err => {
          console.error('AI Agent Error:', err);
          setIsLoading(false);
          isSendingRef.current = false; // エラー時も送信フラグをリセット
          abortStreamRef.current = null; // エラー時もキャンセル関数をクリア
          setAgentStatus({ agent: null, status: 'idle' }); // エラー時もステータスをリセット
          setError(err.message);
          const currentHistory = useAgentStore.getState().history;
          const lastIndex = currentHistory.length - 1;
          const lastMessage = currentHistory[lastIndex];

          if (
            lastMessage &&
            lastMessage.role === 'model' &&
            lastMessage.tempId === tempMessageId
          ) {
            // エラー時もメッセージを更新
            updateMessage(lastIndex, {
              parts: [{ text: `Error: ${err.message}` }],
              isStreaming: false,
            });
          }

          // 通知でエラーを表示
          addNotification({
            type: 'error',
            title: 'AI Agent Error',
            message: 'Failed to send message to AI agent. Please try again.',
            autoClose: false,
            duration: 0,
          });
        },
      }
    );

    // キャンセル関数を保存
    abortStreamRef.current = abortFunction;
  }, [
    currentMessage,
    history,
    isLoading,
    addMessages,
    updateMessage,
    addNotification,
    setAgentStatus,
    createChatSession,
    activeSessionId,
    queryClient,
    setActiveSessionId,
    prevSessionIdRef,
    setCurrentMessage,
    setIsLoading,
  ]);

  return {
    error,
    handleSendMessage,
    handleCancelMessage,
  };
}

