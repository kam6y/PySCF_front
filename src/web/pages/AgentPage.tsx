import React, { useState, useCallback, useRef, useEffect } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeHighlight from 'rehype-highlight';
import { useQueryClient } from '@tanstack/react-query';
import { streamChatWithAgent, executeConfirmedAgentAction } from '../apiClient';
import { useNotificationStore } from '../store/notificationStore';
import { useAgentStore, ChatHistory, AgentStatus } from '../store/agentStore';
import { useChatHistoryStore } from '../store/chatHistoryStore';
import {
  useCreateChatSession,
  useGetChatSessionDetail,
  useUpdateChatSession,
  chatHistoryKeys,
} from '../hooks/useChatHistoryQueries';
import { ConfirmationModal } from '../components/ConfirmationModal';
import { InlineOrbitalViewer } from '../components/InlineOrbitalViewer';
import styles from './AgentPage.module.css';

// 確認リクエストの型定義
type ConfirmationRequest = {
  requires_confirmation: boolean;
  action: string;
  calculation_id: string;
  calculation_name: string;
  message: string;
};

// パラメータパース関数（orbital-viewerコードブロック用）
const parseOrbitalViewerParams = (code: string): Record<string, any> => {
  const params: Record<string, any> = {};
  const lines = code.trim().split('\n');

  lines.forEach(line => {
    const colonIndex = line.indexOf(':');
    if (colonIndex === -1) return;

    const key = line.substring(0, colonIndex).trim();
    const value = line.substring(colonIndex + 1).trim();

    if (!key || !value) return;

    // 数値に変換できる場合は数値として扱う
    const numValue = Number(value);
    params[key] = isNaN(numValue) ? value : numValue;
  });

  return params;
};

export const AgentPage = () => {
  // Zustandストアから会話履歴とエージェントステータスを取得
  const history = useAgentStore(state => state.history);
  const currentAgentStatus = useAgentStore(state => state.currentAgentStatus);
  const addMessage = useAgentStore(state => state.addMessage);
  const addMessages = useAgentStore(state => state.addMessages);
  const updateMessage = useAgentStore(state => state.updateMessage);
  const setHistory = useAgentStore(state => state.setHistory);
  const clearHistory = useAgentStore(state => state.clearHistory);
  const setAgentStatus = useAgentStore(state => state.setAgentStatus);

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

  const [currentMessage, setCurrentMessage] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // タイトル編集用の状態
  const [isEditingTitle, setIsEditingTitle] = useState(false);
  const [editedTitle, setEditedTitle] = useState('');
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const addNotification = useNotificationStore(state => state.addNotification);

  // 重複送信を防止するためのRef（setStateは非同期なので、useRefで同期的に管理）
  const isSendingRef = useRef(false);

  // 前回のセッションIDを追跡（競合状態を回避）
  const prevSessionIdRef = useRef<string | null>(null);

  // Confirmation modal state
  const [confirmationRequest, setConfirmationRequest] =
    useState<ConfirmationRequest | null>(null);
  const [isConfirmationModalOpen, setIsConfirmationModalOpen] = useState(false);
  const [isExecutingAction, setIsExecutingAction] = useState(false);

  // New chat confirmation modal state
  const [isNewChatConfirmationOpen, setIsNewChatConfirmationOpen] =
    useState(false);

  // Parse confirmation request from AI message
  const parseConfirmationRequest = useCallback(
    (text: string): ConfirmationRequest | null => {
      try {
        // Try to extract JSON from markdown code blocks first (```json ... ```)
        const jsonBlockRegex = /```json\s*(\{[\s\S]*?\})\s*```/;
        const jsonBlockMatch = text.match(jsonBlockRegex);
        if (jsonBlockMatch?.[1]) {
          const parsed = JSON.parse(jsonBlockMatch[1]);
          if (
            parsed.requires_confirmation === true &&
            parsed.action &&
            parsed.calculation_id
          ) {
            return parsed as ConfirmationRequest;
          }
        }

        // Try to extract plain JSON object (more robust pattern)
        // This pattern finds JSON objects that contain "requires_confirmation": true
        // and handles multi-line structures with nested content
        const jsonPatternRegex =
          /\{[^{}]*"requires_confirmation"\s*:\s*true[^{}]*(?:\{[^{}]*\}[^{}]*)*\}/g;
        const matches = text.match(jsonPatternRegex);

        if (matches) {
          // Try to parse each match
          for (const match of matches) {
            try {
              const parsed = JSON.parse(match);
              if (
                parsed.requires_confirmation === true &&
                parsed.action &&
                parsed.calculation_id
              ) {
                return parsed as ConfirmationRequest;
              }
            } catch (parseError) {
              continue;
            }
          }
        }

        // Try to find JSON object with a more lenient approach
        // Find all { ... } blocks and try to parse them
        const allJsonRegex = /\{(?:[^{}]|\{[^{}]*\})*\}/g;
        const allMatches = text.match(allJsonRegex);

        if (allMatches) {
          for (const match of allMatches) {
            try {
              const parsed = JSON.parse(match);
              if (
                parsed.requires_confirmation === true &&
                parsed.action &&
                parsed.calculation_id
              ) {
                return parsed as ConfirmationRequest;
              }
            } catch (parseError) {
              // Silently skip invalid JSON
              continue;
            }
          }
        }

        return null;
      } catch (e) {
        console.error(
          '[ConfirmationParser] Error parsing confirmation request:',
          e
        );
        return null;
      }
    },
    []
  );

  // Handle confirmation modal confirm action
  const handleConfirmAction = useCallback(async () => {
    if (!confirmationRequest) return;

    setIsExecutingAction(true);
    try {
      const result = await executeConfirmedAgentAction(
        confirmationRequest.action as 'delete_calculation',
        confirmationRequest.calculation_id
      );

      // Close modal
      setIsConfirmationModalOpen(false);
      setConfirmationRequest(null);

      // Show success notification
      addNotification({
        type: 'success',
        title: 'Action Completed',
        message: result.message || 'The action was completed successfully.',
        autoClose: true,
        duration: 5000,
      });

      // Add AI response to history
      const successMessage: ChatHistory = {
        role: 'model',
        parts: [
          { text: `✅ ${result.message || 'Action completed successfully.'}` },
        ],
      };
      addMessage(successMessage);
    } catch (err: any) {
      // Show error notification
      addNotification({
        type: 'error',
        title: 'Action Failed',
        message: err.message || 'Failed to execute the action.',
        autoClose: false,
        duration: 0,
      });

      // Add error message to history
      const errorMessage: ChatHistory = {
        role: 'model',
        parts: [
          {
            text: `❌ Failed to execute action: ${err.message || 'Unknown error'}`,
          },
        ],
      };
      addMessage(errorMessage);
    } finally {
      setIsExecutingAction(false);
    }
  }, [confirmationRequest, addNotification]);

  // Handle confirmation modal cancel action
  const handleCancelAction = useCallback(() => {
    setIsConfirmationModalOpen(false);
    setConfirmationRequest(null);

    // Add cancellation message to history
    const cancelMessage: ChatHistory = {
      role: 'model',
      parts: [{ text: '❌ Action cancelled by user.' }],
    };
    addMessage(cancelMessage);
  }, [addMessage]);

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
  }, [
    activeSessionId,
    sessionDetailData,
    setHistory,
    clearHistory,
    isLoading,
    queryClient,
  ]);

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
    if (history.length > 0) {
      // 会話がある場合は確認モーダルを表示
      setIsNewChatConfirmationOpen(true);
    } else {
      // 会話がない場合は直接新しいチャットを開始
      handleNewChatConfirm();
    }
  }, [history.length, handleNewChatConfirm, setIsNewChatConfirmationOpen]);

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
    const tempMessageId = Date.now();
    const aiPlaceholder: ChatHistory = {
      role: 'model',
      parts: [{ text: '' }],
      isStreaming: true,
      tempId: tempMessageId, // 一意識別子を追加
    };

    // ユーザーメッセージとAIのプレースホルダーを一括で履歴に追加
    // 一括追加により、状態更新の原子性を保証
    addMessages([userMessage, aiPlaceholder]);

    const aiMessageIndex = history.length + 1; // AIメッセージのインデックス
    setCurrentMessage('');
    setIsLoading(true);
    setError(null);

    streamChatWithAgent(
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

            // Check for confirmation request in the final message
            const messageText = lastMessage.parts[0]?.text || '';
            const confirmation = parseConfirmationRequest(messageText);
            if (confirmation) {
              // Show confirmation modal
              setConfirmationRequest(confirmation);
              setIsConfirmationModalOpen(true);
            }
          }

          // AI応答がデータベースに保存された後、セッション詳細のキャッシュを無効化
          // これにより、次回の履歴確認時に最新データがフェッチされる
          if (sessionIdToUse) {
            queryClient.invalidateQueries({
              queryKey: chatHistoryKeys.sessionDetail(sessionIdToUse),
            });
          }
        },
        onError: err => {
          console.error('AI Agent Error:', err);
          setIsLoading(false);
          isSendingRef.current = false; // エラー時も送信フラグをリセット
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
    setActiveSessionId,
    parseConfirmationRequest,
  ]);

  // Enterで送信（Shift+Enterで改行）
  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        handleSendMessage();
      }
    },
    [handleSendMessage]
  );

  // テキストエリアの高さを自動調整
  const adjustTextareaHeight = useCallback(() => {
    const textarea = textareaRef.current;
    if (textarea) {
      textarea.style.height = 'auto';
      textarea.style.height = `${Math.min(textarea.scrollHeight, 120)}px`;
    }
  }, []);

  useEffect(() => {
    adjustTextareaHeight();
  }, [currentMessage, adjustTextareaHeight]);

  // 新しいメッセージが追加されたときにスクロール
  useEffect(() => {
    if (chatWindowRef.current) {
      chatWindowRef.current.scrollTop = chatWindowRef.current.scrollHeight;
    }
  }, [history]);

  // 空の状態表示
  const renderEmptyState = () => (
    <div className={styles.emptyState}>
      <svg
        className={styles.emptyStateIcon}
        fill="none"
        stroke="currentColor"
        viewBox="0 0 24 24"
        xmlns="http://www.w3.org/2000/svg"
      >
        <path
          strokeLinecap="round"
          strokeLinejoin="round"
          strokeWidth={1.5}
          d="M8 12h.01M12 12h.01M16 12h.01M21 12c0 4.418-4.03 8-9 8a9.863 9.863 0 01-4.255-.949L3 20l1.395-3.72C3.512 15.042 3 13.574 3 12c0-4.418 4.03-8 9-8s9 3.582 9 8z"
        />
      </svg>
      <div className={styles.emptyStateTitle}>Welcome to AI Agent</div>
      <div className={styles.emptyStateText}>
        Start a conversation with the AI agent to get help with molecular
        design, quantum chemistry calculations, and analysis. Ask questions
        about your molecules or request assistance with quantum computational
        chemistry tasks.
      </div>
    </div>
  );

  // タイピングインジケーター
  const renderTypingIndicator = () => (
    <div className={styles.typingIndicator}>
      <div className={styles.typingDots}>
        <div className={styles.typingDot}></div>
        <div className={styles.typingDot}></div>
        <div className={styles.typingDot}></div>
      </div>
    </div>
  );

  // エージェント表示名を取得
  const getAgentDisplayName = (agent: string) => {
    const names: Record<string, string> = {
      supervisor: 'Supervisor',
      quantum_calculation_worker: 'Quantum Calculation Worker',
      research_expert: 'Research Expert',
    };
    return names[agent] || agent;
  };

  return (
    <div className={styles.agentPageContainer}>
      {/* Title and New Chat button container - only show when there's an active session */}
      {activeSessionId && sessionDetailData?.session && (
        <div className={styles.headerContainer}>
          {/* Conversation Title */}
          <div className={styles.titleContainer}>
            {isEditingTitle ? (
              <div className={styles.titleEditMode}>
                <input
                  type="text"
                  className={styles.titleInput}
                  value={editedTitle}
                  onChange={e => setEditedTitle(e.target.value)}
                  onKeyDown={handleTitleKeyDown}
                  autoFocus
                  maxLength={100}
                />
                <div className={styles.titleActionButtons}>
                  <button
                    className={styles.titleSaveButton}
                    onClick={handleSaveTitle}
                    title="Save"
                  >
                    <svg
                      width="16"
                      height="16"
                      fill="none"
                      stroke="currentColor"
                      viewBox="0 0 24 24"
                    >
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        strokeWidth={2}
                        d="M5 13l4 4L19 7"
                      />
                    </svg>
                  </button>
                  <button
                    className={styles.titleCancelButton}
                    onClick={handleCancelEditTitle}
                    title="Cancel"
                  >
                    <svg
                      width="16"
                      height="16"
                      fill="none"
                      stroke="currentColor"
                      viewBox="0 0 24 24"
                    >
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        strokeWidth={2}
                        d="M6 18L18 6M6 6l12 12"
                      />
                    </svg>
                  </button>
                </div>
              </div>
            ) : (
              <div className={styles.titleDisplay}>
                <h2 className={styles.titleText}>
                  {sessionDetailData.session.name}
                </h2>
                <button
                  className={styles.editButton}
                  onClick={handleStartEditTitle}
                  title="Edit title"
                >
                  <svg
                    width="16"
                    height="16"
                    fill="none"
                    stroke="currentColor"
                    viewBox="0 0 24 24"
                  >
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth={2}
                      d="M15.232 5.232l3.536 3.536m-2.036-5.036a2.5 2.5 0 113.536 3.536L6.5 21.036H3v-3.572L16.732 3.732z"
                    />
                  </svg>
                </button>
              </div>
            )}
          </div>

          {/* New Chat button */}
          <div className={styles.clearButtonContainer}>
            <button
              className={styles.clearButton}
              onClick={handleNewChatRequest}
              title="Start a new chat"
            >
              <svg
                width="18"
                height="18"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M12 4v16m8-8H4"
                />
              </svg>
              New Chat
            </button>
          </div>
        </div>
      )}

      <div className={styles.chatWindow} ref={chatWindowRef}>
        {history.length === 0 ? (
          renderEmptyState()
        ) : (
          <>
            {history.map((entry, index) => (
              <div
                key={index}
                className={`${styles.chatMessage} ${styles[entry.role]}`}
              >
                <div className={styles.messageContent}>
                  {entry.role === 'model' ? (
                    <>
                      <ReactMarkdown
                        remarkPlugins={[remarkGfm]}
                        rehypePlugins={[rehypeHighlight]}
                        disallowedElements={[
                          'script',
                          'iframe',
                          'object',
                          'embed',
                        ]}
                        unwrapDisallowed={true}
                        className={styles.markdown}
                        components={{
                          pre: ({ children, ...props }: any) => {
                            // Debug logging for pre blocks
                            console.log('[Pre Component Debug]', {
                              childrenType: typeof children,
                              children: children,
                            });

                            // Check if this pre contains an orbital-viewer code block
                            if (
                              React.isValidElement(children) &&
                              children.props
                            ) {
                              const codeProps = children.props as any;
                              const className = codeProps.className;

                              console.log('[Pre > Code Debug]', {
                                className,
                                hasClassName: !!className,
                                isOrbitalViewer: className?.includes(
                                  'language-orbital-viewer'
                                ),
                              });

                              if (
                                className?.includes('language-orbital-viewer')
                              ) {
                                try {
                                  const codeContent = String(
                                    codeProps.children
                                  ).replace(/\n$/, '');
                                  const params =
                                    parseOrbitalViewerParams(codeContent);

                                  console.log(
                                    '[Orbital Viewer] Parsed params:',
                                    params
                                  );

                                  // Validate required parameters
                                  if (
                                    !params.calculation_id ||
                                    params.orbital_index === undefined
                                  ) {
                                    return (
                                      <div
                                        style={{
                                          padding: '12px',
                                          backgroundColor: '#fef2f2',
                                          border: '1px solid #fecaca',
                                          borderRadius: '8px',
                                          color: '#dc2626',
                                          margin: '12px 0',
                                        }}
                                      >
                                        ❌ Invalid orbital-viewer block: missing
                                        required parameters (calculation_id,
                                        orbital_index)
                                      </div>
                                    );
                                  }

                                  // Render InlineOrbitalViewer component
                                  return (
                                    <InlineOrbitalViewer
                                      calculation_id={params.calculation_id}
                                      orbital_index={params.orbital_index}
                                      grid_size={params.grid_size}
                                      isovalue_pos={params.isovalue_pos}
                                      isovalue_neg={params.isovalue_neg}
                                      onError={error => {
                                        console.error(
                                          'Orbital viewer error:',
                                          error
                                        );
                                        addNotification({
                                          type: 'error',
                                          title: 'Orbital Viewer Error',
                                          message: error,
                                          autoClose: false,
                                          duration: 0,
                                        });
                                      }}
                                    />
                                  );
                                } catch (error) {
                                  console.error(
                                    'Failed to render orbital viewer:',
                                    error
                                  );
                                  return (
                                    <div
                                      style={{
                                        padding: '12px',
                                        backgroundColor: '#fef2f2',
                                        border: '1px solid #fecaca',
                                        borderRadius: '8px',
                                        color: '#dc2626',
                                        margin: '12px 0',
                                      }}
                                    >
                                      ❌ Failed to render orbital viewer:{' '}
                                      {error instanceof Error
                                        ? error.message
                                        : String(error)}
                                    </div>
                                  );
                                }
                              }
                            }

                            // Default pre rendering
                            return <pre {...props}>{children}</pre>;
                          },
                          code: ({ className, children, ...props }: any) => {
                            // Debug logging
                            console.log('[Code Component Debug]', {
                              className,
                              inline: props.inline,
                              node: props.node,
                              childrenType: typeof children,
                              childrenValue: String(children).substring(0, 100),
                            });

                            // Check if this is an orbital-viewer code block
                            const inline = props.inline;
                            if (
                              !inline &&
                              className?.includes('language-orbital-viewer')
                            ) {
                              try {
                                const codeContent = String(children).replace(
                                  /\n$/,
                                  ''
                                );
                                const params =
                                  parseOrbitalViewerParams(codeContent);

                                // Validate required parameters
                                if (
                                  !params.calculation_id ||
                                  params.orbital_index === undefined
                                ) {
                                  return (
                                    <div
                                      style={{
                                        padding: '12px',
                                        backgroundColor: '#fef2f2',
                                        border: '1px solid #fecaca',
                                        borderRadius: '8px',
                                        color: '#dc2626',
                                        margin: '12px 0',
                                      }}
                                    >
                                      ❌ Invalid orbital-viewer block: missing
                                      required parameters (calculation_id,
                                      orbital_index)
                                    </div>
                                  );
                                }

                                // Render InlineOrbitalViewer component
                                return (
                                  <InlineOrbitalViewer
                                    calculation_id={params.calculation_id}
                                    orbital_index={params.orbital_index}
                                    grid_size={params.grid_size}
                                    isovalue_pos={params.isovalue_pos}
                                    isovalue_neg={params.isovalue_neg}
                                    onError={error => {
                                      console.error(
                                        'Orbital viewer error:',
                                        error
                                      );
                                      addNotification({
                                        type: 'error',
                                        title: 'Orbital Viewer Error',
                                        message: error,
                                        autoClose: false,
                                        duration: 0,
                                      });
                                    }}
                                  />
                                );
                              } catch (error) {
                                console.error(
                                  'Failed to render orbital viewer:',
                                  error
                                );
                                return (
                                  <div
                                    style={{
                                      padding: '12px',
                                      backgroundColor: '#fef2f2',
                                      border: '1px solid #fecaca',
                                      borderRadius: '8px',
                                      color: '#dc2626',
                                      margin: '12px 0',
                                    }}
                                  >
                                    ❌ Failed to render orbital viewer:{' '}
                                    {error instanceof Error
                                      ? error.message
                                      : String(error)}
                                  </div>
                                );
                              }
                            }

                            // Default code block rendering
                            return (
                              <code className={className} {...props}>
                                {children}
                              </code>
                            );
                          },
                          a: ({ href, children, ...props }: any) => {
                            // Handle external links - open in default browser
                            const handleClick = async (e: React.MouseEvent) => {
                              if (
                                href &&
                                (href.startsWith('http://') ||
                                  href.startsWith('https://'))
                              ) {
                                e.preventDefault();
                                try {
                                  const result =
                                    await window.electronAPI.openExternalUrl(
                                      href
                                    );
                                  if (!result.success) {
                                    console.error(
                                      'Failed to open URL:',
                                      result.error
                                    );
                                    addNotification({
                                      type: 'error',
                                      title: 'Failed to open link',
                                      message:
                                        result.error ||
                                        'Could not open the URL in your browser.',
                                      autoClose: true,
                                      duration: 5000,
                                    });
                                  }
                                } catch (error) {
                                  console.error(
                                    'Error opening external URL:',
                                    error
                                  );
                                  addNotification({
                                    type: 'error',
                                    title: 'Failed to open link',
                                    message: 'An unexpected error occurred.',
                                    autoClose: true,
                                    duration: 5000,
                                  });
                                }
                              }
                            };

                            return (
                              <a href={href} onClick={handleClick} {...props}>
                                {children}
                              </a>
                            );
                          },
                        }}
                      >
                        {entry.parts[0].text}
                      </ReactMarkdown>
                      {entry.isStreaming && (
                        <span className={styles.cursor}>|</span>
                      )}
                    </>
                  ) : (
                    // ユーザーメッセージはプレーンテキストのまま
                    entry.parts[0].text
                  )}
                </div>
              </div>
            ))}
          </>
        )}
      </div>

      {/* エージェントステータス表示 */}
      {currentAgentStatus.agent && currentAgentStatus.status !== 'idle' && (
        <div className={styles.agentStatusBar}>
          <div className={styles.statusIndicator}>
            {currentAgentStatus.status === 'running' && (
              <>
                <div className={styles.spinner} />
                <span>
                  {getAgentDisplayName(currentAgentStatus.agent)}が動作中...
                </span>
              </>
            )}
            {currentAgentStatus.status === 'responding' && (
              <>
                <div className={styles.typingDots}>
                  <div className={styles.typingDot}></div>
                  <div className={styles.typingDot}></div>
                  <div className={styles.typingDot}></div>
                </div>
                <span>Supervisorが応答を生成中...</span>
              </>
            )}
          </div>
        </div>
      )}

      <div className={styles.inputArea}>
        <textarea
          ref={textareaRef}
          className={styles.inputBox}
          value={currentMessage}
          onChange={e => setCurrentMessage(e.target.value)}
          placeholder="Ask the AI agent about molecular design, calculations, or analysis..."
          onKeyDown={handleKeyDown}
          disabled={isLoading}
          rows={1}
        />
        <button
          className={styles.sendButton}
          onClick={handleSendMessage}
          disabled={isLoading || !currentMessage.trim()}
        >
          {isLoading ? (
            <>
              <svg
                width="16"
                height="16"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
                className={styles.loadingIndicator}
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"
                />
              </svg>
              Sending...
            </>
          ) : (
            <>
              <svg
                width="16"
                height="16"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
                style={{ transform: 'rotate(90deg)' }}
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth={2}
                  d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8"
                />
              </svg>
              Send
            </>
          )}
        </button>
      </div>

      {/* Confirmation Modal */}
      <ConfirmationModal
        isOpen={isConfirmationModalOpen}
        title="Confirm Destructive Action"
        message={
          confirmationRequest?.message || 'Are you sure you want to proceed?'
        }
        confirmButtonText="Confirm"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmAction}
        onCancel={handleCancelAction}
        isLoading={isExecutingAction}
      />

      {/* New Chat Confirmation Modal */}
      <ConfirmationModal
        isOpen={isNewChatConfirmationOpen}
        title="Start New Chat"
        message="Your current conversation will be saved to history. Would you like to start a new chat?"
        confirmButtonText="Start New Chat"
        cancelButtonText="Cancel"
        onConfirm={handleNewChatConfirm}
        onCancel={handleNewChatCancel}
        isLoading={false}
        variant="default"
      />
    </div>
  );
};
