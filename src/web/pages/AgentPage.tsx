import React, { useState, useCallback, useRef, useEffect } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeHighlight from 'rehype-highlight';
import { streamChatWithAgent, executeConfirmedAgentAction } from '../apiClient';
import { useNotificationStore } from '../store/notificationStore';
import { ConfirmationModal } from '../components/ConfirmationModal';
import styles from './AgentPage.module.css';

// 履歴の型定義（生成されたAPI型を使用）
type ChatHistory = {
  role: 'user' | 'model';
  parts: { text: string }[];
  isStreaming?: boolean; // ストリーミング中かを示すフラグ
  tempId?: number; // 一意識別子（ストリーミング中のメッセージ管理用）
};

// 確認リクエストの型定義
type ConfirmationRequest = {
  requires_confirmation: boolean;
  action: string;
  calculation_id: string;
  calculation_name: string;
  message: string;
};

export const AgentPage = () => {
  const [history, setHistory] = useState<ChatHistory[]>([]);
  const [currentMessage, setCurrentMessage] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const addNotification = useNotificationStore(state => state.addNotification);

  // Confirmation modal state
  const [confirmationRequest, setConfirmationRequest] =
    useState<ConfirmationRequest | null>(null);
  const [isConfirmationModalOpen, setIsConfirmationModalOpen] = useState(false);
  const [isExecutingAction, setIsExecutingAction] = useState(false);

  // Parse confirmation request from AI message
  const parseConfirmationRequest = useCallback(
    (text: string): ConfirmationRequest | null => {
      try {
        console.log('[ConfirmationParser] Parsing text for confirmation request');
        console.log('[ConfirmationParser] Text length:', text.length);
        console.log('[ConfirmationParser] Text preview:', text.substring(0, 500));

        // Try to extract JSON from markdown code blocks first (```json ... ```)
        const jsonBlockRegex = /```json\s*(\{[\s\S]*?\})\s*```/;
        const jsonBlockMatch = text.match(jsonBlockRegex);
        if (jsonBlockMatch?.[1]) {
          console.log('[ConfirmationParser] Found JSON in code block');
          const parsed = JSON.parse(jsonBlockMatch[1]);
          if (
            parsed.requires_confirmation === true &&
            parsed.action &&
            parsed.calculation_id
          ) {
            console.log('[ConfirmationParser] Valid confirmation request found in code block');
            return parsed as ConfirmationRequest;
          }
        }

        // Try to extract plain JSON object (more robust pattern)
        // This pattern finds JSON objects that contain "requires_confirmation": true
        // and handles multi-line structures with nested content
        const jsonPatternRegex = /\{[^{}]*"requires_confirmation"\s*:\s*true[^{}]*(?:\{[^{}]*\}[^{}]*)*\}/g;
        const matches = text.match(jsonPatternRegex);
        
        if (matches) {
          console.log('[ConfirmationParser] Found', matches.length, 'potential JSON matches');
          
          // Try to parse each match
          for (const match of matches) {
            try {
              const parsed = JSON.parse(match);
              if (
                parsed.requires_confirmation === true &&
                parsed.action &&
                parsed.calculation_id
              ) {
                console.log('[ConfirmationParser] Valid confirmation request found');
                return parsed as ConfirmationRequest;
              }
            } catch (parseError) {
              console.debug('[ConfirmationParser] Failed to parse match:', parseError);
              continue;
            }
          }
        }

        // Try to find JSON object with a more lenient approach
        // Find all { ... } blocks and try to parse them
        const allJsonRegex = /\{(?:[^{}]|\{[^{}]*\})*\}/g;
        const allMatches = text.match(allJsonRegex);
        
        if (allMatches) {
          console.log('[ConfirmationParser] Trying lenient parsing with', allMatches.length, 'matches');
          
          for (const match of allMatches) {
            try {
              const parsed = JSON.parse(match);
              if (
                parsed.requires_confirmation === true &&
                parsed.action &&
                parsed.calculation_id
              ) {
                console.log('[ConfirmationParser] Valid confirmation request found with lenient parsing');
                return parsed as ConfirmationRequest;
              }
            } catch (parseError) {
              // Silently skip invalid JSON
              continue;
            }
          }
        }

        console.log('[ConfirmationParser] No valid confirmation request found');
        return null;
      } catch (e) {
        console.error('[ConfirmationParser] Error parsing confirmation request:', e);
        return null;
      }
    },
    []
  );;

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
      setHistory(prev => [...prev, successMessage]);
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
      setHistory(prev => [...prev, errorMessage]);
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
    setHistory(prev => [...prev, cancelMessage]);
  }, []);

  // メッセージ送信処理
  const handleSendMessage = useCallback(() => {
    const trimmedMessage = currentMessage.trim();
    if (!trimmedMessage || isLoading) return;

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

    // ユーザーメッセージとAIのプレースホルダーを履歴に追加
    setHistory(prev => [...prev, userMessage, aiPlaceholder]);
    setCurrentMessage('');
    setIsLoading(true);
    setError(null);

    streamChatWithAgent(
      trimmedMessage,
      history, // ストリーム開始前の履歴を渡す
      {
        onMessage: chunk => {
          setHistory(prev => {
            const newHistory = [...prev];
            // 最後のメッセージを安全に取得し、tempIdで確認
            const lastIndex = newHistory.length - 1;
            const lastMessage = newHistory[lastIndex];

            if (
              lastMessage &&
              lastMessage.role === 'model' &&
              lastMessage.tempId === tempMessageId &&
              lastMessage.isStreaming
            ) {
              // 新しいオブジェクトを作成して状態を更新（不変性を保持）
              newHistory[lastIndex] = {
                ...lastMessage,
                parts: [{ text: lastMessage.parts[0].text + chunk }],
              };
            }
            return newHistory;
          });
        },
        onClose: () => {
          setIsLoading(false);
          setHistory(prev => {
            const newHistory = [...prev];
            const lastIndex = newHistory.length - 1;
            const lastMessage = newHistory[lastIndex];

            if (
              lastMessage &&
              lastMessage.role === 'model' &&
              lastMessage.tempId === tempMessageId
            ) {
              // ストリーミング終了時にtempIdを削除し、新しいオブジェクトを作成
              newHistory[lastIndex] = {
                role: lastMessage.role,
                parts: lastMessage.parts,
                isStreaming: false,
              };

              // Check for confirmation request in the final message
              const messageText = lastMessage.parts[0]?.text || '';
              const confirmation = parseConfirmationRequest(messageText);
              if (confirmation) {
                // Show confirmation modal
                setConfirmationRequest(confirmation);
                setIsConfirmationModalOpen(true);
              }
            }
            return newHistory;
          });
        },
        onError: err => {
          setIsLoading(false);
          setError(err.message);
          setHistory(prev => {
            const newHistory = [...prev];
            const lastIndex = newHistory.length - 1;
            const lastMessage = newHistory[lastIndex];

            if (
              lastMessage &&
              lastMessage.role === 'model' &&
              lastMessage.tempId === tempMessageId
            ) {
              // エラー時も新しいオブジェクトを作成
              newHistory[lastIndex] = {
                role: lastMessage.role,
                parts: [{ text: `Error: ${err.message}` }],
                isStreaming: false,
              };
            }
            return newHistory;
          });

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
  }, [currentMessage, history, isLoading, addNotification]);

  // Enter키で送信（Shift+Enterで改行）
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
        about your molecules or request assistance with computational chemistry
        tasks.
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

  return (
    <div className={styles.agentPageContainer}>
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
    </div>
  );
};
