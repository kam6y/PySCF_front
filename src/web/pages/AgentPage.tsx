import React, { useState, useCallback, useRef, useEffect } from 'react';
import { streamChatWithAgent } from '../apiClient';
import { useNotificationStore } from '../store/notificationStore';
import styles from './AgentPage.module.css';

// 履歴の型定義（生成されたAPI型を使用）
type ChatHistory = {
  role: 'user' | 'model';
  parts: { text: string }[];
  isStreaming?: boolean; // ストリーミング中かを示すフラグ
};

export const AgentPage = () => {
  const [history, setHistory] = useState<ChatHistory[]>([]);
  const [currentMessage, setCurrentMessage] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const addNotification = useNotificationStore(state => state.addNotification);

  // メッセージ送信処理
  const handleSendMessage = useCallback(() => {
    const trimmedMessage = currentMessage.trim();
    if (!trimmedMessage || isLoading) return;

    const userMessage: ChatHistory = {
      role: 'user',
      parts: [{ text: trimmedMessage }],
    };

    // ユーザーメッセージとAIのプレースホルダーを履歴に追加
    setHistory(prev => [
      ...prev,
      userMessage,
      { role: 'model', parts: [{ text: '' }], isStreaming: true },
    ]);
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
            const lastMessage = newHistory[newHistory.length - 1];
            if (lastMessage && lastMessage.role === 'model') {
              lastMessage.parts[0].text += chunk;
            }
            return newHistory;
          });
        },
        onClose: () => {
          setIsLoading(false);
          setHistory(prev => {
            const newHistory = [...prev];
            const lastMessage = newHistory[newHistory.length - 1];
            if (lastMessage) {
              lastMessage.isStreaming = false;
            }
            return newHistory;
          });
        },
        onError: err => {
          setIsLoading(false);
          setError(err.message);
          setHistory(prev => {
            const newHistory = [...prev];
            const lastMessage = newHistory[newHistory.length - 1];
            if (lastMessage) {
              lastMessage.isStreaming = false;
              lastMessage.parts[0].text = `Error: ${err.message}`;
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
                  {entry.parts[0].text}
                  {entry.isStreaming && (
                    <span className={styles.cursor}>|</span>
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
    </div>
  );
};
