import React, { useState, useCallback, useRef, useEffect } from 'react';
import { useMutation } from '@tanstack/react-query';
import { chatWithAgent } from '../apiClient';
import { useNotificationStore } from '../store/notificationStore';
import styles from './AgentPage.module.css';

// 履歴の型定義（生成されたAPI型を使用）
type ChatHistory = {
  role: 'user' | 'model';
  parts: { text: string }[];
};

export const AgentPage = () => {
  const [history, setHistory] = useState<ChatHistory[]>([]);
  const [currentMessage, setCurrentMessage] = useState('');
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const addNotification = useNotificationStore(state => state.addNotification);

  // チャットAPI呼び出し用mutation
  const chatMutation = useMutation({
    mutationFn: ({ message, history }: { message: string; history: ChatHistory[] }) =>
      chatWithAgent(message, history),
    onMutate: ({ message }) => {
      // 楽観的更新: ユーザーメッセージを即座に追加
      const newUserMessage: ChatHistory = {
        role: 'user',
        parts: [{ text: message }],
      };
      setHistory(prev => [...prev, newUserMessage]);
      setCurrentMessage('');
    },
    onSuccess: (response) => {
      // AIの応答を履歴に追加
      const agentReply: ChatHistory = {
        role: 'model',
        parts: [{ text: response.reply }],
      };
      setHistory(prev => [...prev, agentReply]);
    },
    onError: (error) => {
      console.error('Failed to send message to agent:', error);
      
      // エラーメッセージを履歴に追加
      const errorReply: ChatHistory = {
        role: 'model',
        parts: [{ text: 'Sorry, an error occurred while processing your message. Please try again.' }],
      };
      setHistory(prev => [...prev, errorReply]);
      
      // 通知でエラーを表示
      addNotification({
        type: 'error',
        title: 'AI Agent Error',
        message: 'Failed to send message to AI agent. Please try again.',
        autoClose: false,
        duration: 0,
      });
    },
  });

  // メッセージ送信処理
  const handleSendMessage = useCallback(() => {
    const trimmedMessage = currentMessage.trim();
    if (!trimmedMessage || chatMutation.isPending) return;

    chatMutation.mutate({
      message: trimmedMessage,
      history: history,
    });
  }, [currentMessage, history, chatMutation]);

  // Enter키で送信（Shift+Enterで改行）
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSendMessage();
    }
  }, [handleSendMessage]);

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
  }, [history, chatMutation.isPending]);

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
        Start a conversation with the AI agent to get help with molecular design,
        quantum chemistry calculations, and analysis. Ask questions about your molecules
        or request assistance with computational chemistry tasks.
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
                </div>
              </div>
            ))}
            {chatMutation.isPending && renderTypingIndicator()}
          </>
        )}
      </div>

      <div className={styles.inputArea}>
        <textarea
          ref={textareaRef}
          className={styles.inputBox}
          value={currentMessage}
          onChange={(e) => setCurrentMessage(e.target.value)}
          placeholder="Ask the AI agent about molecular design, calculations, or analysis..."
          onKeyDown={handleKeyDown}
          disabled={chatMutation.isPending}
          rows={1}
        />
        <button
          className={styles.sendButton}
          onClick={handleSendMessage}
          disabled={chatMutation.isPending || !currentMessage.trim()}
        >
          {chatMutation.isPending ? (
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