import React, { useState, useCallback, useRef, useEffect } from 'react';
import { useAgentStore } from '../store/agentStore';
import { ConfirmationModal } from '../components/ConfirmationModal';
import { ChatMessage } from '../components/ChatMessage';
import { useAgentSession } from '../hooks/useAgentSession';
import { useAgentStream } from './agent/useAgentStream';
import styles from './AgentPage.module.css';

const samplePrompts = [
  {
    text: 'Run a DFT calculation for a water molecule',
    icon: (
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth={2}
        d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z"
      />
    ),
  },
  {
    text: 'Research the structure of photocatalysts',
    icon: (
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth={2}
        d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253"
      />
    ),
  },
  {
    text: 'Extend absorption wavelength based on benzene ring',
    icon: (
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth={2}
        d="M13 10V3L4 14h7v7l9-11h-7z"
      />
    ),
  },
  {
    text: 'Analyze the latest completed calculation',
    icon: (
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth={2}
        d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z"
      />
    ),
  },
];

export const AgentPage = React.memo(() => {
  // Zustandストアから会話履歴とエージェントステータスを取得
  const history = useAgentStore(state => state.history);
  const currentAgentStatus = useAgentStore(state => state.currentAgentStatus);

  const [currentMessage, setCurrentMessage] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const session = useAgentSession({ isLoading });
  const stream = useAgentStream({
    currentMessage,
    setCurrentMessage,
    isLoading,
    setIsLoading,
    activeSessionId: session.activeSessionId,
    setActiveSessionId: session.setActiveSessionId,
    prevSessionIdRef: session.prevSessionIdRef,
    queryClient: session.queryClient,
    createChatSession: session.createChatSession,
  });

  const { handleSendMessage, handleCancelMessage } = stream;

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

  // サンプルプロンプトクリック時のハンドラー
  const handlePromptClick = useCallback(
    (promptText: string) => {
      setCurrentMessage(promptText);
      // テキストエリアにフォーカスを移動
      if (textareaRef.current) {
        textareaRef.current.focus();
      }
    },
    [setCurrentMessage]
  );

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

      {/* Sample prompts */}
      <div className={styles.promptSuggestions}>
        {samplePrompts.map((prompt, index) => (
          <button
            key={index}
            className={styles.promptCard}
            onClick={() => handlePromptClick(prompt.text)}
          >
            <svg
              className={styles.promptIcon}
              fill="none"
              stroke="currentColor"
              viewBox="0 0 24 24"
            >
              {prompt.icon}
            </svg>
            <span className={styles.promptText}>{prompt.text}</span>
          </button>
        ))}
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
    return 'AI Assistant';
  };

  return (
    <div className={styles.agentPageContainer}>
      {/* Title and New Chat button container - only show when there's an active session */}
      {session.activeSessionId && session.sessionDetailData?.session && (
        <div className={styles.headerContainer}>
          {/* Conversation Title */}
          <div className={styles.titleContainer}>
            {session.isEditingTitle ? (
              <div className={styles.titleEditMode}>
                <input
                  type="text"
                  className={styles.titleInput}
                  value={session.editedTitle}
                  onChange={e => session.setEditedTitle(e.target.value)}
                  onKeyDown={session.handleTitleKeyDown}
                  autoFocus
                  maxLength={100}
                />
                <div className={styles.titleActionButtons}>
                  <button
                    className={styles.titleSaveButton}
                    onClick={session.handleSaveTitle}
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
                    onClick={session.handleCancelEditTitle}
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
                  {session.sessionDetailData.session.name}
                </h2>
                <button
                  className={styles.editButton}
                  onClick={session.handleStartEditTitle}
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
              onClick={session.handleNewChatRequest}
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
              <ChatMessage
                key={`${entry.tempId || index}-${entry.role}`}
                entry={entry}
                role={entry.role}
              />
            ))}
          </>
        )}
      </div>

      {/* AIステータス表示 */}
      {currentAgentStatus.agent && currentAgentStatus.status !== 'idle' && (
        <div className={styles.agentStatusBar}>
          <div className={styles.statusIndicator}>
            {currentAgentStatus.status === 'responding' && (
              <>
                <div className={styles.typingDots}>
                  <div className={styles.typingDot}></div>
                  <div className={styles.typingDot}></div>
                  <div className={styles.typingDot}></div>
                </div>
                <span>AI Assistantが応答を生成中...</span>
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
          rows={1}
        />
        <button
          className={styles.sendButton}
          onClick={isLoading ? handleCancelMessage : handleSendMessage}
          disabled={!isLoading && !currentMessage.trim()}
        >
          {isLoading ? (
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
                  d="M6 18L18 6M6 6l12 12"
                />
              </svg>
              Cancel
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

      {/* New Chat Confirmation Modal */}
      <ConfirmationModal
        isOpen={session.isNewChatConfirmationOpen}
        title="Start New Chat"
        message="Your current conversation will be saved to history. Would you like to start a new chat?"
        confirmButtonText="Start New Chat"
        cancelButtonText="Cancel"
        onConfirm={session.handleNewChatConfirm}
        onCancel={session.handleNewChatCancel}
        isLoading={false}
        variant="default"
      />
    </div>
  );
});
