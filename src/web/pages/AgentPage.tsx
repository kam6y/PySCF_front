import React, { useState, useCallback, useRef, useEffect } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeHighlight from 'rehype-highlight';
import { streamChatWithAgent, executeConfirmedAgentAction } from '../apiClient';
import { useNotificationStore } from '../store/notificationStore';
import { useAgentStore, ChatHistory, AgentStatus } from '../store/agentStore';
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
  const updateMessage = useAgentStore(state => state.updateMessage);
  const setHistory = useAgentStore(state => state.setHistory);
  const clearHistory = useAgentStore(state => state.clearHistory);
  const setAgentStatus = useAgentStore(state => state.setAgentStatus);

  const [currentMessage, setCurrentMessage] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const textareaRef = useRef<HTMLTextAreaElement>(null);
  const addNotification = useNotificationStore(state => state.addNotification);

  // 重複送信を防止するためのRef（setStateは非同期なので、useRefで同期的に管理）
  const isSendingRef = useRef(false);

  // Confirmation modal state
  const [confirmationRequest, setConfirmationRequest] =
    useState<ConfirmationRequest | null>(null);
  const [isConfirmationModalOpen, setIsConfirmationModalOpen] = useState(false);
  const [isExecutingAction, setIsExecutingAction] = useState(false);

  // Clear history confirmation modal state
  const [isClearConfirmationOpen, setIsClearConfirmationOpen] = useState(false);

  // Parse confirmation request from AI message
  const parseConfirmationRequest = useCallback(
    (text: string): ConfirmationRequest | null => {
      try {
        console.log(
          '[ConfirmationParser] Parsing text for confirmation request'
        );
        console.log('[ConfirmationParser] Text length:', text.length);
        console.log(
          '[ConfirmationParser] Text preview:',
          text.substring(0, 500)
        );

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
            console.log(
              '[ConfirmationParser] Valid confirmation request found in code block'
            );
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
          console.log(
            '[ConfirmationParser] Found',
            matches.length,
            'potential JSON matches'
          );

          // Try to parse each match
          for (const match of matches) {
            try {
              const parsed = JSON.parse(match);
              if (
                parsed.requires_confirmation === true &&
                parsed.action &&
                parsed.calculation_id
              ) {
                console.log(
                  '[ConfirmationParser] Valid confirmation request found'
                );
                return parsed as ConfirmationRequest;
              }
            } catch (parseError) {
              console.debug(
                '[ConfirmationParser] Failed to parse match:',
                parseError
              );
              continue;
            }
          }
        }

        // Try to find JSON object with a more lenient approach
        // Find all { ... } blocks and try to parse them
        const allJsonRegex = /\{(?:[^{}]|\{[^{}]*\})*\}/g;
        const allMatches = text.match(allJsonRegex);

        if (allMatches) {
          console.log(
            '[ConfirmationParser] Trying lenient parsing with',
            allMatches.length,
            'matches'
          );

          for (const match of allMatches) {
            try {
              const parsed = JSON.parse(match);
              if (
                parsed.requires_confirmation === true &&
                parsed.action &&
                parsed.calculation_id
              ) {
                console.log(
                  '[ConfirmationParser] Valid confirmation request found with lenient parsing'
                );
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

  // Handle clear history confirmation
  const handleClearHistoryRequest = useCallback(() => {
    setIsClearConfirmationOpen(true);
  }, []);

  const handleClearHistoryConfirm = useCallback(() => {
    clearHistory();
    setIsClearConfirmationOpen(false);
    addNotification({
      type: 'success',
      title: 'History Cleared',
      message: 'Conversation history has been cleared.',
      autoClose: true,
      duration: 3000,
    });
  }, [clearHistory, addNotification]);

  const handleClearHistoryCancel = useCallback(() => {
    setIsClearConfirmationOpen(false);
  }, []);

  // メッセージ送信処理
  const handleSendMessage = useCallback(() => {
    const trimmedMessage = currentMessage.trim();
    
    // 重複送信を防止（同期的にチェック）
    if (!trimmedMessage || isLoading || isSendingRef.current) {
      return;
    }

    // 送信中フラグを即座に設定（同期的）
    isSendingRef.current = true;

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
    addMessage(userMessage);
    addMessage(aiPlaceholder);
    const aiMessageIndex = history.length + 1; // AIメッセージのインデックス
    setCurrentMessage('');
    setIsLoading(true);
    setError(null);

    streamChatWithAgent(
      trimmedMessage,
      history, // ストリーム開始前の履歴を渡す
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
            status: status === 'running' || status === 'responding' ? status : 'idle'
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
        },
        onError: err => {
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
  }, [currentMessage, history, isLoading, addNotification, setAgentStatus]);

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

  // エージェント表示名を取得
  const getAgentDisplayName = (agent: string) => {
    const names: Record<string, string> = {
      'supervisor': 'Supervisor',
      'quantum_calculation_worker': 'Quantum Calculation Worker',
      'research_expert': 'Research Expert'
    };
    return names[agent] || agent;
  };

  return (
    <div className={styles.agentPageContainer}>
      {/* Clear button */}
      {history.length > 0 && (
        <div className={styles.clearButtonContainer}>
          <button
            className={styles.clearButton}
            onClick={handleClearHistoryRequest}
            title="Clear conversation history"
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
                d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16"
              />
            </svg>
            Clear History
          </button>
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
                            if (React.isValidElement(children) && children.props) {
                              const codeProps = children.props as any;
                              const className = codeProps.className;
                              
                              console.log('[Pre > Code Debug]', {
                                className,
                                hasClassName: !!className,
                                isOrbitalViewer: className?.includes('language-orbital-viewer'),
                              });
                              
                              if (className?.includes('language-orbital-viewer')) {
                                try {
                                  const codeContent = String(codeProps.children).replace(/\n$/, '');
                                  const params = parseOrbitalViewerParams(codeContent);
                                  
                                  console.log('[Orbital Viewer] Parsed params:', params);
                                  
                                  // Validate required parameters
                                  if (!params.calculation_id || params.orbital_index === undefined) {
                                    return (
                                      <div style={{
                                        padding: '12px',
                                        backgroundColor: '#fef2f2',
                                        border: '1px solid #fecaca',
                                        borderRadius: '8px',
                                        color: '#dc2626',
                                        margin: '12px 0'
                                      }}>
                                        ❌ Invalid orbital-viewer block: missing required parameters (calculation_id, orbital_index)
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
                                      onError={(error) => {
                                        console.error('Orbital viewer error:', error);
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
                                  console.error('Failed to render orbital viewer:', error);
                                  return (
                                    <div style={{
                                      padding: '12px',
                                      backgroundColor: '#fef2f2',
                                      border: '1px solid #fecaca',
                                      borderRadius: '8px',
                                      color: '#dc2626',
                                      margin: '12px 0'
                                    }}>
                                      ❌ Failed to render orbital viewer: {error instanceof Error ? error.message : String(error)}
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
                            if (!inline && className?.includes('language-orbital-viewer')) {
                              try {
                                const codeContent = String(children).replace(/\n$/, '');
                                const params = parseOrbitalViewerParams(codeContent);
                                
                                // Validate required parameters
                                if (!params.calculation_id || params.orbital_index === undefined) {
                                  return (
                                    <div style={{
                                      padding: '12px',
                                      backgroundColor: '#fef2f2',
                                      border: '1px solid #fecaca',
                                      borderRadius: '8px',
                                      color: '#dc2626',
                                      margin: '12px 0'
                                    }}>
                                      ❌ Invalid orbital-viewer block: missing required parameters (calculation_id, orbital_index)
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
                                    onError={(error) => {
                                      console.error('Orbital viewer error:', error);
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
                                console.error('Failed to render orbital viewer:', error);
                                return (
                                  <div style={{
                                    padding: '12px',
                                    backgroundColor: '#fef2f2',
                                    border: '1px solid #fecaca',
                                    borderRadius: '8px',
                                    color: '#dc2626',
                                    margin: '12px 0'
                                  }}>
                                    ❌ Failed to render orbital viewer: {error instanceof Error ? error.message : String(error)}
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
                              if (href && (href.startsWith('http://') || href.startsWith('https://'))) {
                                e.preventDefault();
                                try {
                                  const result = await window.electronAPI.openExternalUrl(href);
                                  if (!result.success) {
                                    console.error('Failed to open URL:', result.error);
                                    addNotification({
                                      type: 'error',
                                      title: 'Failed to open link',
                                      message: result.error || 'Could not open the URL in your browser.',
                                      autoClose: true,
                                      duration: 5000,
                                    });
                                  }
                                } catch (error) {
                                  console.error('Error opening external URL:', error);
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
                <span>{getAgentDisplayName(currentAgentStatus.agent)}が動作中...</span>
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

      {/* Clear History Confirmation Modal */}
      <ConfirmationModal
        isOpen={isClearConfirmationOpen}
        title="Clear Conversation History"
        message="Are you sure you want to clear all conversation history? This action cannot be undone."
        confirmButtonText="Clear"
        cancelButtonText="Cancel"
        onConfirm={handleClearHistoryConfirm}
        onCancel={handleClearHistoryCancel}
        isLoading={false}
      />
    </div>
  );
};
