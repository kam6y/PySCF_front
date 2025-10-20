import React from 'react';
import { useGetChatSessions } from '../hooks/useChatHistoryQueries';
import { useChatHistoryStore } from '../store/chatHistoryStore';
import styles from './ChatHistoryList.module.css';

interface ChatHistoryListProps {
  onSessionSelect: (sessionId: string) => void;
  onRequestDelete: (sessionId: string, sessionName: string) => void;
}

export const ChatHistoryList: React.FC<ChatHistoryListProps> = ({
  onSessionSelect,
  onRequestDelete
}) => {
  const { data, isLoading, error, refetch } = useGetChatSessions();
  const activeSessionId = useChatHistoryStore(state => state.activeSessionId);

  const handleSessionClick = (sessionId: string) => {
    onSessionSelect(sessionId);
  };

  const handleDeleteSession = (sessionId: string, sessionName: string, e: React.MouseEvent) => {
    e.stopPropagation();
    onRequestDelete(sessionId, sessionName);
  };

  if (isLoading) {
    return (
      <div className={styles.loadingContainer}>
        <p>Loading chat history...</p>
      </div>
    );
  }

  if (error) {
    return (
      <div className={styles.errorContainer}>
        <p>Error loading chat history</p>
        <p className={styles.errorMessage}>{(error as Error).message}</p>
        <button
          className={styles.retryButton}
          onClick={() => refetch()}
        >
          Retry
        </button>
      </div>
    );
  }

  if (!data || data.sessions.length === 0) {
    return (
      <div className={styles.emptyContainer}>
        <p>No chat history yet</p>
        <p className={styles.emptyHint}>Start a conversation to create your first chat</p>
      </div>
    );
  }

  return (
    <div className={styles.chatHistoryList}>
      {data.sessions.map(session => (
        <div
          key={session.id}
          className={`${styles.chatHistoryCard} ${
            session.id === activeSessionId ? styles.active : ''
          }`}
          onClick={() => handleSessionClick(session.id)}
        >
          <div className={styles.chatInfo}>
            <div className={styles.chatName}>{session.name}</div>
            <div className={styles.chatMeta}>
              <div className={styles.chatDate}>
                {new Date(session.updated_at).toLocaleString('en-US', {
                  year: 'numeric',
                  month: '2-digit',
                  day: '2-digit',
                  hour: '2-digit',
                  minute: '2-digit',
                })}
              </div>
              <div className={styles.chatMessageCount}>
                {session.message_count} {session.message_count === 1 ? 'message' : 'messages'}
              </div>
            </div>
          </div>
          <div className={styles.chatActions}>
            <button
              className={styles.deleteBtn}
              onClick={(e) => handleDeleteSession(session.id, session.name, e)}
              aria-label="Delete chat"
            >
              <svg
                width="14"
                height="14"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              >
                <polyline points="3,6 5,6 21,6"></polyline>
                <path d="m5,6 1,14 c0,1 1,2 2,2 h8 c1,0 2,-1 2,-2 l1,-14"></path>
                <path d="m10,11 v6"></path>
                <path d="m14,11 v6"></path>
                <path d="m7,6 V4 c0,-1 1,-2 2,-2 h6 c1,0 2,1 2,2 v2"></path>
              </svg>
            </button>
          </div>
        </div>
      ))}
    </div>
  );
};
