import React, { useMemo } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeHighlight from 'rehype-highlight';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import 'katex/dist/katex.min.css';
import { ChatHistory } from '../store/agentStore';
import { useNotificationStore } from '../store/notificationStore';
import styles from './ChatMessage.module.css';

interface ChatMessageProps {
  entry: ChatHistory;
  role: 'user' | 'model';
}

export const ChatMessage: React.FC<ChatMessageProps> = React.memo(
  ({ entry, role }) => {
    const addNotification = useNotificationStore(
      state => state.addNotification
    );

    // ReactMarkdownのcomponentsをメモ化してリロードを防ぐ
    const markdownComponents = useMemo(
      () => ({
        a: ({ href, children, ...props }: any) => {
          // Handle external links - open in default browser
          const handleClick = async (e: React.MouseEvent) => {
            if (
              href &&
              (href.startsWith('http://') || href.startsWith('https://'))
            ) {
              e.preventDefault();
              try {
                const result = await window.electronAPI.openExternalUrl(href);
                if (!result.success) {
                  console.error('Failed to open URL:', result.error);
                  addNotification({
                    type: 'error',
                    title: 'Failed to open link',
                    message:
                      result.error || 'Could not open the URL in your browser.',
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
      }),
      [addNotification]
    );

    return (
      <div className={`${styles.chatMessage} ${styles[role]}`}>
        {role === 'model' ? (
          <>
            <ReactMarkdown
              remarkPlugins={[remarkGfm, remarkMath]}
              rehypePlugins={[
                rehypeHighlight,
                [
                  rehypeKatex,
                  {
                    throwOnError: false,
                    trust: false,
                    strict: 'warn',
                  },
                ],
              ]}
              disallowedElements={['script', 'iframe', 'object', 'embed']}
              unwrapDisallowed={true}
              className={styles.markdown}
              components={markdownComponents}
            >
              {entry.parts[0].text}
            </ReactMarkdown>
            {entry.isStreaming && <span className={styles.cursor}>|</span>}
          </>
        ) : (
          // ユーザーメッセージはプレーンテキストのまま
          entry.parts[0].text
        )}
      </div>
    );
  },
  (prevProps, nextProps) => {
    // カスタム比較関数: メッセージの内容とストリーミング状態が変わらなければ再レンダリングしない
    return (
      prevProps.entry.parts[0].text === nextProps.entry.parts[0].text &&
      prevProps.entry.isStreaming === nextProps.entry.isStreaming &&
      prevProps.role === nextProps.role
    );
  }
);
