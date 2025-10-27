import React, { useMemo } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeHighlight from 'rehype-highlight';
import { ChatHistory } from '../store/agentStore';
import { InlineOrbitalViewer } from './InlineOrbitalViewer';
import { InlineIRSpectrumViewer } from './InlineIRSpectrumViewer';
import { InlineMullikenChargeViewer } from './InlineMullikenChargeViewer';
import { LazyViewer } from './LazyViewer';
import { useNotificationStore } from '../store/notificationStore';
import styles from './ChatMessage.module.css';

// パラメータパース関数（汎用的なコードブロックパラメータ用）
const parseCodeBlockParams = (code: string): Record<string, any> => {
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
        pre: ({ children, ...props }: any) => {
          // Check if this pre contains an orbital-viewer code block
          if (React.isValidElement(children) && children.props) {
            const codeProps = children.props as any;
            const className = codeProps.className;

            // Handle orbital-viewer code block
            if (className?.includes('language-orbital-viewer')) {
              try {
                const codeContent = String(codeProps.children).replace(
                  /\n$/,
                  ''
                );
                const params = parseCodeBlockParams(codeContent);

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
                      ❌ Invalid orbital-viewer block: missing required
                      parameters (calculation_id, orbital_index)
                    </div>
                  );
                }

                // Render InlineOrbitalViewer component
                return (
                  <LazyViewer>
                    <InlineOrbitalViewer
                      calculation_id={params.calculation_id}
                      orbital_index={params.orbital_index}
                      grid_size={params.grid_size}
                      isovalue_pos={params.isovalue_pos}
                      isovalue_neg={params.isovalue_neg}
                      onError={error => {
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
                  </LazyViewer>
                );
              } catch (error) {
                console.error('Failed to render orbital viewer:', error);
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
                    {error instanceof Error ? error.message : String(error)}
                  </div>
                );
              }
            }

            // Handle ir-spectrum code block
            if (className?.includes('language-ir-spectrum')) {
              try {
                const codeContent = String(codeProps.children).replace(
                  /\n$/,
                  ''
                );
                const params = parseCodeBlockParams(codeContent);

                // Validate required parameters
                if (!params.calculation_id) {
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
                      ❌ Invalid ir-spectrum block: missing required parameter
                      (calculation_id)
                    </div>
                  );
                }

                // Render InlineIRSpectrumViewer component
                return (
                  <LazyViewer>
                    <InlineIRSpectrumViewer
                      calculation_id={params.calculation_id}
                      broadening_fwhm={params.broadening_fwhm}
                      x_min={params.x_min}
                      x_max={params.x_max}
                      show_peaks={params.show_peaks}
                      height={params.height}
                      onError={error => {
                        console.error('IR spectrum viewer error:', error);
                        addNotification({
                          type: 'error',
                          title: 'IR Spectrum Viewer Error',
                          message: error,
                          autoClose: false,
                          duration: 0,
                        });
                      }}
                    />
                  </LazyViewer>
                );
              } catch (error) {
                console.error('Failed to render IR spectrum viewer:', error);
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
                    ❌ Failed to render IR spectrum viewer:{' '}
                    {error instanceof Error ? error.message : String(error)}
                  </div>
                );
              }
            }

            // Handle mulliken-charges code block
            if (className?.includes('language-mulliken-charges')) {
              try {
                const codeContent = String(codeProps.children).replace(
                  /\n$/,
                  ''
                );
                const params = parseCodeBlockParams(codeContent);

                // Validate required parameters
                if (!params.calculation_id) {
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
                      ❌ Invalid mulliken-charges block: missing required
                      parameter (calculation_id)
                    </div>
                  );
                }

                // Render InlineMullikenChargeViewer component
                return (
                  <LazyViewer>
                    <InlineMullikenChargeViewer
                      calculation_id={params.calculation_id}
                      height={params.height}
                      onError={error => {
                        console.error('Mulliken charges viewer error:', error);
                        addNotification({
                          type: 'error',
                          title: 'Mulliken Charges Viewer Error',
                          message: error,
                          autoClose: false,
                          duration: 0,
                        });
                      }}
                    />
                  </LazyViewer>
                );
              } catch (error) {
                console.error(
                  'Failed to render Mulliken charges viewer:',
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
                    ❌ Failed to render Mulliken charges viewer:{' '}
                    {error instanceof Error ? error.message : String(error)}
                  </div>
                );
              }
            }
          }

          // Default pre rendering
          return <pre {...props}>{children}</pre>;
        },
        code: ({ className, children, ...props }: any) => {
          // Check if this is a special code block
          const inline = props.inline;

          // Handle orbital-viewer code block
          if (!inline && className?.includes('language-orbital-viewer')) {
            try {
              const codeContent = String(children).replace(/\n$/, '');
              const params = parseCodeBlockParams(codeContent);

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
                    ❌ Invalid orbital-viewer block: missing required parameters
                    (calculation_id, orbital_index)
                  </div>
                );
              }

              // Render InlineOrbitalViewer component
              return (
                <LazyViewer>
                  <InlineOrbitalViewer
                    calculation_id={params.calculation_id}
                    orbital_index={params.orbital_index}
                    grid_size={params.grid_size}
                    isovalue_pos={params.isovalue_pos}
                    isovalue_neg={params.isovalue_neg}
                    onError={error => {
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
                </LazyViewer>
              );
            } catch (error) {
              console.error('Failed to render orbital viewer:', error);
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
                  {error instanceof Error ? error.message : String(error)}
                </div>
              );
            }
          }

          // Handle ir-spectrum code block
          if (!inline && className?.includes('language-ir-spectrum')) {
            try {
              const codeContent = String(children).replace(/\n$/, '');
              const params = parseCodeBlockParams(codeContent);

              // Validate required parameters
              if (!params.calculation_id) {
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
                    ❌ Invalid ir-spectrum block: missing required parameter
                    (calculation_id)
                  </div>
                );
              }

              // Render InlineIRSpectrumViewer component
              return (
                <LazyViewer>
                  <InlineIRSpectrumViewer
                    calculation_id={params.calculation_id}
                    broadening_fwhm={params.broadening_fwhm}
                    x_min={params.x_min}
                    x_max={params.x_max}
                    show_peaks={params.show_peaks}
                    height={params.height}
                    onError={error => {
                      console.error('IR spectrum viewer error:', error);
                      addNotification({
                        type: 'error',
                        title: 'IR Spectrum Viewer Error',
                        message: error,
                        autoClose: false,
                        duration: 0,
                      });
                    }}
                  />
                </LazyViewer>
              );
            } catch (error) {
              console.error('Failed to render IR spectrum viewer:', error);
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
                  ❌ Failed to render IR spectrum viewer:{' '}
                  {error instanceof Error ? error.message : String(error)}
                </div>
              );
            }
          }

          // Handle mulliken-charges code block
          if (!inline && className?.includes('language-mulliken-charges')) {
            try {
              const codeContent = String(children).replace(/\n$/, '');
              const params = parseCodeBlockParams(codeContent);

              // Validate required parameters
              if (!params.calculation_id) {
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
                    ❌ Invalid mulliken-charges block: missing required
                    parameter (calculation_id)
                  </div>
                );
              }

              // Render InlineMullikenChargeViewer component
              return (
                <LazyViewer>
                  <InlineMullikenChargeViewer
                    calculation_id={params.calculation_id}
                    height={params.height}
                    onError={error => {
                      console.error('Mulliken charges viewer error:', error);
                      addNotification({
                        type: 'error',
                        title: 'Mulliken Charges Viewer Error',
                        message: error,
                        autoClose: false,
                        duration: 0,
                      });
                    }}
                  />
                </LazyViewer>
              );
            } catch (error) {
              console.error('Failed to render Mulliken charges viewer:', error);
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
                  ❌ Failed to render Mulliken charges viewer:{' '}
                  {error instanceof Error ? error.message : String(error)}
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
        <div className={styles.messageContent}>
          {role === 'model' ? (
            <>
              <ReactMarkdown
                remarkPlugins={[remarkGfm]}
                rehypePlugins={[rehypeHighlight]}
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
