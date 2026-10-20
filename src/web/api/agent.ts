import { fetchEventSource } from '@microsoft/fetch-event-source';
import { components } from '../types/generated-api';
import { getApiBaseUrl } from './core';

type AgentChatRequest = components['schemas']['AgentChatRequest'];

const isAbortError = (error: unknown, signal: AbortSignal): boolean => {
  if (signal.aborted) {
    return true;
  }

  if (typeof error !== 'object' || error === null || !('name' in error)) {
    return false;
  }

  return error.name === 'AbortError';
};

export const streamChatWithAgent = (
  message: string,
  history: AgentChatRequest['history'],
  sessionId: string | null,
  callbacks: {
    onMessage: (chunk: string) => void;
    onClose: () => void;
    onError: (error: Error) => void;
    onAgentStatus?: (
      status: 'running' | 'completed' | 'responding',
      agent: string
    ) => void;
  }
) => {
  const ctrl = new AbortController();
  let isStreamClosed = false;

  (async () => {
    try {
      const headers: HeadersInit = {
        'Content-Type': 'application/json',
        Accept: 'text/event-stream',
      };

      await fetchEventSource(`${getApiBaseUrl()}/api/agent/chat`, {
        method: 'POST',
        headers,
        body: JSON.stringify({ message, history, session_id: sessionId }),
        signal: ctrl.signal,

        onopen: async response => {
          if (!response.ok) {
            const errorText = await response.text();
            if (!isStreamClosed) {
              isStreamClosed = true;
              callbacks.onError(
                new Error(`Failed to connect: ${response.status} ${errorText}`)
              );
            }
            ctrl.abort();
          }
        },

        onmessage(event) {
          if (isStreamClosed) {
            return;
          }

          try {
            const parsedData = JSON.parse(event.data);

            if (parsedData.type === 'chunk' && parsedData.payload?.text) {
              callbacks.onMessage(parsedData.payload.text);
            } else if (
              parsedData.type === 'agent_status' &&
              parsedData.payload
            ) {
              if (callbacks.onAgentStatus) {
                callbacks.onAgentStatus(
                  parsedData.payload.status,
                  parsedData.payload.agent
                );
              }
            } else if (parsedData.type === 'done') {
              if (!isStreamClosed) {
                isStreamClosed = true;
                callbacks.onClose();
              }
              ctrl.abort();
            } else if (parsedData.type === 'error') {
              if (!isStreamClosed) {
                isStreamClosed = true;
                callbacks.onError(
                  new Error(
                    parsedData.payload?.message ||
                      'An unknown stream error occurred.'
                  )
                );
              }
              ctrl.abort();
            }
          } catch (e) {
            if (!isStreamClosed) {
              isStreamClosed = true;
              callbacks.onError(
                new Error('Failed to parse message from stream.')
              );
              ctrl.abort();
            }
          }
        },

        onclose() {
          if (!isStreamClosed) {
            isStreamClosed = true;
            callbacks.onClose();
          }
        },

        onerror(err) {
          if (isAbortError(err, ctrl.signal)) {
            isStreamClosed = true;
            return;
          }

          if (!isStreamClosed) {
            isStreamClosed = true;
            callbacks.onError(
              err instanceof Error ? err : new Error(String(err))
            );
          }
          throw err;
        },
      });
    } catch (error) {
      if (isAbortError(error, ctrl.signal)) {
        isStreamClosed = true;
        return;
      }

      if (!isStreamClosed) {
        isStreamClosed = true;
        callbacks.onError(
          error instanceof Error ? error : new Error(String(error))
        );
      }
    }
  })();

  return () => {
    isStreamClosed = true;
    ctrl.abort();
  };
};
