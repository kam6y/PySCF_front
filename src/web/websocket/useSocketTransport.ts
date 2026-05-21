import { useCallback, useEffect, useRef, useState } from 'react';
import { io, Socket } from 'socket.io-client';

export interface SocketEventHandlers {
  onConnect?: (socket: Socket) => void;
  onDisconnect?: (reason: string) => void;
  onReconnect?: (socket: Socket, attemptNumber: number) => void;
  onConnectError?: (error: Error) => void;
  onError?: (errorData: unknown) => void;
  dataListeners?: Record<string, (data: unknown) => void>;
}

export interface UseSocketTransportOptions {
  handlers: SocketEventHandlers;
}

export interface UseSocketTransportReturn {
  isConnected: boolean;
  reconnect: () => void;
  disconnect: () => void;
  emit: (event: string, data?: unknown) => void;
}

export const useSocketTransport = ({
  handlers,
}: UseSocketTransportOptions): UseSocketTransportReturn => {
  const socketRef = useRef<Socket | null>(null);
  const disconnectTimerRef = useRef<NodeJS.Timeout | null>(null);
  const isConnectingRef = useRef<boolean>(false);
  const isMountedRef = useRef<boolean>(false);
  const [isConnected, setIsConnected] = useState<boolean>(false);

  // handlers を ref 化して安定化
  const onConnectRef = useRef(handlers.onConnect);
  onConnectRef.current = handlers.onConnect;
  const onDisconnectRef = useRef(handlers.onDisconnect);
  onDisconnectRef.current = handlers.onDisconnect;
  const onReconnectRef = useRef(handlers.onReconnect);
  onReconnectRef.current = handlers.onReconnect;
  const onConnectErrorRef = useRef(handlers.onConnectError);
  onConnectErrorRef.current = handlers.onConnectError;
  const onErrorRef = useRef(handlers.onError);
  onErrorRef.current = handlers.onError;
  const dataListenersRef = useRef(handlers.dataListeners);
  dataListenersRef.current = handlers.dataListeners;

  const clearDisconnectTimer = useCallback(() => {
    if (disconnectTimerRef.current) {
      clearTimeout(disconnectTimerRef.current);
      disconnectTimerRef.current = null;
    }
  }, []);

  const disconnect = useCallback(() => {
    clearDisconnectTimer();

    // 接続処理中の場合はフラグをリセット
    if (isConnectingRef.current) {
      isConnectingRef.current = false;
    }

    const socket = socketRef.current;
    if (socket) {
      try {
        socket.disconnect();
      } catch (error) {
        console.error('[UnifiedWebSocket] Error disconnecting:', error);
      }
      socketRef.current = null;
    }

    if (isMountedRef.current) {
      setIsConnected(false);
    }
  }, [clearDisconnectTimer]);

  const connect = useCallback(() => {
    // 既存の接続が有効な場合はスキップ
    if (socketRef.current?.connected) {
      console.log(
        '[UnifiedWebSocket] Already connected, skipping reconnection'
      );
      return;
    }

    // 接続処理中の場合はスキップ
    if (isConnectingRef.current) {
      console.log(
        '[UnifiedWebSocket] Connection already in progress, skipping'
      );
      return;
    }

    // 既存の接続があれば切断
    disconnect();

    // 接続処理を開始
    isConnectingRef.current = true;

    void (async () => {
      const port = window.flaskPort;
      if (!port) {
        console.error('[UnifiedWebSocket] Backend port not set. Cannot connect.');
        isConnectingRef.current = false;
        return;
      }
      const serverUrl = `http://127.0.0.1:${port}`;
      console.log(`[UnifiedWebSocket] Connecting to ${serverUrl}`);

      try {
        const authToken = await window.electronAPI?.getAuthToken?.();
        if (!isConnectingRef.current) {
          return;
        }

        const socket = io(serverUrl, {
          auth: { token: authToken ?? undefined },
          transports: ['websocket', 'polling'],
          timeout: 15000,
          reconnection: true,
          reconnectionDelay: 1000,
          reconnectionDelayMax: 5000,
          reconnectionAttempts: 5,
          randomizationFactor: 0.5,
          forceNew: true,
          upgrade: true,
          rememberUpgrade: false,
          autoConnect: true,
          withCredentials: false,
          extraHeaders: {
            Accept: 'application/json',
            'Cache-Control': 'no-cache',
          },
        });

        socketRef.current = socket;

        socket.on('connect', () => {
          socketRef.current = socket;
          isConnectingRef.current = false;
          if (isMountedRef.current) {
            setIsConnected(true);
          }
          onConnectRef.current?.(socket);
        });

        socket.on('disconnect', (reason: string) => {
          console.log(`[UnifiedWebSocket] Disconnected, reason: ${reason}`);
          if (isMountedRef.current) {
            setIsConnected(false);
          }
          onDisconnectRef.current?.(reason);
          if (!socket.active && socketRef.current === socket) {
            socketRef.current = null;
          }
        });

        socket.io.on('reconnect', attemptNumber => {
          const handler = onReconnectRef.current;
          if (!handler) return;
          Promise.resolve(handler(socket, attemptNumber)).catch(error => {
            console.error('[UnifiedWebSocket] Reconnect handler failed:', error);
          });
        });

        socket.io.on('reconnect_error', (error: Error) => {
          console.error('[UnifiedWebSocket] Reconnection failed:', error);
        });

        socket.on('connect_error', (error: Error) => {
          console.error('[UnifiedWebSocket] Connection error:', error);
          isConnectingRef.current = false;

          // ネットワーク切断やサーバー一時停止など、自動再接続されるエラーは通知しない
          const isTransientError =
            error.message.includes('websocket error') ||
            error.message.includes('502') ||
            error.message.includes('503');

          if (isTransientError) {
            console.log(
              '[UnifiedWebSocket] Transient error, auto-reconnecting...'
            );
            return;
          }

          onConnectErrorRef.current?.(error);
        });

        socket.on('error', (errorData: unknown) => {
          onErrorRef.current?.(errorData);
        });

        const initialListeners = dataListenersRef.current;
        if (initialListeners) {
          for (const eventName of Object.keys(initialListeners)) {
            socket.on(eventName, (data: unknown) => {
              const latestListener = dataListenersRef.current?.[eventName];
              latestListener?.(data);
            });
          }
        }
      } catch (error) {
        console.error('[UnifiedWebSocket] Failed to create connection:', error);
        isConnectingRef.current = false;

        const connectError =
          error instanceof Error ? error : new Error(String(error));
        onConnectErrorRef.current?.(connectError);
      }
    })();
  }, [disconnect]);

  // WebSocket接続の管理（StrictMode考慮）
  useEffect(() => {
    isMountedRef.current = true;

    // StrictModeでの二重実行を考慮して、少し遅延させる
    const connectTimer = setTimeout(() => {
      if (isMountedRef.current) {
        connect();
      }
    }, 50);

    return () => {
      isMountedRef.current = false;
      clearTimeout(connectTimer);

      // StrictModeでの即座の切断を防ぐため、少し遅延させる
      disconnectTimerRef.current = setTimeout(() => {
        if (!isMountedRef.current) {
          disconnect();
        }
      }, 100);
    };
  }, [connect, disconnect]);

  const emit = useCallback((event: string, data?: unknown) => {
    const socket = socketRef.current;
    if (!socket?.connected) return;
    try {
      socket.emit(event, data);
    } catch (error) {
      console.error(`[UnifiedWebSocket] Failed to emit ${event}:`, error);
    }
  }, []);

  return {
    isConnected,
    reconnect: connect,
    disconnect,
    emit,
  };
};
