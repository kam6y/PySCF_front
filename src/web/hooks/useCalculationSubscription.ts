// src/web/hooks/useCalculationSubscription.ts

import { useEffect, useRef, useCallback } from 'react';
import { io, Socket } from 'socket.io-client';
import { CalculationInstance } from '../types/api-types';
import { showErrorNotification } from '../store/notificationStore';

export interface UseCalculationSubscriptionOptions {
  calculationId: string | null;
  status?: 'pending' | 'running' | 'completed' | 'error' | 'waiting';
  onUpdate: (calculation: CalculationInstance) => void;
  onError?: (error: string) => void;
}

export const useCalculationSubscription = ({
  calculationId,
  status,
  onUpdate,
  onError,
}: UseCalculationSubscriptionOptions) => {
  const socketRef = useRef<Socket | null>(null);
  const onUpdateRef = useRef(onUpdate);
  const onErrorRef = useRef(onError);
  const currentCalculationIdRef = useRef<string | null>(null);
  const disconnectTimerRef = useRef<NodeJS.Timeout | null>(null);

  // コールバック関数の参照を最新に保つ
  useEffect(() => {
    onUpdateRef.current = onUpdate;
    onErrorRef.current = onError;
  }, [onUpdate, onError]);

  // 遅延切断タイマーをクリアするヘルパー関数
  const clearDisconnectTimer = useCallback(() => {
    if (disconnectTimerRef.current) {
      clearTimeout(disconnectTimerRef.current);
      disconnectTimerRef.current = null;
    }
  }, []);

  const disconnect = useCallback(() => {
    clearDisconnectTimer();
    if (socketRef.current) {
      const socket = socketRef.current;
      try {
        // Leave calculation room if we have an active calculation
        if (currentCalculationIdRef.current) {
          socket.emit('leave_calculation', {
            calculation_id: currentCalculationIdRef.current,
          });
        }

        // Disconnect from Socket.IO server
        socket.disconnect();
      } catch (error) {
        console.error('Error disconnecting Socket.IO:', error);
      }
      socketRef.current = null;
    }
    currentCalculationIdRef.current = null;
  }, [clearDisconnectTimer]);

  const connect = useCallback(
    (calcId: string) => {
      // 一時的IDの場合は接続を試行しない
      if (calcId.startsWith('new-calculation-')) {
        return;
      }

      // 既に同じIDで接続中の場合は何もしない
      if (
        currentCalculationIdRef.current === calcId &&
        socketRef.current &&
        socketRef.current.connected
      ) {
        return;
      }

      // 既存接続をクリーンアップ
      disconnect();

      try {
        // Get the Flask-SocketIO server URL (always HTTP in Electron)
        const serverUrl = `http://127.0.0.1:${window.flaskPort}`;

        // Create Socket.IO connection
        const socket = io(serverUrl, {
          transports: ['websocket', 'polling'],
          timeout: 5000,
        });

        socketRef.current = socket;
        currentCalculationIdRef.current = calcId;

        socket.on('connect', () => {
          console.log('Socket.IO connected, joining calculation room');
          // Join the calculation room to receive updates
          socket.emit('join_calculation', { calculation_id: calcId });
        });

        socket.on(
          'calculation_update',
          (updatedCalculation: CalculationInstance) => {
            try {
              if (onUpdateRef.current) {
                onUpdateRef.current(updatedCalculation);
              }

              // 計算が失敗した場合、詳細エラーをトースト通知で表示
              if (updatedCalculation.status === 'error') {
                const calculationName =
                  updatedCalculation.name || updatedCalculation.id;

                // 複数のエラー情報ソースを確認（バックエンドとの不整合に対応）
                const errorMessage =
                  updatedCalculation.errorMessage ||
                  (updatedCalculation as any).error || // バックエンドが設定するフィールド
                  updatedCalculation.results?.error ||
                  'Detailed error information is not available.';

                showErrorNotification(
                  `Calculation "${calculationName}" failed`,
                  errorMessage
                );
              }

              // 計算が完了または失敗した場合、少し待ってから接続を切断
              if (
                updatedCalculation.status === 'completed' ||
                updatedCalculation.status === 'error'
              ) {
                // 既存のタイマーをクリア
                clearDisconnectTimer();

                // 3秒後に切断
                disconnectTimerRef.current = setTimeout(() => {
                  disconnect();
                }, 3000);
              }
            } catch (e) {
              console.error('Failed to process calculation update:', e);
              if (onErrorRef.current) {
                onErrorRef.current(
                  'Unable to parse response from server. Calculation status updates may be temporarily stopped.'
                );
              }
            }
          }
        );

        socket.on('error', (errorData: any) => {
          console.error('Socket.IO error:', errorData);
          if (onErrorRef.current) {
            const errorMessage =
              errorData?.error || 'A problem occurred while monitoring calculation status.';
            onErrorRef.current(errorMessage);
          }
        });

        socket.on('connect_error', (error: Error) => {
          console.error('Socket.IO connection error:', error);
          if (onErrorRef.current) {
            let errorMessage =
              'Failed to connect to server. Please check your network connection.';

            if (error.message.includes('timeout')) {
              errorMessage =
                'Connection timed out. The server may not be running.';
            } else if (error.message.includes('xhr poll error')) {
              errorMessage =
                'Communication with server was interrupted. Calculation continues, but status updates have stopped.';
            }

            onErrorRef.current(errorMessage);
          }
        });

        socket.on('disconnect', (reason: string) => {
          console.log(
            `Socket.IO disconnected for calculation ${calcId}, reason:`,
            reason
          );
          if (socketRef.current === socket) {
            socketRef.current = null;
          }
          if (currentCalculationIdRef.current === calcId) {
            currentCalculationIdRef.current = null;
          }

          // Auto-reconnect for certain disconnect reasons
          if (reason === 'io server disconnect') {
            // Server initiated disconnect, likely planned shutdown
            if (onErrorRef.current) {
              onErrorRef.current(
                'Connection to server was disconnected. Please refresh the page to check status.'
              );
            }
          }
        });
      } catch (error) {
        console.error('Failed to create Socket.IO connection:', error);
        if (onErrorRef.current) {
          let errorMessage = 'Failed to start real-time monitoring.';

          if (error instanceof Error) {
            if (error.message.includes('network')) {
              errorMessage =
                'Unable to connect to server due to network error.';
            } else if (error.message.includes('security')) {
              errorMessage =
                'Connection blocked by security settings. HTTPS may be required.';
            }
          }

          onErrorRef.current(errorMessage);
        }
      }
    },
    [disconnect, clearDisconnectTimer]
  );

  useEffect(() => {
    // 計算IDが有効で一時的IDでない場合に接続
    // running/pending状態では即座に接続
    // completed/error状態でも接続を維持（自動切断タイマーで後で切断）
    if (calculationId && !calculationId.startsWith('new-calculation-')) {
      // running/pending状態では即座に接続
      if (status === 'running' || status === 'pending') {
        connect(calculationId);
      } else if (status === 'completed' || status === 'error') {
        // completed/error状態でも、まだ接続していない場合は接続する
        // （既に接続している場合は、イベントハンドラで自動切断タイマーが設定される）
        if (!socketRef.current || !socketRef.current.connected) {
          connect(calculationId);
        }
      } else {
        // その他の状態（unknown等）では切断
        disconnect();
      }
    } else {
      disconnect();
    }

    // クリーンアップ
    return () => {
      disconnect();
    };
  }, [calculationId, status, connect, disconnect]);
};
