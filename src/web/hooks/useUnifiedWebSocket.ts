import { useEffect, useRef, useCallback } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import { io, Socket } from 'socket.io-client';
import { CalculationInstance } from '../types/api-types';
import {
  showErrorNotification,
  showSuccessNotification,
  showInfoNotification,
  showResourceInsufficientErrorNotification,
} from '../store/notificationStore';

export interface UseUnifiedWebSocketOptions {
  activeCalculationId: string | null;
}

/**
 * 統合WebSocketフック - すべてのWebSocket通信を一元管理
 *
 * 機能:
 * - グローバル計算監視 (global_updates ルーム)
 * - アクティブ計算の個別監視 (calculation_{id} ルーム)
 * - 統一されたReact Queryキャッシュ更新
 * - 重複更新防止
 * - 統一された通知システム
 */
export const useUnifiedWebSocket = ({
  activeCalculationId,
}: UseUnifiedWebSocketOptions) => {
  const queryClient = useQueryClient();
  const socketRef = useRef<Socket | null>(null);
  const lastUpdateTimestamps = useRef<Map<string, number>>(new Map());
  const currentActiveCalculationId = useRef<string | null>(null);
  const disconnectTimerRef = useRef<NodeJS.Timeout | null>(null);

  // 重複更新防止チェック
  const isDuplicateUpdate = useCallback((calculationId: string): boolean => {
    const currentTimestamp = Date.now();
    const lastUpdateTime = lastUpdateTimestamps.current.get(calculationId);

    if (lastUpdateTime && currentTimestamp - lastUpdateTime < 100) {
      console.log(
        `[UnifiedWebSocket] Skipping duplicate update for calculation ${calculationId} (within 100ms)`
      );
      return true;
    }

    lastUpdateTimestamps.current.set(calculationId, currentTimestamp);
    return false;
  }, []);

  // 統一されたキャッシュ更新ロジック
  const updateCalculationCache = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;

      try {
        // 1. 個別計算詳細のキャッシュを更新
        queryClient.setQueryData(['calculation', calculationId], {
          calculation: updatedCalculation,
        });

        // 2. 計算リストのキャッシュを更新（楽観的更新）
        queryClient.setQueryData(['calculations'], (oldData: any) => {
          if (!oldData?.calculations) {
            console.warn(
              '[UnifiedWebSocket] No calculations data found for cache update'
            );
            return oldData;
          }

          const calculationIndex = oldData.calculations.findIndex(
            (calc: any) => calc.id === calculationId
          );

          let updatedCalculations;
          if (calculationIndex === -1) {
            // 計算がリストに存在しない場合は新規追加
            const newCalculationItem = {
              id: calculationId,
              name: updatedCalculation.name || calculationId,
              status: updatedCalculation.status,
              date:
                updatedCalculation.updatedAt ||
                updatedCalculation.createdAt ||
                new Date().toISOString(),
            };

            updatedCalculations = [newCalculationItem, ...oldData.calculations];
            console.log(
              `[UnifiedWebSocket] Added new calculation ${calculationId} to cache with status ${updatedCalculation.status}`
            );
          } else {
            // 既存の計算の場合は更新
            updatedCalculations = oldData.calculations.map((calc: any) =>
              calc.id === calculationId
                ? {
                    ...calc,
                    status: updatedCalculation.status,
                    date: updatedCalculation.updatedAt || calc.date,
                    name: updatedCalculation.name || calc.name,
                  }
                : calc
            );
            console.log(
              `[UnifiedWebSocket] Updated existing calculation ${calculationId} with status ${updatedCalculation.status}`
            );
          }

          return {
            ...oldData,
            calculations: updatedCalculations,
          };
        });
      } catch (error) {
        console.error(
          `[UnifiedWebSocket] Failed to update cache for calculation ${calculationId}:`,
          error
        );
        // フォールバック：エラーが発生した場合は関連クエリを無効化
        queryClient.invalidateQueries({
          queryKey: ['calculation', calculationId],
        });
        queryClient.invalidateQueries({ queryKey: ['calculations'] });
      }
    },
    [queryClient]
  );

  // 通知処理
  const handleNotifications = useCallback(
    (
      updatedCalculation: CalculationInstance,
      previousStatus: string | undefined
    ) => {
      const calculationId = updatedCalculation.id;
      const isActiveCalculation =
        calculationId === currentActiveCalculationId.current;
      const molecularName = updatedCalculation.name || 'Unknown';
      const currentErrorMessage =
        updatedCalculation.errorMessage || updatedCalculation.results?.error;

      // 計算完了通知: running -> completed
      if (
        previousStatus === 'running' &&
        updatedCalculation.status === 'completed'
      ) {
        if (!isActiveCalculation) {
          showSuccessNotification(
            'Calculation completed',
            `Calculation for ${molecularName} has been completed.`,
            calculationId
          );
        }
      }

      // 計算エラー通知: running -> error または pending -> error
      if (
        (previousStatus === 'running' || previousStatus === 'pending') &&
        updatedCalculation.status === 'error'
      ) {
        const errorMessage =
          currentErrorMessage || 'Detailed error information is not available.';

        // リソース不足エラーの判定
        const isResourceInsufficientError =
          errorMessage &&
          (errorMessage.toLowerCase().includes('cpu usage') ||
            errorMessage.toLowerCase().includes('memory usage') ||
            errorMessage.toLowerCase().includes('system cpu usage') ||
            errorMessage.toLowerCase().includes('system memory usage') ||
            errorMessage.toLowerCase().includes('no active calculations'));

        const title =
          previousStatus === 'running'
            ? `Calculation "${molecularName}" failed`
            : `Calculation "${molecularName}" failed to start`;

        if (isResourceInsufficientError) {
          showResourceInsufficientErrorNotification(
            errorMessage,
            calculationId
          );
        } else {
          showErrorNotification(title, errorMessage, calculationId);
        }
      }

      // 計算開始通知: pending -> running
      if (
        previousStatus === 'pending' &&
        updatedCalculation.status === 'running'
      ) {
        if (!isActiveCalculation) {
          showInfoNotification(
            'Calculation started',
            `Calculation for ${molecularName} has been started.`,
            calculationId
          );
        }
      }

      // 待機からの開始通知: waiting -> running
      if (
        previousStatus === 'waiting' &&
        updatedCalculation.status === 'running'
      ) {
        if (!isActiveCalculation) {
          showInfoNotification(
            'Queued calculation started',
            `Calculation for ${molecularName} has started from queued status.`,
            calculationId
          );
        }
      }

      // 待機状態への移行通知: pending -> waiting
      if (
        previousStatus === 'pending' &&
        updatedCalculation.status === 'waiting'
      ) {
        if (!isActiveCalculation) {
          const waitingReason =
            updatedCalculation.waitingReason ||
            'Please wait for available system resources or execution slots.';
          showInfoNotification(
            'Calculation is waiting',
            `Calculation for ${molecularName} is currently waiting. Reason: ${waitingReason}`,
            calculationId
          );
        }
      }
    },
    []
  );

  // 統一された計算更新処理
  const handleCalculationUpdate = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;

      if (!updatedCalculation || !calculationId) {
        console.warn(
          '[UnifiedWebSocket] Received invalid calculation update data:',
          updatedCalculation
        );
        return;
      }

      // 重複更新チェック
      if (isDuplicateUpdate(calculationId)) {
        return;
      }

      console.log(
        `[UnifiedWebSocket] Processing update for calculation ${calculationId}: ${updatedCalculation.status}`
      );

      // 前のステータスを取得
      const previousDetailData = queryClient.getQueryData([
        'calculation',
        calculationId,
      ]) as { calculation: CalculationInstance } | undefined;
      const previousStatus = previousDetailData?.calculation?.status;

      // 前のエラーメッセージを取得
      const previousErrorMessage =
        previousDetailData?.calculation?.errorMessage ||
        previousDetailData?.calculation?.results?.error;
      const currentErrorMessage =
        updatedCalculation.errorMessage || updatedCalculation.results?.error;

      // 変更があった場合のみ処理
      const hasStatusChanged = previousStatus !== updatedCalculation.status;
      const hasErrorMessageChanged =
        currentErrorMessage !== previousErrorMessage;
      const shouldUpdate =
        hasStatusChanged || hasErrorMessageChanged || !previousStatus;

      if (!shouldUpdate) {
        console.log(
          `[UnifiedWebSocket] Skipping update for ${calculationId}: no meaningful changes detected`
        );
        return;
      }

      console.log(
        `[UnifiedWebSocket] Status change detected for ${calculationId}: ${previousStatus || 'pending'} -> ${updatedCalculation.status}`
      );

      // キャッシュ更新
      updateCalculationCache(updatedCalculation);

      // 通知処理
      handleNotifications(updatedCalculation, previousStatus);
    },
    [
      isDuplicateUpdate,
      queryClient,
      updateCalculationCache,
      handleNotifications,
    ]
  );

  // エラーハンドリング
  const handleWebSocketError = useCallback(
    (error: string) => {
      console.error('[UnifiedWebSocket] Error:', error);
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
      showErrorNotification('Real-time monitoring error', error);
    },
    [queryClient]
  );

  // 遅延切断タイマーをクリア
  const clearDisconnectTimer = useCallback(() => {
    if (disconnectTimerRef.current) {
      clearTimeout(disconnectTimerRef.current);
      disconnectTimerRef.current = null;
    }
  }, []);

  // WebSocket切断
  const disconnect = useCallback(() => {
    clearDisconnectTimer();
    if (socketRef.current) {
      const socket = socketRef.current;
      try {
        // Leave all rooms
        socket.emit('leave_global_updates');
        if (currentActiveCalculationId.current) {
          socket.emit('leave_calculation', {
            calculation_id: currentActiveCalculationId.current,
          });
        }
        socket.disconnect();
      } catch (error) {
        console.error('[UnifiedWebSocket] Error disconnecting:', error);
      }
      socketRef.current = null;
    }
    currentActiveCalculationId.current = null;
  }, [clearDisconnectTimer]);

  // WebSocket接続
  const connect = useCallback(() => {
    // 既存の接続が有効な場合はスキップ
    if (socketRef.current?.connected) {
      console.log(
        '[UnifiedWebSocket] Already connected, skipping reconnection'
      );
      return;
    }

    // 既存の接続があれば切断
    disconnect();

    try {
      const serverUrl = `http://127.0.0.1:${window.flaskPort}`;
      console.log(`[UnifiedWebSocket] Connecting to ${serverUrl}`);

      const socket = io(serverUrl, {
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
        console.log(
          '[UnifiedWebSocket] Connected, joining global_updates room'
        );
        socket.emit('join_global_updates');
      });

      socket.on('calculation_update', handleCalculationUpdate);

      socket.on('error', (errorData: any) => {
        console.error('[UnifiedWebSocket] Socket error:', errorData);
        const errorMessage =
          errorData?.error || 'A problem occurred with calculation monitoring.';
        handleWebSocketError(errorMessage);
      });

      socket.on('connect_error', (error: Error) => {
        console.error('[UnifiedWebSocket] Connection error:', error);
        let errorMessage = 'Failed to connect to monitoring server.';

        if (error.message.includes('timeout')) {
          errorMessage = 'Connection timed out. The server may not be running.';
        } else if (error.message.includes('xhr poll error')) {
          errorMessage = 'Server communication was interrupted.';
        } else if (error.message.includes('websocket error')) {
          errorMessage =
            'WebSocket connection failed. Falling back to polling mode.';
        } else if (error.message.includes('400')) {
          errorMessage =
            'Server rejected connection (HTTP 400). Check CORS settings.';
        } else if (error.message.includes('403')) {
          errorMessage =
            'Connection forbidden (HTTP 403). Check server authentication.';
        } else if (
          error.message.includes('502') ||
          error.message.includes('503')
        ) {
          errorMessage =
            'Server is temporarily unavailable. Retrying connection...';
        }

        if (!error.message.includes('502') && !error.message.includes('503')) {
          handleWebSocketError(errorMessage);
        }
      });

      socket.on('disconnect', (reason: string) => {
        console.log(`[UnifiedWebSocket] Disconnected, reason: ${reason}`);
        if (socketRef.current === socket) {
          socketRef.current = null;
        }
      });

      socket.on('reconnect', () => {
        console.log('[UnifiedWebSocket] Reconnected successfully');
      });

      socket.on('reconnect_error', (error: Error) => {
        console.error('[UnifiedWebSocket] Reconnection failed:', error);
      });
    } catch (error) {
      console.error('[UnifiedWebSocket] Failed to create connection:', error);
      let errorMessage = 'Failed to start monitoring.';

      if (error instanceof Error) {
        if (error.message.includes('network')) {
          errorMessage = 'Unable to connect to server due to network error.';
        } else if (error.message.includes('security')) {
          errorMessage = 'Connection blocked by security settings.';
        }
      }

      handleWebSocketError(errorMessage);
    }
  }, [disconnect, handleCalculationUpdate, handleWebSocketError]);

  // アクティブ計算ルームの管理
  const manageActiveCalculationRoom = useCallback(
    (newActiveId: string | null) => {
      const socket = socketRef.current;
      if (!socket?.connected) return;

      const previousId = currentActiveCalculationId.current;

      // 前のルームから退出
      if (previousId && !previousId.startsWith('new-calculation-')) {
        console.log(
          `[UnifiedWebSocket] Leaving calculation room: ${previousId}`
        );
        socket.emit('leave_calculation', { calculation_id: previousId });
      }

      // 新しいルームに参加
      if (newActiveId && !newActiveId.startsWith('new-calculation-')) {
        console.log(
          `[UnifiedWebSocket] Joining calculation room: ${newActiveId}`
        );
        socket.emit('join_calculation', { calculation_id: newActiveId });
      }

      currentActiveCalculationId.current = newActiveId;
    },
    []
  );

  // WebSocket接続の管理
  useEffect(() => {
    connect();
    return () => disconnect();
  }, [connect, disconnect]);

  // アクティブ計算IDの変更監視
  useEffect(() => {
    manageActiveCalculationRoom(activeCalculationId);
  }, [activeCalculationId, manageActiveCalculationRoom]);

  return {
    isConnected: socketRef.current?.connected || false,
    reconnect: connect,
    disconnect,
  };
};
