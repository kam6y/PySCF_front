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

/**
 * グローバルWebSocket監視フック
 *
 * 全ての計算の状態変化をリアルタイムで監視し、以下の機能を提供：
 * - 非アクティブな計算の完了・エラー通知
 * - サイドバー計算リストのリアルタイム更新
 * - React Queryキャッシュの統一的な更新
 */
export const useGlobalCalculationWebSocket = (
  activeCalculationId: string | null
) => {
  const queryClient = useQueryClient();
  const socketRef = useRef<Socket | null>(null);
  const disconnectTimerRef = useRef<NodeJS.Timeout | null>(null);
  const lastUpdateTimestamps = useRef<Map<string, number>>(new Map());

  // 計算更新の統一処理ロジック
  const handleCalculationUpdate = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;
      const isActiveCalculation = calculationId === activeCalculationId;
      const currentTimestamp = Date.now();

      // 重複更新防止：同じ計算の更新が短時間内に複数回来た場合はスキップ
      const lastUpdateTime = lastUpdateTimestamps.current.get(calculationId);
      if (lastUpdateTime && currentTimestamp - lastUpdateTime < 100) {
        // 100ms以内の重複更新は無視
        console.log(
          `[WebSocket] Skipping duplicate update for calculation ${calculationId} (within 100ms)`
        );
        return;
      }
      lastUpdateTimestamps.current.set(calculationId, currentTimestamp);

      console.log(
        `[WebSocket] Processing update for calculation ${calculationId}: ${updatedCalculation.status}`
      );

      // 前のステータスを取得してステータス変化を検出
      const previousDetailData = queryClient.getQueryData([
        'calculation',
        calculationId,
      ]) as { calculation: CalculationInstance } | undefined;
      const previousStatus = previousDetailData?.calculation?.status;

      // 計算リストから前のステータスも確認
      const calculationsData = queryClient.getQueryData(['calculations']) as
        | { calculations: any[] }
        | undefined;
      const previousListItem = calculationsData?.calculations?.find(
        calc => calc.id === calculationId
      );
      const previousListStatus = previousListItem?.status;

      // ステータス・エラーメッセージ変更の検証：ステータス変化またはエラーメッセージの変化があった場合に更新
      const currentErrorMessage =
        updatedCalculation.errorMessage || updatedCalculation.results?.error;
      const previousErrorMessage =
        previousDetailData?.calculation?.errorMessage ||
        previousDetailData?.calculation?.results?.error;

      const hasStatusChanged =
        previousStatus !== updatedCalculation.status ||
        previousListStatus !== updatedCalculation.status;

      const hasErrorMessageChanged =
        currentErrorMessage !== previousErrorMessage;

      const shouldUpdate =
        hasStatusChanged || hasErrorMessageChanged || !previousStatus;

      if (!shouldUpdate) {
        console.log(
          `[WebSocket] Skipping update for ${calculationId}: no meaningful changes detected`
        );
        return;
      }

      console.log(
        `[WebSocket] Status change detected for ${calculationId}: ${previousStatus || 'pending'} -> ${updatedCalculation.status}`
      );

      try {
        // 1. 個別計算詳細のキャッシュを更新
        queryClient.setQueryData(['calculation', calculationId], {
          calculation: updatedCalculation,
        });

        // 2. 計算リストのキャッシュを更新（楽観的更新）
        queryClient.setQueryData(['calculations'], (oldData: any) => {
          if (!oldData?.calculations) {
            console.warn(
              '[WebSocket] No calculations data found for cache update'
            );
            return oldData;
          }

          const calculationIndex = oldData.calculations.findIndex(
            (calc: any) => calc.id === calculationId
          );

          let updatedCalculations;
          if (calculationIndex === -1) {
            // 計算がリストに存在しない場合は新規追加（楽観的更新）
            const newCalculationItem = {
              id: calculationId,
              name: updatedCalculation.name || calculationId,
              status: updatedCalculation.status,
              date:
                updatedCalculation.updatedAt ||
                updatedCalculation.createdAt ||
                new Date().toISOString(),
              // 最小限必要な情報のみを設定
            };

            // リストの先頭に新しい計算を追加
            updatedCalculations = [newCalculationItem, ...oldData.calculations];
            console.log(
              `[WebSocket] Added new calculation ${calculationId} to cache with status ${updatedCalculation.status}`
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
              `[WebSocket] Updated existing calculation ${calculationId} with status ${updatedCalculation.status}`
            );
          }

          return {
            ...oldData,
            calculations: updatedCalculations,
          };
        });
      } catch (error) {
        console.error(
          `[WebSocket] Failed to update cache for calculation ${calculationId}:`,
          error
        );
        // フォールバック：エラーが発生した場合は関連クエリを無効化
        queryClient.invalidateQueries({
          queryKey: ['calculation', calculationId],
        });
        queryClient.invalidateQueries({ queryKey: ['calculations'] });
        return;
      }

      // 3. 状態変化に基づく通知処理
      const molecularName = updatedCalculation.name || 'Unknown';

      // 計算完了通知: running -> completed の変化を検出
      if (
        previousStatus === 'running' &&
        updatedCalculation.status === 'completed'
      ) {
        // アクティブな計算の場合は既に個別のWebSocketで通知済みの可能性があるため、
        // 非アクティブな計算のみ通知
        if (!isActiveCalculation) {
          showSuccessNotification(
            'Calculation completed',
            `Calculation for ${molecularName} has been completed.`,
            calculationId
          );
        }
      }

      // 計算エラー通知: running -> error の変化を検出
      if (
        previousStatus === 'running' &&
        updatedCalculation.status === 'error'
      ) {
        const errorMessage =
          updatedCalculation.errorMessage ||
          updatedCalculation.results?.error ||
          'Detailed error information is not available.';

        // リソース不足エラーの判定
        const isResourceInsufficientError =
          errorMessage &&
          (errorMessage.toLowerCase().includes('cpu usage') ||
            errorMessage.toLowerCase().includes('memory usage') ||
            errorMessage.toLowerCase().includes('system cpu usage') ||
            errorMessage.toLowerCase().includes('system memory usage') ||
            errorMessage.toLowerCase().includes('no active calculations'));

        if (isResourceInsufficientError) {
          showResourceInsufficientErrorNotification(
            errorMessage,
            calculationId
          );
        } else {
          showErrorNotification(
            `Calculation "${molecularName}" failed`,
            errorMessage,
            calculationId
          );
        }
      }

      // 計算エラー通知: pending -> error の変化を検出（開始時のエラー）
      if (
        previousStatus === 'pending' &&
        updatedCalculation.status === 'error'
      ) {
        const errorMessage =
          updatedCalculation.errorMessage ||
          updatedCalculation.results?.error ||
          'Detailed error information is not available.';

        // リソース不足エラーの判定
        const isResourceInsufficientError =
          errorMessage &&
          (errorMessage.toLowerCase().includes('cpu usage') ||
            errorMessage.toLowerCase().includes('memory usage') ||
            errorMessage.toLowerCase().includes('system cpu usage') ||
            errorMessage.toLowerCase().includes('system memory usage') ||
            errorMessage.toLowerCase().includes('no active calculations'));

        if (isResourceInsufficientError) {
          showResourceInsufficientErrorNotification(
            errorMessage,
            calculationId
          );
        } else {
          showErrorNotification(
            `Calculation "${molecularName}" failed to start`,
            errorMessage,
            calculationId
          );
        }
      }

      // 計算開始通知: pending -> running の変化を検出
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

      // 待機からの開始通知: waiting -> running の変化を検出
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

      // 待機状態への移行通知: pending -> waiting の変化を検出
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

      // レースコンディション回避のため、手動キャッシュ更新後の無効化は削除
      // WebSocketでリアルタイム更新されるため、追加のサーバー取得は不要
    },
    [queryClient, activeCalculationId]
  );

  // エラーハンドリング（改善版）
  const handleWebSocketError = useCallback(
    (error: string) => {
      console.error('[WebSocket] Global monitoring error:', error);
      // WebSocketエラーの場合は、キャッシュ無効化でフォールバック
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
      showErrorNotification('Real-time monitoring error', error);
    },
    [queryClient]
  );

  // 遅延切断タイマーをクリアするヘルパー関数
  const clearDisconnectTimer = useCallback(() => {
    if (disconnectTimerRef.current) {
      clearTimeout(disconnectTimerRef.current);
      disconnectTimerRef.current = null;
    }
  }, []);

  // WebSocket接続の切断
  const disconnect = useCallback(() => {
    clearDisconnectTimer();
    if (socketRef.current) {
      const socket = socketRef.current;
      try {
        // Leave global updates room
        socket.emit('leave_global_updates');

        // Disconnect from Socket.IO server
        socket.disconnect();
      } catch (error) {
        console.error('Error disconnecting from global WebSocket:', error);
      }
      socketRef.current = null;
    }
  }, [clearDisconnectTimer]);

  // WebSocket接続の確立
  const connect = useCallback(() => {
    // 既存の接続が有効な場合はスキップ
    if (socketRef.current?.connected) {
      console.log(
        '[WebSocket] Already connected to global monitoring, skipping reconnection'
      );
      return;
    }

    // 既存の接続があれば切断
    disconnect();

    try {
      // Get the Flask-SocketIO server URL (always HTTP in Electron)
      const serverUrl = `http://127.0.0.1:${window.flaskPort}`;
      console.log(
        `[WebSocket] Connecting to global monitoring at ${serverUrl}`
      );

      // Create Socket.IO connection with improved settings
      const socket = io(serverUrl, {
        transports: ['websocket', 'polling'],
        timeout: 10000, // 接続タイムアウトを延長
        reconnection: true, // 自動再接続を有効化
        reconnectionDelay: 2000,
        reconnectionAttempts: 3,
        forceNew: true, // 新しい接続を強制
      });

      socketRef.current = socket;

      socket.on('connect', () => {
        console.log(
          '[WebSocket] Global Socket.IO connected, joining global_updates room'
        );
        // Join the global updates room to receive all calculation updates
        socket.emit('join_global_updates');
      });

      socket.on(
        'calculation_update',
        (updatedCalculation: CalculationInstance) => {
          try {
            if (!updatedCalculation || !updatedCalculation.id) {
              console.warn(
                '[WebSocket] Received invalid calculation update data:',
                updatedCalculation
              );
              return;
            }

            console.log(
              `[WebSocket] Global received update for: ${updatedCalculation.id}, status: ${updatedCalculation.status}`
            );
            handleCalculationUpdate(updatedCalculation);
          } catch (e) {
            console.error(
              '[WebSocket] Failed to process global calculation update:',
              e
            );
            handleWebSocketError(
              'Unable to parse response from server. Calculation status updates may be temporarily stopped.'
            );
          }
        }
      );

      socket.on('error', (errorData: any) => {
        console.error('[WebSocket] Global Socket.IO error:', errorData);
        const errorMessage =
          errorData?.error ||
          'A problem occurred with global calculation monitoring.';
        handleWebSocketError(errorMessage);
      });

      socket.on('connect_error', (error: Error) => {
        console.error('[WebSocket] Global Socket.IO connection error:', error);
        let errorMessage = 'Failed to connect to global monitoring server.';

        if (error.message.includes('timeout')) {
          errorMessage =
            'Global monitoring connection timed out. The server may not be running.';
        } else if (error.message.includes('xhr poll error')) {
          errorMessage =
            'Global monitoring server communication was interrupted.';
        }

        handleWebSocketError(errorMessage);
      });

      socket.on('disconnect', (reason: string) => {
        console.log(
          `[WebSocket] Global Socket.IO disconnected, reason: ${reason}`
        );
        if (socketRef.current === socket) {
          socketRef.current = null;
        }

        // Auto-reconnect for certain disconnect reasons
        if (reason === 'io server disconnect') {
          // Server initiated disconnect, likely planned shutdown
          console.log(
            '[WebSocket] Server initiated disconnect, scheduling reconnection attempt'
          );
        } else if (
          reason === 'transport close' ||
          reason === 'transport error'
        ) {
          // Connection issues, let Socket.IO handle automatic reconnection
          console.log(
            '[WebSocket] Connection lost, automatic reconnection will be attempted'
          );
        }
      });

      socket.on('reconnect', () => {
        console.log('[WebSocket] Global monitoring reconnected successfully');
        // Note: 再接続時のキャッシュ無効化を削除
        // WebSocketがリアルタイムデータを提供するため、自動的な
        // サーバーフェッチは競合状態を引き起こす可能性がある
        // 必要に応じて、WebSocket更新が失敗した場合のみフォールバックとして使用
      });

      socket.on('reconnect_error', (error: Error) => {
        console.error('[WebSocket] Reconnection failed:', error);
      });
    } catch (error) {
      console.error('Failed to create global Socket.IO connection:', error);
      let errorMessage = 'Failed to start global monitoring.';

      if (error instanceof Error) {
        if (error.message.includes('network')) {
          errorMessage =
            'Unable to connect to global monitoring server due to network error.';
        } else if (error.message.includes('security')) {
          errorMessage =
            'Global monitoring connection blocked by security settings.';
        }
      }

      handleWebSocketError(errorMessage);
    }
  }, [disconnect, handleCalculationUpdate, handleWebSocketError]);

  // エフェクト: WebSocket接続の管理
  useEffect(() => {
    // グローバル監視を開始
    connect();

    // クリーンアップ
    return () => {
      disconnect();
    };
  }, [connect, disconnect]);

  // 外部インターフェース
  return {
    isConnected: socketRef.current?.connected || false,
    reconnect: connect,
    disconnect,
  };
};
