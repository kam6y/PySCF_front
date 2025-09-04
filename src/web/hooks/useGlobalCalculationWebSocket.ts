import { useEffect, useRef, useCallback } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import { io, Socket } from 'socket.io-client';
import { CalculationInstance } from '../types/api-types';
import {
  showErrorNotification,
  showSuccessNotification,
  showInfoNotification,
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

  // 計算更新の統一処理ロジック
  const handleCalculationUpdate = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;
      const isActiveCalculation = calculationId === activeCalculationId;

      // 前のステータスを取得してステータス変化を検出
      const previousData = queryClient.getQueryData([
        'calculation',
        calculationId,
      ]) as { calculation: CalculationInstance } | undefined;
      const previousStatus = previousData?.calculation?.status;

      // 1. 個別計算詳細のキャッシュを更新
      queryClient.setQueryData(['calculation', calculationId], {
        calculation: updatedCalculation,
      });

      // 2. 計算リストのキャッシュを更新（即座に反映）
      queryClient.setQueryData(['calculations'], (oldData: any) => {
        if (!oldData?.calculations) return oldData;

        const updatedCalculations = oldData.calculations.map((calc: any) =>
          calc.id === calculationId
            ? {
                ...calc,
                status: updatedCalculation.status,
                date: updatedCalculation.updatedAt,
                name: updatedCalculation.name,
              }
            : calc
        );

        return {
          ...oldData,
          calculations: updatedCalculations,
        };
      });

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
        if (!isActiveCalculation) {
          const errorMessage =
            updatedCalculation.errorMessage ||
            updatedCalculation.results?.error ||
            'Detailed error information is not available.';

          showErrorNotification(
            `Calculation "${molecularName}" failed`,
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

      // 4. バックグラウンドでの最新データ取得用（念のため）
      queryClient.invalidateQueries({
        queryKey: ['calculations'],
      });
    },
    [queryClient, activeCalculationId]
  );

  // エラーハンドリング
  const handleWebSocketError = useCallback((error: string) => {
    console.error('Global WebSocket error:', error);
    showErrorNotification('Real-time monitoring error', error);
  }, []);

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
    // 既存の接続があれば切断
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

      socket.on('connect', () => {
        console.log('Global Socket.IO connected, joining global_updates room');
        // Join the global updates room to receive all calculation updates
        socket.emit('join_global_updates');
      });

      socket.on(
        'calculation_update',
        (updatedCalculation: CalculationInstance) => {
          try {
            console.log(
              'Global WebSocket received update for:',
              updatedCalculation.id,
              'status:',
              updatedCalculation.status
            );
            handleCalculationUpdate(updatedCalculation);
          } catch (e) {
            console.error('Failed to process global calculation update:', e);
            handleWebSocketError(
              'Unable to parse response from server. Calculation status updates may be temporarily stopped.'
            );
          }
        }
      );

      socket.on('error', (errorData: any) => {
        console.error('Global Socket.IO error:', errorData);
        const errorMessage =
          errorData?.error || 'A problem occurred with global calculation monitoring.';
        handleWebSocketError(errorMessage);
      });

      socket.on('connect_error', (error: Error) => {
        console.error('Global Socket.IO connection error:', error);
        let errorMessage = 'Failed to connect to global monitoring server.';

        if (error.message.includes('timeout')) {
          errorMessage =
            'Global monitoring connection timed out. The server may not be running.';
        } else if (error.message.includes('xhr poll error')) {
          errorMessage = 'Global monitoring server communication was interrupted.';
        }

        handleWebSocketError(errorMessage);
      });

      socket.on('disconnect', (reason: string) => {
        console.log('Global Socket.IO disconnected, reason:', reason);
        if (socketRef.current === socket) {
          socketRef.current = null;
        }

        // Auto-reconnect for certain disconnect reasons
        if (reason === 'io server disconnect') {
          // Server initiated disconnect, likely planned shutdown
          handleWebSocketError(
            'Connection to global monitoring server was disconnected. Please refresh the page to check status.'
          );
        }
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
