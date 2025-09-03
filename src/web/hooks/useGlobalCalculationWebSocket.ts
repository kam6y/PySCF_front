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
            '計算が完了しました',
            `${molecularName}の計算が完了しました。`,
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
            '詳細なエラー情報が利用できません。';

          showErrorNotification(
            `計算「${molecularName}」が失敗しました`,
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
            '計算を開始しました',
            `${molecularName}の計算を開始しました。`,
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
    showErrorNotification('リアルタイム監視エラー', error);
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
              'サーバーからの応答が解析できませんでした。計算状況の更新が一時的に停止する可能性があります。'
            );
          }
        }
      );

      socket.on('error', (errorData: any) => {
        console.error('Global Socket.IO error:', errorData);
        const errorMessage =
          errorData?.error || 'グローバル計算監視に問題が発生しました。';
        handleWebSocketError(errorMessage);
      });

      socket.on('connect_error', (error: Error) => {
        console.error('Global Socket.IO connection error:', error);
        let errorMessage = 'グローバル監視用サーバーへの接続に失敗しました。';

        if (error.message.includes('timeout')) {
          errorMessage =
            'グローバル監視接続がタイムアウトしました。サーバーが起動していない可能性があります。';
        } else if (error.message.includes('xhr poll error')) {
          errorMessage = 'グローバル監視のサーバー通信が中断されました。';
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
            'グローバル監視サーバーとの接続が切断されました。ページを更新して状況を確認してください。'
          );
        }
      });
    } catch (error) {
      console.error('Failed to create global Socket.IO connection:', error);
      let errorMessage = 'グローバル監視の開始に失敗しました。';

      if (error instanceof Error) {
        if (error.message.includes('network')) {
          errorMessage =
            'ネットワークエラーによりグローバル監視サーバーに接続できませんでした。';
        } else if (error.message.includes('security')) {
          errorMessage =
            'セキュリティ設定によりグローバル監視接続できませんでした。';
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
