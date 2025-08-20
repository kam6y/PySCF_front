// src/web/hooks/useCalculationSubscription.ts

import { useEffect, useRef, useCallback } from 'react';
import { CalculationInstance } from '../types/api-types';
import { getWebSocketUrl } from '../apiClient';

export interface UseCalculationSubscriptionOptions {
  calculationId: string | null;
  status?: 'pending' | 'running' | 'completed' | 'error';
  onUpdate: (calculation: CalculationInstance) => void;
  onError?: (error: string) => void;
}

export const useCalculationSubscription = ({
  calculationId,
  status,
  onUpdate,
  onError,
}: UseCalculationSubscriptionOptions) => {
  const wsRef = useRef<WebSocket | null>(null);
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
    if (wsRef.current) {
      const ws = wsRef.current;
      try {
        // WebSocketの状態に関係なく、強制的に接続を閉じる
        if (ws.readyState !== WebSocket.CLOSED) {
          ws.close();
        }
      } catch (error) {
        console.error('Error closing WebSocket:', error);
      }
      wsRef.current = null;
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
        wsRef.current &&
        (wsRef.current.readyState === WebSocket.CONNECTING ||
          wsRef.current.readyState === WebSocket.OPEN)
      ) {
        return;
      }

      // 既存接続をクリーンアップ
      disconnect();

      try {
        const wsUrl = getWebSocketUrl(`/ws/calculations/${calcId}`);
        const ws = new WebSocket(wsUrl);
        wsRef.current = ws;
        currentCalculationIdRef.current = calcId;

        ws.onopen = () => {
          // WebSocket接続確立
        };

        ws.onmessage = event => {
          try {
            const updatedCalculation = JSON.parse(
              event.data
            ) as CalculationInstance;
            
            if (onUpdateRef.current) {
              onUpdateRef.current(updatedCalculation);
            }
            
            // 計算が完了または失敗した場合、少し待ってから接続を切断
            if (updatedCalculation.status === 'completed' || updatedCalculation.status === 'error') {
              // 既存のタイマーをクリア
              clearDisconnectTimer();
              
              // 3秒後に切断
              disconnectTimerRef.current = setTimeout(() => {
                disconnect();
              }, 3000);
            }
          } catch (e) {
            console.error('Failed to parse WebSocket message:', e);
            if (onErrorRef.current) {
              onErrorRef.current(
                'サーバーからの応答が解析できませんでした。計算状況の更新が一時的に停止する可能性があります。'
              );
            }
          }
        };

        ws.onerror = event => {
          console.error('WebSocket error:', event);
          if (onErrorRef.current) {
            // WebSocketの状態に基づいてより具体的なエラーメッセージを提供
            let errorMessage = '計算状況の監視に問題が発生しました。';

            if (ws.readyState === WebSocket.CONNECTING) {
              errorMessage =
                'サーバーへの接続に失敗しました。ネットワーク接続を確認してください。';
            } else if (ws.readyState === WebSocket.OPEN) {
              errorMessage =
                'サーバーとの通信が中断されました。計算は継続中ですが、状況の更新が停止しています。';
            } else if (ws.readyState === WebSocket.CLOSED) {
              errorMessage =
                'サーバーとの接続が切断されました。ページを更新して状況を確認してください。';
            }

            onErrorRef.current(errorMessage);
          }
        };

        ws.onclose = event => {
          console.log(
            `WebSocket disconnected for calculation ${calcId}`,
            event.code,
            event.reason
          );
          if (wsRef.current === ws) {
            wsRef.current = null;
          }
          if (currentCalculationIdRef.current === calcId) {
            currentCalculationIdRef.current = null;
          }
        };
      } catch (error) {
        console.error('Failed to create WebSocket connection:', error);
        if (onErrorRef.current) {
          let errorMessage = 'リアルタイム監視の開始に失敗しました。';

          if (error instanceof Error) {
            // 特定のエラータイプに基づいてメッセージを調整
            if (error.name === 'SecurityError') {
              errorMessage =
                'セキュリティ設定により接続できませんでした。HTTPSが必要な可能性があります。';
            } else if (error.name === 'SyntaxError') {
              errorMessage =
                'WebSocketのURL形式が無効です。設定を確認してください。';
            } else if (error.message.includes('network')) {
              errorMessage =
                'ネットワークエラーによりサーバーに接続できませんでした。';
            }
          }

          onErrorRef.current(errorMessage);
        }
      }
    },
    [disconnect]
  );

  useEffect(() => {
    // 計算IDが有効で一時的IDでない場合に接続
    // running/pending状態では即座に接続
    // completed/error状態でも接続を維持（自動切断タイマーで後で切断）
    if (
      calculationId &&
      !calculationId.startsWith('new-calculation-')
    ) {
      // running/pending状態では即座に接続
      if (status === 'running' || status === 'pending') {
        connect(calculationId);
      } else if (status === 'completed' || status === 'error') {
        // completed/error状態でも、まだ接続していない場合は接続する
        // （既に接続している場合は、onmessageハンドラで自動切断タイマーが設定される）
        if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
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
