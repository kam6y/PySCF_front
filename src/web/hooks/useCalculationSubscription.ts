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

  // コールバック関数の参照を最新に保つ
  useEffect(() => {
    onUpdateRef.current = onUpdate;
    onErrorRef.current = onError;
  }, [onUpdate, onError]);

  const disconnect = useCallback(() => {
    if (wsRef.current) {
      const ws = wsRef.current;
      if (
        ws.readyState === WebSocket.OPEN ||
        ws.readyState === WebSocket.CONNECTING
      ) {
        ws.close();
      }
      wsRef.current = null;
    }
    currentCalculationIdRef.current = null;
  }, []);

  const connect = useCallback(
    (calcId: string) => {
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
          console.log(`WebSocket connected for calculation ${calcId}`);
        };

        ws.onmessage = event => {
          try {
            const updatedCalculation = JSON.parse(
              event.data
            ) as CalculationInstance;
            if (onUpdateRef.current) {
              onUpdateRef.current(updatedCalculation);
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
    // 計算が実行中または保留中でかつIDが有効な場合のみ接続
    // completed/errorになるまで接続を維持することで、最終状態の更新を確実に受信
    if (calculationId && (status === 'running' || status === 'pending')) {
      connect(calculationId);
    } else {
      disconnect();
    }

    // クリーンアップ
    return () => {
      disconnect();
    };
  }, [calculationId, status, connect, disconnect]);
};
