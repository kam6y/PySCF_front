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
              onErrorRef.current('Failed to parse server response.');
            }
          }
        };

        ws.onerror = event => {
          console.error('WebSocket error:', event);
          if (onErrorRef.current) {
            onErrorRef.current(
              'Connection to the calculation status server failed.'
            );
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
          onErrorRef.current('Failed to establish WebSocket connection.');
        }
      }
    },
    [disconnect]
  );

  useEffect(() => {
    // 計算が実行中でかつIDが有効な場合のみ接続
    if (calculationId && status === 'running') {
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
