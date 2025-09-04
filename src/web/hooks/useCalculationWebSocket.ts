import { useRef } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import { useCalculationSubscription } from './useCalculationSubscription';
import { CalculationInstance } from '../types/api-types';
import {
  showErrorNotification,
  showSuccessNotification,
} from '../store/notificationStore';

/**
 * WebSocket経由の計算更新を処理する専用フック
 *
 * React Queryキャッシュの更新ロジックを統一し、
 * App.tsxから複雑なキャッシュ操作を分離
 */
export const useCalculationWebSocket = (
  calculationId: string | null,
  status?: 'pending' | 'running' | 'completed' | 'error' | 'waiting'
) => {
  const queryClient = useQueryClient();
  const lastUpdateTimestampRef = useRef<number>(0);

  // 統一されたキャッシュ更新ロジック
  const handleCalculationUpdate = (updatedCalculation: CalculationInstance) => {
    const currentTimestamp = Date.now();
    
    // 重複更新防止：グローバルWebSocketとの競合を回避
    if (lastUpdateTimestampRef.current && (currentTimestamp - lastUpdateTimestampRef.current) < 100) {
      console.log(
        `[Individual WebSocket] Skipping duplicate update for calculation ${updatedCalculation.id} (within 100ms)`
      );
      return;
    }
    lastUpdateTimestampRef.current = currentTimestamp;

    // 前のステータスを取得してステータス変化を検出
    const previousData = queryClient.getQueryData([
      'calculation',
      updatedCalculation.id,
    ]) as { calculation: CalculationInstance } | undefined;
    const previousStatus = previousData?.calculation?.status;

    // 計算完了通知: running -> completed の変化を検出
    if (
      previousStatus === 'running' &&
      updatedCalculation.status === 'completed'
    ) {
      const molecularName = updatedCalculation.name || updatedCalculation.id;
      showSuccessNotification(
        'Calculation completed',
        `Calculation for ${molecularName} has been completed.`
      );
    }

    // 1. 個別計算詳細のキャッシュを更新
    queryClient.setQueryData(['calculation', updatedCalculation.id], {
      calculation: updatedCalculation,
    });

    // Note: 計算リスト['calculations']の更新はuseGlobalCalculationWebSocketに委任
    // 責任分離により、グローバルWebSocketが全ての計算リストを統一管理し、
    // 個別WebSocketは詳細データのみに集中することで競合状態を回避
  };

  // エラーハンドリング
  const handleWebSocketError = (error: string) => {
    console.error('WebSocket error:', error);

    // WebSocketエラーをトースト通知で表示
    showErrorNotification('Real-time monitoring error', error);

    // フォールバック: WebSocketエラー時のみ、個別計算の詳細データを再取得
    // グローバルWebSocketが利用可能な場合は計算リストは更新されるため、
    // ここでは個別計算の詳細のみをフォールバックとして取得
    if (calculationId && !calculationId.startsWith('new-calculation-')) {
      queryClient.invalidateQueries({
        queryKey: ['calculation', calculationId],
      });
    }
  };

  // WebSocket接続の管理
  useCalculationSubscription({
    calculationId,
    status,
    onUpdate: handleCalculationUpdate,
    onError: handleWebSocketError,
  });

  // 外部インターフェース（現在は内部で処理完結しているが、将来拡張可能）
  return {
    // 必要に応じて、接続状態やエラー状態を返すことも可能
    isConnected: calculationId && !calculationId.startsWith('new-calculation-'),
  };
};
