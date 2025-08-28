import { useQueryClient } from '@tanstack/react-query';
import { useCalculationSubscription } from './useCalculationSubscription';
import { CalculationInstance } from '../types/api-types';
import { showErrorNotification, showSuccessNotification } from '../store/notificationStore';

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

  // 統一されたキャッシュ更新ロジック
  const handleCalculationUpdate = (updatedCalculation: CalculationInstance) => {
    // 前のステータスを取得してステータス変化を検出
    const previousData = queryClient.getQueryData(['calculation', updatedCalculation.id]) as { calculation: CalculationInstance } | undefined;
    const previousStatus = previousData?.calculation?.status;

    // 計算完了通知: running -> completed の変化を検出
    if (previousStatus === 'running' && updatedCalculation.status === 'completed') {
      showSuccessNotification('計算が完了しました', `${updatedCalculation.name}の計算が完了しました。`);
    }

    // 1. 個別計算詳細のキャッシュを更新
    queryClient.setQueryData(['calculation', updatedCalculation.id], {
      calculation: updatedCalculation,
    });

    // 2. 計算リストのキャッシュを更新（即座に反映）
    queryClient.setQueryData(['calculations'], (oldData: any) => {
      if (!oldData?.calculations) return oldData;

      const updatedCalculations = oldData.calculations.map((calc: any) =>
        calc.id === updatedCalculation.id
          ? {
              ...calc,
              status: updatedCalculation.status,
              date: updatedCalculation.updatedAt,
            }
          : calc
      );

      return {
        ...oldData,
        calculations: updatedCalculations,
      };
    });

    // 3. バックグラウンドでの最新データ取得用（念のため）
    queryClient.invalidateQueries({
      queryKey: ['calculations'],
    });
  };

  // エラーハンドリング
  const handleWebSocketError = (error: string) => {
    console.error('WebSocket error:', error);

    // WebSocketエラーをトースト通知で表示
    showErrorNotification('リアルタイム監視エラー', error);
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
