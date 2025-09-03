import { useQueryClient } from '@tanstack/react-query';
import {
  useStartCalculation,
  useUpdateCalculationName,
  useDeleteCalculation,
} from './useCalculationQueries';
import {
  CalculationInstance,
  QuantumCalculationRequest,
  ApiError,
} from '../types/api-types';
import { useCalculationStore } from '../store/calculationStore';
import {
  showErrorNotification,
  showInfoNotification,
} from '../store/notificationStore';

export const useCalculationActions = () => {
  const queryClient = useQueryClient();
  const startCalculationMutation = useStartCalculation();
  const updateCalculationNameMutation = useUpdateCalculationName();
  const deleteCalculationMutation = useDeleteCalculation();
  const { clearStagedCalculation, setActiveCalculationId } =
    useCalculationStore();

  const handleApiError = (error: unknown, defaultMessage: string) => {
    console.error(defaultMessage, error);

    if (error instanceof ApiError) {
      showErrorNotification(defaultMessage, error.getUserMessage());
    } else if (error instanceof Error) {
      showErrorNotification(defaultMessage, error.message);
    } else {
      showErrorNotification(defaultMessage, '不明なエラーが発生しました。');
    }
  };

  const handleStartCalculation = async (
    calculationParams: QuantumCalculationRequest
  ): Promise<CalculationInstance> => {
    try {
      const response =
        await startCalculationMutation.mutateAsync(calculationParams);
      const runningCalculation = response.calculation;

      // 新規計算作成時の後処理
      clearStagedCalculation();
      setActiveCalculationId(runningCalculation.id);

      // ステータスに基づく通知
      if (runningCalculation.status === 'waiting') {
        const waitingReason =
          runningCalculation.waitingReason ||
          'システムリソースの空きまたは実行スロットの空きをお待ちください。';
        showInfoNotification(
          '計算を待機中です',
          `${calculationParams.name}の計算は現在待機中です。理由: ${waitingReason}`
        );
      } else if (runningCalculation.status === 'running') {
        showInfoNotification(
          '計算を開始しました',
          `${calculationParams.name}の計算を開始しました。`
        );
      } else if (runningCalculation.status === 'pending') {
        showInfoNotification(
          '計算を準備中です',
          `${calculationParams.name}の計算を準備しています。まもなく開始されます。`
        );
      } else {
        // その他のステータス（error等）の場合は汎用メッセージ
        showInfoNotification(
          '計算をリクエストしました',
          `${calculationParams.name}の計算をリクエストしました。ステータス: ${runningCalculation.status}`
        );
      }

      return runningCalculation;
    } catch (error) {
      handleApiError(error, '計算の開始に失敗しました');
      throw error;
    }
  };

  const handleCalculationRename = async (id: string, newName: string) => {
    try {
      await updateCalculationNameMutation.mutateAsync({ id, newName });
    } catch (error) {
      console.error('Failed to rename calculation:', error);
      throw error;
    }
  };

  const handleCalculationDelete = async (calculationId: string) => {
    try {
      await deleteCalculationMutation.mutateAsync(calculationId);
      // 削除された計算がアクティブだった場合の後処理は呼び出し元で処理
    } catch (error) {
      handleApiError(error, '計算の削除に失敗しました');
      throw error;
    }
  };

  const handleCalculationUpdate = (updatedCalculation: CalculationInstance) => {
    // React Queryキャッシュを直接更新
    queryClient.setQueryData(['calculation', updatedCalculation.id], {
      calculation: updatedCalculation,
    });
    queryClient.invalidateQueries({ queryKey: ['calculations'] });
  };

  return {
    handleStartCalculation,
    handleCalculationRename,
    handleCalculationDelete,
    handleCalculationUpdate,

    // Loading states
    isStarting: startCalculationMutation.isPending,
    isRenaming: updateCalculationNameMutation.isPending,
    isDeleting: deleteCalculationMutation.isPending,
  };
};
