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
  showResourceInsufficientErrorNotification,
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
      showErrorNotification(defaultMessage, 'An unknown error occurred.');
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
          'Please wait for available system resources or execution slots.';
        showInfoNotification(
          'Calculation is waiting',
          `Calculation for ${calculationParams.name} is currently waiting. Reason: ${waitingReason}`
        );
      } else if (runningCalculation.status === 'running') {
        showInfoNotification(
          'Calculation started',
          `Calculation for ${calculationParams.name} has been started.`
        );
      } else if (runningCalculation.status === 'pending') {
        showInfoNotification(
          'Preparing calculation',
          `Preparing calculation for ${calculationParams.name}. It will start soon.`
        );
      } else if (runningCalculation.status === 'error') {
        // エラーステータスの場合はエラーメッセージを確認してリソース不足エラーかを判定
        const errorMessage = runningCalculation.errorMessage || runningCalculation.results?.error;
        
        if (errorMessage) {
          // リソース不足エラーの判定
          const isResourceInsufficientError = 
            errorMessage.toLowerCase().includes('cpu usage') ||
            errorMessage.toLowerCase().includes('memory usage') ||
            errorMessage.toLowerCase().includes('system cpu usage') ||
            errorMessage.toLowerCase().includes('system memory usage') ||
            errorMessage.toLowerCase().includes('no active calculations');

          if (isResourceInsufficientError) {
            showResourceInsufficientErrorNotification(
              errorMessage,
              runningCalculation.id
            );
          } else {
            showErrorNotification(
              `Calculation "${calculationParams.name}" failed`,
              errorMessage,
              runningCalculation.id
            );
          }
        } else {
          showErrorNotification(
            `Calculation "${calculationParams.name}" failed`,
            'Detailed error information is not available.',
            runningCalculation.id
          );
        }
      } else {
        // その他のステータスの場合は汎用メッセージ
        showInfoNotification(
          'Calculation requested',
          `Calculation for ${calculationParams.name} has been requested. Status: ${runningCalculation.status}`
        );
      }

      return runningCalculation;
    } catch (error) {
      handleApiError(error, 'Failed to start calculation');
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
      handleApiError(error, 'Failed to delete calculation');
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
