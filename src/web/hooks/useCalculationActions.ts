import { useQueryClient } from '@tanstack/react-query';
import {
  calculationQueryKeys,
  useStartCalculation,
  useUpdateCalculationName,
  useDeleteCalculation,
  usePauseCalculation,
  useResumeCalculation,
} from './useCalculationQueries';
import {
  CalculationInstance,
  QuantumCalculationRequest,
} from '../types/api-types';
import { useCalculationStore } from '../store/calculationStore';
import { showInfoNotification } from '../store/notificationStore';
import { handleError } from '../utils/errorHandler';
import { isResourceInsufficientError } from '../utils/errorClassifier';

export const useCalculationActions = () => {
  const queryClient = useQueryClient();
  const startCalculationMutation = useStartCalculation();
  const updateCalculationNameMutation = useUpdateCalculationName();
  const deleteCalculationMutation = useDeleteCalculation();
  const pauseCalculationMutation = usePauseCalculation();
  const resumeCalculationMutation = useResumeCalculation();
  const { clearStagedCalculation, setActiveCalculationId } =
    useCalculationStore();

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
        const errorMessage =
          runningCalculation.error || runningCalculation.results?.error;

        if (errorMessage) {
          const isResourceError = isResourceInsufficientError(errorMessage);

          if (isResourceError) {
            handleError(new Error(errorMessage), 'Calculation failed');
          } else {
            handleError(
              new Error(errorMessage),
              `Calculation "${calculationParams.name}" failed`
            );
          }
        } else {
          handleError(
            new Error('Detailed error information is not available.'),
            `Calculation "${calculationParams.name}" failed`
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
      handleError(error, 'Failed to start calculation');
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
      handleError(error, 'Failed to delete calculation');
      throw error;
    }
  };

  const handleCalculationPause = async (calculationId: string) => {
    try {
      await pauseCalculationMutation.mutateAsync(calculationId);
      showInfoNotification(
        'Calculation pausing',
        'The calculation will pause after the current iteration completes.'
      );
    } catch (error) {
      handleError(error, 'Failed to pause calculation');
      throw error;
    }
  };

  const handleCalculationResume = async (calculationId: string) => {
    try {
      await resumeCalculationMutation.mutateAsync(calculationId);
      showInfoNotification(
        'Calculation resuming',
        'The calculation is resuming from where it was paused.'
      );
    } catch (error) {
      handleError(error, 'Failed to resume calculation');
      throw error;
    }
  };

  const handleCalculationUpdate = (updatedCalculation: CalculationInstance) => {
    // React Queryキャッシュを直接更新
    queryClient.setQueryData(
      calculationQueryKeys.detail(updatedCalculation.id),
      {
        calculation: updatedCalculation,
      }
    );
    queryClient.invalidateQueries({
      queryKey: calculationQueryKeys.list(),
    });
  };

  return {
    handleStartCalculation,
    handleCalculationRename,
    handleCalculationDelete,
    handleCalculationPause,
    handleCalculationResume,
    handleCalculationUpdate,

    // Loading states
    isStarting: startCalculationMutation.isPending,
    isRenaming: updateCalculationNameMutation.isPending,
    isDeleting: deleteCalculationMutation.isPending,
    isPausing: pauseCalculationMutation.isPending,
    isResuming: resumeCalculationMutation.isPending,
  };
};
