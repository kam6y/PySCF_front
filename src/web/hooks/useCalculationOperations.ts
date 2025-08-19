import {
  useStartCalculation,
  useUpdateCalculationName,
  useDeleteCalculation,
} from './useCalculationQueries';
import { useCalculationStore } from '../store/calculationStore';
import {
  QuantumCalculationRequest,
  CalculationInstance,
  ApiError,
} from '../types/api-types';

export interface CalculationOperations {
  handleStartCalculation: (
    params: QuantumCalculationRequest
  ) => Promise<CalculationInstance>;
  handleCalculationRename: (
    calculationId: string,
    newName: string
  ) => Promise<void>;
  handleCalculationDelete: (calculationId: string) => Promise<void>;
}

export const useCalculationOperations = (
  setCurrentPage: (page: 'calculation-settings') => void
): CalculationOperations => {
  const startCalculationMutation = useStartCalculation();
  const updateCalculationNameMutation = useUpdateCalculationName();
  const deleteCalculationMutation = useDeleteCalculation();

  const {
    activeCalculationId,
    setActiveCalculationId,
    clearStagedCalculation,
  } = useCalculationStore();

  const handleApiError = (error: unknown, defaultMessage: string) => {
    console.error(defaultMessage, error);

    if (error instanceof ApiError) {
      alert(error.getUserMessage());
    } else if (error instanceof Error) {
      alert(`${defaultMessage}: ${error.message}`);
    } else {
      alert(`${defaultMessage}: 不明なエラーが発生しました。`);
    }
  };

  const handleStartCalculation = async (
    params: QuantumCalculationRequest
  ): Promise<CalculationInstance> => {
    try {
      const response = await startCalculationMutation.mutateAsync(params);
      const runningCalculation = response.calculation;

      clearStagedCalculation();
      setActiveCalculationId(runningCalculation.id);

      return runningCalculation;
    } catch (error) {
      handleApiError(error, '計算の開始に失敗しました');
      throw error;
    }
  };

  const handleCalculationRename = async (
    calculationId: string,
    newName: string
  ): Promise<void> => {
    try {
      await updateCalculationNameMutation.mutateAsync({
        id: calculationId,
        newName,
      });
    } catch (error) {
      handleApiError(error, '計算名の変更に失敗しました');
    }
  };

  const handleCalculationDelete = async (
    calculationId: string
  ): Promise<void> => {
    try {
      await deleteCalculationMutation.mutateAsync(calculationId);
      if (activeCalculationId === calculationId) {
        setActiveCalculationId(null);
        clearStagedCalculation();
        setCurrentPage('calculation-settings');
      }
    } catch (error) {
      handleApiError(error, '計算の削除に失敗しました');
    }
  };

  return {
    handleStartCalculation,
    handleCalculationRename,
    handleCalculationDelete,
  };
};
