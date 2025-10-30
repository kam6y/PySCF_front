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
import { showErrorNotification } from '../store/notificationStore';

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
      showErrorNotification(defaultMessage, error.getUserMessage());
    } else if (error instanceof Error) {
      showErrorNotification(defaultMessage, error.message);
    } else {
      showErrorNotification(defaultMessage, 'An unknown error occurred.');
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
      handleApiError(error, 'Failed to start calculation');
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
      handleApiError(error, 'Failed to change calculation name');
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
      handleApiError(error, 'Failed to delete calculation');
    }
  };

  return {
    handleStartCalculation,
    handleCalculationRename,
    handleCalculationDelete,
  };
};
