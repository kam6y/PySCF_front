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
import { handleError } from '../utils/errorHandler';

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
      handleError(error, 'Failed to start calculation');
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
      handleError(error, 'Failed to change calculation name');
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
      handleError(error, 'Failed to delete calculation');
    }
  };

  return {
    handleStartCalculation,
    handleCalculationRename,
    handleCalculationDelete,
  };
};
