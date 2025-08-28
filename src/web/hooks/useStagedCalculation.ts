import { useCalculationStore } from '../store/calculationStore';
import {
  CalculationInstance,
  QuantumCalculationRequest,
} from '../types/api-types';

export const useStagedCalculation = () => {
  const {
    stagedCalculation,
    setStagedCalculation,
    clearStagedCalculation,
    setActiveCalculationId,
  } = useCalculationStore();

  const createNewCalculation = (
    setCurrentPage: (page: 'calculation-settings') => void,
    handleSidebarClose: () => void
  ) => {
    const defaultParams: QuantumCalculationRequest = {
      calculation_method: 'DFT',
      basis_function: '6-31G(d)',
      exchange_correlation: 'B3LYP',
      charges: 0,
      spin_multiplicity: 0,
      solvent_method: 'none',
      solvent: '-',
      xyz: '',
      name: '',
      tddft_nstates: 10,
      tddft_method: 'TDDFT',
      tddft_analyze_nto: false,
    };

    const newId = `new-calculation-${Date.now()}`;
    const newCalculation: CalculationInstance = {
      id: newId,
      name: '',
      status: 'pending',
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      parameters: defaultParams,
      results: undefined,
    };

    setStagedCalculation(newCalculation);
    setActiveCalculationId(newId);
    setCurrentPage('calculation-settings');
    handleSidebarClose();
  };

  const createNewFromExisting = (
    originalCalc: CalculationInstance,
    newParams: QuantumCalculationRequest
  ) => {
    const newId = `new-calculation-${Date.now()}`;
    const newCalculation: CalculationInstance = {
      ...originalCalc,
      id: newId,
      status: 'pending',
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      parameters: newParams,
      results: undefined,
      workingDirectory: undefined,
      errorMessage: undefined,
    };

    setStagedCalculation(newCalculation);
    setActiveCalculationId(newId);
  };

  const updateStagedCalculation = (updatedCalculation: CalculationInstance) => {
    if (stagedCalculation && updatedCalculation.id === stagedCalculation.id) {
      setStagedCalculation(updatedCalculation);
    }
  };

  const clearStaged = () => {
    clearStagedCalculation();
  };

  return {
    stagedCalculation,
    createNewCalculation,
    createNewFromExisting,
    updateStagedCalculation,
    clearStaged,
    isStagedCalculation: (id: string | null) =>
      id ? id.startsWith('new-calculation-') : false,
  };
};
