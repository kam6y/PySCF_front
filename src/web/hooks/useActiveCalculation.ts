import { useEffect } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import {
  useGetCalculations,
  useGetCalculationDetails,
} from './useCalculationQueries';
import { useCalculationStore } from '../store/calculationStore';
import {
  CalculationInstance,
  QuantumCalculationRequest,
} from '../types/api-types';

export interface ActiveCalculation {
  calculations: CalculationInstance[];
  sidebarCalculations: CalculationInstance[];
  activeCalculation: CalculationInstance | undefined;
  calculationsLoading: boolean;
  calculationsError: Error | null;
  handleCalculationSelect: (calculationId: string) => void;
  handleNewCalculation: () => void;
  handleActiveCalculationUpdate: (
    updatedCalculation: CalculationInstance
  ) => void;
  handleCreateNewFromExisting: (
    originalCalc: CalculationInstance,
    newParams: QuantumCalculationRequest
  ) => void;
}

export const useActiveCalculation = (
  setCurrentPage: (page: 'calculation-settings') => void,
  handleSidebarClose: () => void
): ActiveCalculation => {
  const queryClient = useQueryClient();

  const {
    activeCalculationId,
    stagedCalculation,
    setActiveCalculationId,
    setStagedCalculation,
    clearStagedCalculation,
  } = useCalculationStore();

  const {
    data: calculationsData,
    isLoading: calculationsLoading,
    error: calculationsError,
  } = useGetCalculations();
  const { data: activeCalculationDetails } =
    useGetCalculationDetails(activeCalculationId);

  // 計算リストのデータ変換
  const calculations =
    calculationsData?.calculations.map(c => ({
      id: c.id,
      name: c.name,
      status: c.status,
      createdAt: new Date(c.date).toISOString(),
      updatedAt: new Date(c.date).toISOString(),
      parameters: {} as any,
      results: undefined,
    })) || [];

  // アクティブな計算（詳細データがあればそれを使用、なければリストから取得）
  const activeCalculation =
    stagedCalculation ||
    activeCalculationDetails?.calculation ||
    calculations.find(c => c.id === activeCalculationId);

  // サイドバー用の計算リストのフィルタリング
  const sidebarCalculations = calculations.filter(
    calc =>
      calc.status !== 'pending' ||
      (calc.parameters &&
        calc.parameters.xyz &&
        calc.parameters.xyz.trim() !== '')
  );

  // 最初の計算を自動選択（新規計算作成中は除く）
  useEffect(() => {
    if (!activeCalculationId && calculations.length > 0 && !stagedCalculation) {
      setActiveCalculationId(calculations[0].id);
    }
  }, [
    calculations,
    activeCalculationId,
    setActiveCalculationId,
    stagedCalculation,
  ]);

  const handleCalculationSelect = (calculationId: string) => {
    setActiveCalculationId(calculationId);
    clearStagedCalculation();
    handleSidebarClose();
  };

  const handleNewCalculation = () => {
    const defaultParams: QuantumCalculationRequest = {
      calculation_method: 'DFT',
      basis_function: '6-31G(d)',
      exchange_correlation: 'B3LYP',
      charges: 0,
      spin_multiplicity: 1,
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

  const handleActiveCalculationUpdate = (
    updatedCalculation: CalculationInstance
  ) => {
    if (stagedCalculation && updatedCalculation.id === stagedCalculation.id) {
      setStagedCalculation(updatedCalculation);
    } else {
      // Queryキャッシュを更新
      queryClient.setQueryData(['calculation', updatedCalculation.id], {
        calculation: updatedCalculation,
      });
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    }
  };

  const handleCreateNewFromExisting = (
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

  return {
    calculations,
    sidebarCalculations,
    activeCalculation,
    calculationsLoading,
    calculationsError,
    handleCalculationSelect,
    handleNewCalculation,
    handleActiveCalculationUpdate,
    handleCreateNewFromExisting,
  };
};
