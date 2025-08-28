import { useEffect } from 'react';
import { useCalculationStore } from '../store/calculationStore';
import { useGetCalculations } from './useCalculationQueries';

export const useActiveCalculationId = () => {
  const {
    activeCalculationId,
    stagedCalculation,
    setActiveCalculationId,
  } = useCalculationStore();

  const { data: calculationsData } = useGetCalculations();

  // 自動選択ロジック：最初の計算を選択（新規計算作成中は除く）
  useEffect(() => {
    if (!activeCalculationId && 
        calculationsData?.calculations && 
        calculationsData.calculations.length > 0 && 
        !stagedCalculation) {
      setActiveCalculationId(calculationsData.calculations[0].id);
    }
  }, [
    calculationsData,
    activeCalculationId,
    setActiveCalculationId,
    stagedCalculation,
  ]);

  const selectCalculation = (calculationId: string | null) => {
    setActiveCalculationId(calculationId);
  };

  return {
    activeCalculationId,
    selectCalculation,
  };
};