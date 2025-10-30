import { useCalculationData } from './useCalculationData';
import { useCalculationStore } from '../store/calculationStore';

/**
 * アクティブ計算への簡素化されたアクセス
 *
 * 複雑なロジックはuseCalculationDataに移動し、
 * このフックは主要なアクティブ計算情報への簡単なアクセスを提供
 */
export const useActiveCalculation = () => {
  const calculationData = useCalculationData();
  const { selectCalculation } = useCalculationStore();

  return {
    // メイン状態（useCalculationDataから取得）
    activeCalculation: calculationData.activeCalculation,
    activeCalculationId: calculationData.activeCalculationId,
    calculationType: calculationData.calculationType,
    isStagedCalculation: calculationData.isStagedCalculation,

    // アクション
    selectCalculation,

    // ローディングとエラー状態
    isLoading: calculationData.isLoading,
    detailsLoading: calculationData.detailsLoading,
    calculationsLoading: calculationData.calculationsLoading,
    calculationsError: calculationData.calculationsError,

    // 補助データ
    sidebarCalculations: calculationData.sidebarCalculations,
  };
};
