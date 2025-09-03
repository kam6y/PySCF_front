import { useMemo } from 'react';
import { 
  useGetCalculations, 
  useGetCalculationDetails 
} from './useCalculationQueries';
import { useCalculationStore } from '../store/calculationStore';

/**
 * 計算データ取得専用フック
 * 
 * TanStack Queryを使用してサーバーからの計算データを管理
 * データフェッチングのロジックを単一の場所に集約
 */
export const useCalculationData = () => {
  const { activeCalculationId, stagedCalculation } = useCalculationStore();

  // サーバーデータの取得
  const {
    data: calculationsData,
    isLoading: calculationsLoading,
    error: calculationsError,
  } = useGetCalculations();

  const { 
    data: detailsData, 
    isLoading: detailsLoading 
  } = useGetCalculationDetails(activeCalculationId);

  // アクティブ計算の決定（優先順位付き）
  const activeCalculation = useMemo(() => {
    // 1. ステージド計算（新規作成中）が最優先
    if (stagedCalculation) {
      return stagedCalculation;
    }

    // 2. サーバーからの詳細データ
    if (detailsData?.calculation) {
      return detailsData.calculation;
    }

    // 3. フォールバック: null
    return null;
  }, [stagedCalculation, detailsData?.calculation]);

  // 計算の種類を判定
  const calculationType = useMemo(() => {
    if (stagedCalculation) return 'staged';
    if (detailsData?.calculation) return 'detailed';
    return 'none';
  }, [stagedCalculation, detailsData?.calculation]);

  // ステージド計算かどうか
  const isStagedCalculation = useMemo(
    () => calculationType === 'staged',
    [calculationType]
  );

  // ローディング状態
  const isLoading = calculationsLoading || 
    (activeCalculationId && !stagedCalculation && detailsLoading);

  // サイドバー用の計算リスト
  const sidebarCalculations = useMemo(() => {
    if (!calculationsData?.calculations) return [];

    return calculationsData.calculations.filter(
      calc => calc.status !== 'pending' || true // 必要に応じてフィルタリング条件を調整
    );
  }, [calculationsData]);

  return {
    // メイン状態
    activeCalculation,
    activeCalculationId,
    calculationType,
    isStagedCalculation,

    // ローディングとエラー状態
    isLoading,
    detailsLoading,
    calculationsLoading,
    calculationsError,

    // 補助データ
    sidebarCalculations,
    
    // 生のデータ（必要に応じて）
    calculationsData,
    detailsData,
  };
};