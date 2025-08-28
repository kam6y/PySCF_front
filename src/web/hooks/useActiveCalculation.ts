import { useMemo } from 'react';
import { useCalculationStore } from '../store/calculationStore';
import { useGetCalculations, useGetCalculationDetails } from './useCalculationQueries';
import { CalculationInstance } from '../types/api-types';

/**
 * アクティブ計算の統一された状態管理フック
 * 
 * 複雑な導出ロジックを内部で処理し、外部には単一のactiveCalculationインターフェースを提供
 * 
 * 優先順位:
 * 1. ステージド計算（新規作成中）
 * 2. サーバーからの詳細データ（永続化済み）
 * 3. 計算リストからの基本データ（フォールバック）
 */
export const useActiveCalculation = () => {
  const { activeCalculationId, stagedCalculation, setActiveCalculationId } = useCalculationStore();
  
  // サーバーデータの取得
  const { data: calculationsData, isLoading: calculationsLoading, error: calculationsError } = useGetCalculations();
  const { data: detailsData, isLoading: detailsLoading } = useGetCalculationDetails(activeCalculationId);

  // 計算リストからの基本データは詳細データの代替としては使用しない
  // （型の不整合を避けるため、詳細データまたはステージドデータのみを使用）
  const basicCalculation = null;

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
    
    // 3. フォールバック: リストからの基本データ
    return basicCalculation;
  }, [stagedCalculation, detailsData?.calculation, basicCalculation]);

  // 計算の種類を判定
  const calculationType = useMemo(() => {
    if (stagedCalculation) return 'staged';
    if (detailsData?.calculation) return 'detailed';
    if (basicCalculation) return 'basic';
    return 'none';
  }, [stagedCalculation, detailsData?.calculation, basicCalculation]);

  // ステージド計算かどうか
  const isStagedCalculation = useMemo(() => 
    calculationType === 'staged', 
    [calculationType]
  );

  // ローディング状態
  const isLoading = calculationsLoading || (activeCalculationId && !stagedCalculation && detailsLoading);

  // サイドバー用の計算リスト（サーバーデータをそのまま使用）
  const sidebarCalculations = useMemo(() => {
    if (!calculationsData?.calculations) return [];
    
    // サーバーから来るデータをフィルタリングのみ行い、変換は行わない
    return calculationsData.calculations.filter(calc => 
      calc.status !== 'pending' ||
      // pending状態でも表示する条件があればここに追加
      true
    );
  }, [calculationsData]);

  // 計算選択関数
  const selectCalculation = (calculationId: string | null) => {
    setActiveCalculationId(calculationId);
  };

  return {
    // メイン状態
    activeCalculation,
    activeCalculationId,
    calculationType,
    isStagedCalculation,
    
    // アクション
    selectCalculation,
    
    // ローディングとエラー状態
    isLoading,
    detailsLoading,
    calculationsLoading,
    calculationsError,
    
    // 補助データ
    sidebarCalculations,
  };
};