import { useMemo } from 'react';
import {
  useGetCalculations,
  useGetCalculationDetails,
} from './useCalculationQueries';
import { CalculationInstance } from '../types/api-types';

export const useCalculationData = (activeCalculationId: string | null) => {
  const {
    data: calculationsData,
    isLoading: calculationsLoading,
    error: calculationsError,
  } = useGetCalculations();

  const { data: activeCalculationDetails, isLoading: detailsLoading } =
    useGetCalculationDetails(activeCalculationId);

  // 計算リストのデータ変換（必要な場合のみ）
  const calculations = useMemo(() => 
    calculationsData?.calculations?.map(c => ({
      id: c.id,
      name: c.name,
      status: c.status,
      createdAt: new Date(c.date).toISOString(),
      updatedAt: new Date(c.date).toISOString(),
      parameters: {} as any,
      results: undefined,
    })) || [], [calculationsData]);

  // サイドバー用の計算リスト（フィルタリング）
  const sidebarCalculations = useMemo(() => 
    calculations.filter(
      calc =>
        calc.status !== 'pending' ||
        (calc.parameters &&
          calc.parameters.xyz &&
          calc.parameters.xyz.trim() !== '')
    ), [calculations]);

  // アクティブ計算の詳細データ（詳細が取得できている場合はそれを使用）
  const activeCalculationDetail = activeCalculationDetails?.calculation;

  // アクティブ計算の基本情報（リストから取得）
  const activeCalculationBasic = calculations.find(c => c.id === activeCalculationId);

  return {
    // Raw data from React Query
    calculations,
    sidebarCalculations,
    activeCalculationDetail,
    activeCalculationBasic,
    
    // Loading states
    calculationsLoading,
    detailsLoading,
    
    // Errors
    calculationsError,
  };
};