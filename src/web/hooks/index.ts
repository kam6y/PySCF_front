// React Query ベースのフック
export * from './useCalculationQueries';

export { useUnifiedWebSocket } from './useUnifiedWebSocket';
export type { UseUnifiedWebSocketOptions } from './useUnifiedWebSocket';

// 新しい状態管理フック（Zustandベース）
export { useAppState } from './useAppState';
export { useCalculationData } from './useCalculationData';
export { useActiveCalculation } from './useActiveCalculation';
export { useCalculationActions } from './useCalculationActions';

// その他のフック
export { useActiveCalculationId } from './useActiveCalculationId';

// 設定管理フック
export {
  useAppSettings,
  useGetSettings,
  useUpdateSettings,
} from './useAppSettings';

// 互換性のため残存（非推奨）
export { useCalculationOperations } from './useCalculationOperations';
export type { CalculationOperations } from './useCalculationOperations';

// パフォーマンス最適化フック
export { useProcessedCalculationResults } from './useProcessedCalculationResults';
