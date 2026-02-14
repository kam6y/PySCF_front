// React Query ベースのフック
export * from './useCalculationQueries';
export { useMethodDefaults } from './useMethodDefaults';

export { useUnifiedWebSocket } from './useUnifiedWebSocket';
export type { UseUnifiedWebSocketOptions } from './useUnifiedWebSocket';

// 新しい状態管理フック（Zustandベース）
export { useAppState } from './useAppState';
export { useCalculationData } from './useCalculationData';
export { useActiveCalculation } from './useActiveCalculation';
export { useCalculationActions } from './useCalculationActions';

// 設定管理フック
export {
  useAppSettings,
  useGetSettings,
  useUpdateSettings,
} from './useAppSettings';

// GPU4PySCF hooks
export {
  useGpu4Pyscf,
  useGpu4PyscfStatus,
  useInstallGpu4Pyscf,
} from './useGpu4Pyscf';

// パフォーマンス最適化フック
export { useProcessedCalculationResults } from './useProcessedCalculationResults';
