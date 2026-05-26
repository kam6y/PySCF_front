// src/web/hooks/useCalculationQueries.ts

import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import * as quantumApi from '../api/quantum';
import {
  CalculationInstance,
  CalculationListResponseData,
  CalculationSummary,
  QuantumCalculationRequest,
} from '../types/api-types';

export const calculationQueryKeys = {
  all: ['calculations'] as const,
  list: () => [...calculationQueryKeys.all, 'list'] as const,
  detail: (id: string) => [...calculationQueryKeys.all, 'detail', id] as const,
  orbitals: (id: string) =>
    [...calculationQueryKeys.all, 'orbitals', id] as const,
  orbitalCube: (id: string, idx: number, opts?: object) =>
    [...calculationQueryKeys.all, 'orbital-cube', id, idx, opts] as const,
  supportedParams: () =>
    [...calculationQueryKeys.all, 'supported-parameters'] as const,
};

const toCalculationSummary = (
  calculation: CalculationInstance,
  previousSummary: CalculationSummary
): CalculationSummary => ({
  ...previousSummary,
  id: calculation.id,
  name: calculation.name,
  status: calculation.status,
  date: previousSummary.date || calculation.createdAt,
});

// 計算リストを取得するQuery
export const useGetCalculations = () => {
  return useQuery({
    queryKey: calculationQueryKeys.list(),
    queryFn: quantumApi.getCalculations,

    // リストは頻繁に変更される可能性があるため、staleTimeを短めに
    staleTime: 30 * 1000, // 30秒
    gcTime: 5 * 60 * 1000, // 5分

    // リストはフォーカス時に再フェッチすると便利
    refetchOnWindowFocus: true,
    refetchOnReconnect: true,
  });
};

// 特定の計算詳細を取得するQuery
export const useGetCalculationDetails = (id: string | null) => {
  return useQuery({
    queryKey: calculationQueryKeys.detail(id ?? ''),
    queryFn: () => quantumApi.getCalculationDetails(id!),
    enabled: !!id && !id.startsWith('new-calculation-'), // idが存在し、一時IDでない場合にのみ実行

    // SSE realtime update streamが更新を提供するため、ポーリングは不要
    staleTime: 60 * 1000, // 1分 - realtime update streamが主な更新メカニズム
    gcTime: 10 * 60 * 1000, // 10分 - 詳細データを長めに保持

    // ウィンドウフォーカス時の再フェッチを無効化（realtime update streamが更新を管理）
    refetchOnWindowFocus: false,

    // ネットワーク復帰時の再フェッチは有効（realtime update streamより先に復帰する可能性）
    refetchOnReconnect: true,

    // コンポーネント再マウント時は同期（realtime update stream切断時の不整合を防ぐ）
    refetchOnMount: true,
  });
};

// 計算を開始するMutation
export const useStartCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (params: QuantumCalculationRequest) =>
      quantumApi.startCalculation(params),
    onSuccess: () => {
      // 成功したら計算リストのキャッシュを無効化して再取得させる
      queryClient.invalidateQueries({ queryKey: calculationQueryKeys.list() });
    },
  });
};

// 計算を削除するMutation
export const useDeleteCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (id: string) => quantumApi.deleteCalculation(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: calculationQueryKeys.list() });
    },
  });
};

// 計算を一時停止するMutation
export const usePauseCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (id: string) => quantumApi.pauseCalculation(id),
    onSuccess: (_data, id) => {
      // 成功したら関連するキャッシュを更新
      queryClient.invalidateQueries({ queryKey: calculationQueryKeys.list() });
      queryClient.invalidateQueries({
        queryKey: calculationQueryKeys.detail(id),
      });
    },
  });
};

// 計算を再開するMutation
export const useResumeCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (id: string) => quantumApi.resumeCalculation(id),
    onSuccess: (data, id) => {
      // サーバーレスポンスを即座にキャッシュに反映
      queryClient.setQueryData(calculationQueryKeys.detail(id), {
        calculation: data.calculation,
      });

      // リストキャッシュも更新
      queryClient.setQueryData(
        calculationQueryKeys.list(),
        (oldData: CalculationListResponseData | undefined) => {
          if (!oldData?.calculations) return oldData;
          return {
            ...oldData,
            calculations: oldData.calculations.map(
              (calc: CalculationListResponseData['calculations'][number]) =>
                calc.id === id
                  ? toCalculationSummary(data.calculation, calc)
                  : calc
            ),
          };
        }
      );
    },
  });
};

// 計算名を更新するMutation
export const useUpdateCalculationName = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: ({ id, newName }: { id: string; newName: string }) =>
      quantumApi.updateCalculationName(id, newName),
    onSuccess: (_data, variables) => {
      // 成功したら関連するキャッシュを更新
      queryClient.invalidateQueries({ queryKey: calculationQueryKeys.list() });
      queryClient.invalidateQueries({
        queryKey: calculationQueryKeys.detail(variables.id),
      });
    },
  });
};

// 軌道情報を取得するQuery
export const useGetOrbitals = (calculationId: string | null) => {
  return useQuery({
    queryKey: calculationQueryKeys.orbitals(calculationId ?? ''),
    queryFn: () => quantumApi.getOrbitals(calculationId!),
    enabled: !!calculationId && !calculationId.startsWith('new-calculation-'), // idが存在し、一時IDでない場合にのみ実行
  });
};

// 軌道のCUBEファイルを取得するQuery
export const useGetOrbitalCube = (
  calculationId: string | null,
  orbitalIndex: number | null,
  options?: {
    gridSize?: number;
    isovaluePos?: number;
    isovalueNeg?: number;
  }
) => {
  return useQuery({
    queryKey: calculationQueryKeys.orbitalCube(
      calculationId ?? '',
      orbitalIndex ?? -1,
      options
    ),
    queryFn: () =>
      quantumApi.getOrbitalCube(calculationId!, orbitalIndex!, options),
    enabled:
      !!calculationId &&
      orbitalIndex !== null &&
      orbitalIndex >= 0 &&
      !calculationId.startsWith('new-calculation-'),
    staleTime: 5 * 60 * 1000, // 5分間キャッシュを保持（計算切り替え時の更新を確保）
    gcTime: 24 * 60 * 60 * 1000, // 24時間メモリに保持（永続化されたファイルアクセス用）
    refetchOnWindowFocus: false, // ウィンドウフォーカス時の再取得を無効化
    refetchOnMount: false, // マウント時の再取得を無効化（staleTimeを優先してキャッシュを活用）
  });
};

// サポートされているパラメータを取得するQuery
export const useSupportedParameters = () => {
  return useQuery({
    queryKey: calculationQueryKeys.supportedParams(),
    queryFn: quantumApi.getSupportedParameters,
    staleTime: 24 * 60 * 60 * 1000, // 24時間キャッシュを保持（パラメータは頻繁に変更されない）
    gcTime: 24 * 60 * 60 * 1000, // 24時間メモリに保持
    refetchOnWindowFocus: false, // ウィンドウフォーカス時の再取得を無効化
    refetchOnMount: false, // マウント時の再取得を無効化（キャッシュを優先）
  });
};
