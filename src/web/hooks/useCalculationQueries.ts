// src/web/hooks/useCalculationQueries.ts

import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import * as quantumApi from '../api/quantum';
import { searchPubChem, convertSmilesToXyz } from '../api/molecule';
import { QuantumCalculationRequest } from '../types/api-types';

// 計算リストを取得するQuery
export const useGetCalculations = () => {
  return useQuery({
    queryKey: ['calculations'],
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
    queryKey: ['calculation', id],
    queryFn: () => quantumApi.getCalculationDetails(id!),
    enabled: !!id && !id.startsWith('new-calculation-'), // idが存在し、一時IDでない場合にのみ実行

    // WebSocketがリアルタイム更新を提供するため、ポーリングは不要
    staleTime: 60 * 1000, // 1分 - WebSocketが主な更新メカニズム
    gcTime: 10 * 60 * 1000, // 10分 - 詳細データを長めに保持

    // ウィンドウフォーカス時の再フェッチを無効化（WebSocketが更新を管理）
    refetchOnWindowFocus: false,

    // ネットワーク復帰時の再フェッチは有効（WebSocketより先に復帰する可能性）
    refetchOnReconnect: true,

    // コンポーネント再マウント時は同期（WebSocket切断時の不整合を防ぐ）
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
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    },
  });
};

// 計算を削除するMutation
export const useDeleteCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (id: string) => quantumApi.deleteCalculation(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    },
  });
};

// 計算を一時停止するMutation
export const usePauseCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (id: string) => quantumApi.pauseCalculation(id),
    onSuccess: (data, id) => {
      // 成功したら関連するキャッシュを更新
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
      queryClient.invalidateQueries({ queryKey: ['calculation', id] });
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
      queryClient.setQueryData(['calculation', id], {
        calculation: data.calculation,
      });

      // リストキャッシュも更新
      queryClient.setQueryData(['calculations'], (oldData: any) => {
        if (!oldData?.calculations) return oldData;
        return {
          ...oldData,
          calculations: oldData.calculations.map((calc: any) =>
            calc.id === id ? data.calculation : calc
          ),
        };
      });
    },
  });
};

// 計算名を更新するMutation
export const useUpdateCalculationName = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: ({ id, newName }: { id: string; newName: string }) =>
      quantumApi.updateCalculationName(id, newName),
    onSuccess: (data, variables) => {
      // 成功したら関連するキャッシュを更新
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
      queryClient.invalidateQueries({
        queryKey: ['calculation', variables.id],
      });
    },
  });
};

// PubChem検索Mutation
export const useSearchPubChem = () => {
  return useMutation({
    mutationFn: ({
      query,
      searchType,
    }: {
      query: string;
      searchType: 'name' | 'cid';
    }) => searchPubChem(query, searchType),
  });
};

// SMILES変換Mutation
export const useConvertSmilesToXyz = () => {
  return useMutation({
    mutationFn: (smiles: string) => convertSmilesToXyz(smiles),
  });
};

// 軌道情報を取得するQuery
export const useGetOrbitals = (calculationId: string | null) => {
  return useQuery({
    queryKey: ['orbitals', calculationId],
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
    queryKey: ['orbital-cube', calculationId, orbitalIndex, options],
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

// 軌道のCUBEファイルを生成するMutation（再生成が必要な場合）
export const useGenerateOrbitalCube = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: ({
      calculationId,
      orbitalIndex,
      options,
    }: {
      calculationId: string;
      orbitalIndex: number;
      options?: {
        gridSize?: number;
        isovaluePos?: number;
        isovalueNeg?: number;
      };
    }) => quantumApi.getOrbitalCube(calculationId, orbitalIndex, options),
    onSuccess: (data, variables) => {
      // 成功したら該当するキャッシュを更新
      queryClient.setQueryData(
        [
          'orbital-cube',
          variables.calculationId,
          variables.orbitalIndex,
          variables.options,
        ],
        data
      );
    },
  });
};

// CUBE files management
export const useListCubeFiles = (calculationId: string | null) => {
  return useQuery({
    queryKey: ['cube-files', calculationId],
    queryFn: () => quantumApi.listCubeFiles(calculationId!),
    enabled: !!calculationId && !calculationId.startsWith('new-calculation-'),
  });
};

export const useDeleteCubeFiles = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: ({
      calculationId,
      orbitalIndex,
    }: {
      calculationId: string;
      orbitalIndex?: number;
    }) => quantumApi.deleteCubeFiles(calculationId, orbitalIndex),
    onSuccess: (data, variables) => {
      // Invalidate related queries
      queryClient.invalidateQueries({
        queryKey: ['cube-files', variables.calculationId],
      });
      queryClient.invalidateQueries({
        queryKey: ['orbital-cube', variables.calculationId],
      });
    },
  });
};

// サポートされているパラメータを取得するQuery
export const useSupportedParameters = () => {
  return useQuery({
    queryKey: ['supported-parameters'],
    queryFn: quantumApi.getSupportedParameters,
    staleTime: 24 * 60 * 60 * 1000, // 24時間キャッシュを保持（パラメータは頻繁に変更されない）
    gcTime: 24 * 60 * 60 * 1000, // 24時間メモリに保持
    refetchOnWindowFocus: false, // ウィンドウフォーカス時の再取得を無効化
    refetchOnMount: false, // マウント時の再取得を無効化（キャッシュを優先）
  });
};
