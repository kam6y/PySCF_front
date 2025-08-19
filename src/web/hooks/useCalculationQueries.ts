// src/web/hooks/useCalculationQueries.ts

import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import * as apiClient from '../apiClient';
import { QuantumCalculationRequest } from '../types/api-types';

// 計算リストを取得するQuery
export const useGetCalculations = () => {
  return useQuery({
    queryKey: ['calculations'],
    queryFn: apiClient.getCalculations,
  });
};

// 特定の計算詳細を取得するQuery
export const useGetCalculationDetails = (id: string | null) => {
  return useQuery({
    queryKey: ['calculation', id],
    queryFn: () => apiClient.getCalculationDetails(id!),
    enabled: !!id && !id.startsWith('new-calculation-'), // idが存在し、一時IDでない場合にのみ実行
  });
};

// 計算を開始するMutation
export const useStartCalculation = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (params: QuantumCalculationRequest) =>
      apiClient.startCalculation(params),
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
    mutationFn: (id: string) => apiClient.deleteCalculation(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    },
  });
};

// 計算名を更新するMutation
export const useUpdateCalculationName = () => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: ({ id, newName }: { id: string; newName: string }) =>
      apiClient.updateCalculationName(id, newName),
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
    }) => apiClient.searchPubChem(query, searchType),
  });
};

// SMILES変換Mutation
export const useConvertSmilesToXyz = () => {
  return useMutation({
    mutationFn: (smiles: string) => apiClient.convertSmilesToXyz(smiles),
  });
};

// 軌道情報を取得するQuery
export const useGetOrbitals = (calculationId: string | null) => {
  return useQuery({
    queryKey: ['orbitals', calculationId],
    queryFn: () => apiClient.getOrbitals(calculationId!),
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
      apiClient.getOrbitalCube(calculationId!, orbitalIndex!, options),
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
    }) => apiClient.getOrbitalCube(calculationId, orbitalIndex, options),
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
    queryFn: () => apiClient.listCubeFiles(calculationId!),
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
    }) => apiClient.deleteCubeFiles(calculationId, orbitalIndex),
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
