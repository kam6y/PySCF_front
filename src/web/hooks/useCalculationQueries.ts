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
