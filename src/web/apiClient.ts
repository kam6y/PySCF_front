// src/web/apiClient.ts

import {
  CalculationDetailsResponse,
  CalculationListResponse,
  CalculationParameters,
} from './types/calculation';

const API_BASE_URL = 'http://127.0.0.1:5000';

// --- Response Type Definitions ---
interface PubChemCompoundInfo {
    cid: number;
    iupac_name: string;
    molecular_formula: string;
    molecular_weight: number;
    synonyms: string[];
}

export interface PubChemSearchResponse {
    xyz: string;
    compound_info: PubChemCompoundInfo;
}

export interface SmilesConvertResponse {
    xyz: string;
}


/**
 * APIリクエストを処理する汎用関数
 */
const request = async <T>(endpoint: string, options: RequestInit = {}): Promise<T> => {
  const response = await fetch(`${API_BASE_URL}${endpoint}`, {
    headers: {
      'Content-Type': 'application/json',
      ...options.headers,
    },
    ...options,
  });

  const data = await response.json();

  if (!response.ok || !data.success) {
    throw new Error(data.error || 'APIとの通信中に不明なエラーが発生しました。');
  }

  return data;
};

// --- API Functions ---

export const getCalculations = (): Promise<CalculationListResponse> => {
  return request<CalculationListResponse>('/api/quantum/calculations', { method: 'GET' });
};

export const getCalculationDetails = (id: string): Promise<CalculationDetailsResponse> => {
  return request<CalculationDetailsResponse>(`/api/quantum/calculations/${id}`, { method: 'GET' });
};

export const startCalculation = (params: CalculationParameters): Promise<{ calculation_id: string; calculation_results: any; }> => {
    return request<{ data: { calculation_id: string; calculation_results: any; } }>('/api/quantum/calculate', {
        method: 'POST',
        body: JSON.stringify(params),
    }).then(res => res.data);
};

export const searchPubChem = (query: string, searchType: 'name' | 'cid'): Promise<PubChemSearchResponse> => {
    return request<{ data: PubChemSearchResponse }>('/api/pubchem/search', {
        method: 'POST',
        body: JSON.stringify({ query, search_type: searchType }),
    }).then(res => res.data);
};

export const convertSmilesToXyz = (smiles: string): Promise<SmilesConvertResponse> => {
    return request<{ data: SmilesConvertResponse }>('/api/smiles/convert', {
        method: 'POST',
        body: JSON.stringify({ smiles }),
    }).then(res => res.data);
};