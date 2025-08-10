// src/web/apiClient.ts

import {
  CalculationDetailsResponse,
  CalculationListResponse,
  CalculationParameters,
  DeleteResponse,
  RenameResponse,
  CalculationInstance,
} from './types/calculation';

let API_BASE_URL = 'http://127.0.0.1:5000'; // Default, will be updated

export const setApiBaseUrl = (port: number) => {
    API_BASE_URL = `http://127.0.0.1:${port}`;
    console.log(`API base URL set to: ${API_BASE_URL}`);
};

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

// API Response for starting a calculation (now async)
export interface StartCalculationResponse {
    calculation: CalculationInstance;
}

type ApiResponse<T> = {
    success: boolean;
    data: T;
    error?: string;
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

  const apiResponse = await response.json() as ApiResponse<T>;

  if (!response.ok || !apiResponse.success) {
    throw new Error(apiResponse.error || 'APIとの通信中に不明なエラーが発生しました。');
  }

  return apiResponse.data; // Return the 'data' property directly
};

// --- API Functions ---

export const getCalculations = (): Promise<CalculationListResponse['data']> => {
  return request<CalculationListResponse['data']>('/api/quantum/calculations', { method: 'GET' });
};

export const getCalculationDetails = (id: string): Promise<CalculationDetailsResponse['data']> => {
  if (!id || id === 'undefined' || id === 'null') {
    return Promise.reject(new Error('Invalid calculation ID provided.'));
  }
  return request<CalculationDetailsResponse['data']>(`/api/quantum/calculations/${id}`, { method: 'GET' });
};

export const startCalculation = (params: CalculationParameters): Promise<StartCalculationResponse> => {
    return request<StartCalculationResponse>('/api/quantum/calculate', {
        method: 'POST',
        body: JSON.stringify(params),
    });
};

export const updateCalculationName = (id: string, newName: string): Promise<RenameResponse> => {
    return request<RenameResponse>(`/api/quantum/calculations/${id}`, {
      method: 'PUT',
      body: JSON.stringify({ name: newName }),
    });
};
  
export const deleteCalculation = (id: string): Promise<DeleteResponse> => {
    return request<DeleteResponse>(`/api/quantum/calculations/${id}`, {
      method: 'DELETE',
    });
};

export const searchPubChem = (query: string, searchType: 'name' | 'cid'): Promise<PubChemSearchResponse> => {
    return request<PubChemSearchResponse>('/api/pubchem/search', {
        method: 'POST',
        body: JSON.stringify({ query, search_type: searchType }),
    });
};

export const convertSmilesToXyz = (smiles: string): Promise<SmilesConvertResponse> => {
    return request<SmilesConvertResponse>('/api/smiles/convert', {
        method: 'POST',
        body: JSON.stringify({ smiles }),
    });
};