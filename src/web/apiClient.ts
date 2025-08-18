// src/web/apiClient.ts

import {
  CalculationDetailsResponseData,
  CalculationListResponseData,
  QuantumCalculationRequest,
  CalculationUpdateResponseData,
  CalculationDeleteResponseData,
  PubChemSearchResponseData,
  SMILESConvertResponseData,
  StartCalculationResponseData,
} from './types/api-types';

let API_BASE_URL = 'http://127.0.0.1:5000'; // Default, will be updated

export const setApiBaseUrl = (port: number) => {
  API_BASE_URL = `http://127.0.0.1:${port}`;
  console.log(`API base URL set to: ${API_BASE_URL}`);
};

// ===== この関数を追加 =====
export const getWebSocketUrl = (path: string): string => {
  // API_BASE_URL (http://127.0.0.1:xxxx) から ws://... を生成
  const wsUrl = API_BASE_URL.replace(/^http/, 'ws');
  return `${wsUrl}${path}`;
};

// Re-export the generated type for backward compatibility
export type StartCalculationResponse = StartCalculationResponseData;

type ApiResponse<T> = {
  success: boolean;
  data: T;
  error?: string;
};

/**
 * APIリクエストを処理する汎用関数
 */
const request = async <T>(
  endpoint: string,
  options: RequestInit = {}
): Promise<T> => {
  const response = await fetch(`${API_BASE_URL}${endpoint}`, {
    headers: {
      'Content-Type': 'application/json',
      ...options.headers,
    },
    ...options,
  });

  const apiResponse = (await response.json()) as ApiResponse<T>;

  if (!response.ok || !apiResponse.success) {
    throw new Error(
      apiResponse.error || 'APIとの通信中に不明なエラーが発生しました。'
    );
  }

  return apiResponse.data; // Return the 'data' property directly
};

// --- API Functions ---

export const getCalculations = (): Promise<CalculationListResponseData> => {
  return request<CalculationListResponseData>('/api/quantum/calculations', {
    method: 'GET',
  });
};

export const getCalculationDetails = (
  id: string
): Promise<CalculationDetailsResponseData> => {
  if (!id || id === 'undefined' || id === 'null') {
    return Promise.reject(new Error('Invalid calculation ID provided.'));
  }
  return request<CalculationDetailsResponseData>(
    `/api/quantum/calculations/${id}`,
    { method: 'GET' }
  );
};

export const startCalculation = (
  params: QuantumCalculationRequest
): Promise<StartCalculationResponse> => {
  return request<StartCalculationResponse>('/api/quantum/calculate', {
    method: 'POST',
    body: JSON.stringify(params),
  });
};

export const updateCalculationName = (
  id: string,
  newName: string
): Promise<CalculationUpdateResponseData> => {
  return request<CalculationUpdateResponseData>(
    `/api/quantum/calculations/${id}`,
    {
      method: 'PUT',
      body: JSON.stringify({ name: newName }),
    }
  );
};

export const deleteCalculation = (
  id: string
): Promise<CalculationDeleteResponseData> => {
  return request<CalculationDeleteResponseData>(
    `/api/quantum/calculations/${id}`,
    {
      method: 'DELETE',
    }
  );
};

export const searchPubChem = (
  query: string,
  searchType: 'name' | 'cid'
): Promise<PubChemSearchResponseData> => {
  return request<PubChemSearchResponseData>('/api/pubchem/search', {
    method: 'POST',
    body: JSON.stringify({ query, search_type: searchType }),
  });
};

export const convertSmilesToXyz = (
  smiles: string
): Promise<SMILESConvertResponseData> => {
  return request<SMILESConvertResponseData>('/api/smiles/convert', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });
};
