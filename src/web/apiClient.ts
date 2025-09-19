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
  OrbitalsResponseData,
  OrbitalCubeResponseData,
  CubeFilesListResponseData,
  CubeFilesDeleteResponseData,
  SupportedParametersResponseData,
  IRSpectrumResponseData,
} from './types/api-types';
import { components } from './types/generated-api';

let API_BASE_URL = 'http://127.0.0.1:5000'; // Default, will be updated

export const setApiBaseUrl = (port: number) => {
  API_BASE_URL = `http://127.0.0.1:${port}`;
  console.log(`API base URL set to: ${API_BASE_URL}`);
};

// Re-export the generated type for backward compatibility
export type StartCalculationResponse = StartCalculationResponseData;

// Settings API types
type AppSettings = components['schemas']['AppSettings'];
type SettingsResponse = components['schemas']['SettingsResponse'];

// System resource API types
type SystemResourceResponse = components['schemas']['SystemResourceResponse'];
type SystemResourceSummary = components['schemas']['SystemResourceSummary'];

type ApiResponse<T> = {
  success: boolean;
  data: T;
  error?: string;
};

/**
 * HTTPステータスコードやレスポンス詳細を保持するカスタムエラークラス
 */
export class ApiError extends Error {
  public readonly status: number;
  public readonly statusText: string;
  public readonly url: string;
  public readonly response?: any;
  public readonly isNetworkError: boolean;

  constructor(
    message: string,
    status: number,
    statusText: string,
    url: string,
    response?: any,
    isNetworkError = false
  ) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.statusText = statusText;
    this.url = url;
    this.response = response;
    this.isNetworkError = isNetworkError;

    // Maintains proper stack trace for where our error was thrown (only available on V8)
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, ApiError);
    }
  }

  /**
   * エラーのタイプを判定するヘルパーメソッド
   */
  get errorType(): 'network' | 'client' | 'server' | 'unknown' {
    if (this.isNetworkError) return 'network';
    if (this.status >= 400 && this.status < 500) return 'client';
    if (this.status >= 500) return 'server';
    return 'unknown';
  }

  /**
   * ユーザーフレンドリーなエラーメッセージを生成
   */
  getUserMessage(): string {
    switch (this.errorType) {
      case 'network':
        return 'A network connection error occurred. Please check your internet connection.';
      case 'client':
        if (this.status === 404) {
          return 'The requested resource was not found.';
        }
        if (this.status === 400) {
          return 'There is an issue with the request. Please check your input.';
        }
        if (this.status === 401) {
          return 'Authentication is required.';
        }
        if (this.status === 403) {
          return 'You do not have permission to access this resource.';
        }
        return 'A request error occurred.';
      case 'server':
        if (this.status === 503) {
          return 'The server is temporarily unavailable. Please try again later.';
        }
        return 'A server error occurred. Please contact the administrator.';
      default:
        return this.message || 'An unknown error occurred.';
    }
  }
}

/**
 * APIリクエストを処理する汎用関数
 */
const request = async <T>(
  endpoint: string,
  options: RequestInit = {}
): Promise<T> => {
  const url = `${API_BASE_URL}${endpoint}`;

  try {
    const response = await fetch(url, {
      headers: {
        'Content-Type': 'application/json',
        ...options.headers,
      },
      ...options,
    });

    let apiResponse: ApiResponse<T>;

    try {
      apiResponse = await response.json();
    } catch (jsonError) {
      // JSON解析エラー（レスポンスがJSONでない場合）
      throw new ApiError(
        'An invalid response was returned from the server.',
        response.status,
        response.statusText,
        url,
        null,
        false
      );
    }

    if (!response.ok || !apiResponse.success) {
      const errorMessage =
        apiResponse.error ||
        `HTTPエラー: ${response.status} ${response.statusText}`;

      throw new ApiError(
        errorMessage,
        response.status,
        response.statusText,
        url,
        apiResponse,
        false
      );
    }

    return apiResponse.data; // Return the 'data' property directly
  } catch (error) {
    // fetchそのものが失敗した場合（ネットワークエラーなど）
    if (error instanceof ApiError) {
      throw error; // 既にApiErrorの場合はそのまま投げる
    }

    // ネットワークエラーやその他のfetchエラー
    throw new ApiError(
      'A network error occurred. Unable to connect to the server.',
      0, // ネットワークエラーの場合はステータスコード0
      'Network Error',
      url,
      null,
      true
    );
  }
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
    return Promise.reject(
      new ApiError(
        'Invalid calculation ID provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${id}`,
        null,
        false
      )
    );
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
      method: 'PATCH',
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

/**
 * Get molecular orbital information for a calculation
 */
export const getOrbitals = (
  calculationId: string
): Promise<OrbitalsResponseData> => {
  if (
    !calculationId ||
    calculationId === 'undefined' ||
    calculationId === 'null'
  ) {
    return Promise.reject(
      new ApiError(
        'Invalid calculation ID provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${calculationId}/orbitals`,
        null,
        false
      )
    );
  }
  return request<OrbitalsResponseData>(
    `/api/quantum/calculations/${calculationId}/orbitals`,
    { method: 'GET' }
  );
};

/**
 * Generate and get CUBE file for a specific molecular orbital
 */
export const getOrbitalCube = (
  calculationId: string,
  orbitalIndex: number,
  options?: {
    gridSize?: number;
    isovaluePos?: number;
    isovalueNeg?: number;
  }
): Promise<OrbitalCubeResponseData> => {
  if (
    !calculationId ||
    calculationId === 'undefined' ||
    calculationId === 'null'
  ) {
    return Promise.reject(
      new ApiError(
        'Invalid calculation ID provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${calculationId}/orbitals/${orbitalIndex}/cube`,
        null,
        false
      )
    );
  }

  if (orbitalIndex < 0 || !Number.isInteger(orbitalIndex)) {
    return Promise.reject(
      new ApiError(
        'Invalid orbital index provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${calculationId}/orbitals/${orbitalIndex}/cube`,
        null,
        false
      )
    );
  }

  // Build query parameters
  const queryParams = new URLSearchParams();
  if (options?.gridSize !== undefined) {
    queryParams.append('gridSize', options.gridSize.toString());
  }
  if (options?.isovaluePos !== undefined) {
    queryParams.append('isovaluePos', options.isovaluePos.toString());
  }
  if (options?.isovalueNeg !== undefined) {
    queryParams.append('isovalueNeg', options.isovalueNeg.toString());
  }

  const queryString = queryParams.toString();
  const endpoint = `/api/quantum/calculations/${calculationId}/orbitals/${orbitalIndex}/cube${
    queryString ? `?${queryString}` : ''
  }`;

  return request<OrbitalCubeResponseData>(endpoint, { method: 'GET' });
};

/**
 * List all CUBE files for a calculation
 */
export const listCubeFiles = (
  calculationId: string
): Promise<CubeFilesListResponseData> => {
  if (
    !calculationId ||
    calculationId === 'undefined' ||
    calculationId === 'null'
  ) {
    return Promise.reject(
      new ApiError(
        'Invalid calculation ID provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${calculationId}/orbitals/cube-files`,
        null,
        false
      )
    );
  }

  return request<CubeFilesListResponseData>(
    `/api/quantum/calculations/${calculationId}/orbitals/cube-files`,
    { method: 'GET' }
  );
};

/**
 * Delete CUBE files for a calculation
 */
export const deleteCubeFiles = (
  calculationId: string,
  orbitalIndex?: number
): Promise<CubeFilesDeleteResponseData> => {
  if (
    !calculationId ||
    calculationId === 'undefined' ||
    calculationId === 'null'
  ) {
    return Promise.reject(
      new ApiError(
        'Invalid calculation ID provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${calculationId}/orbitals/cube-files`,
        null,
        false
      )
    );
  }

  // Build query parameters
  const queryParams = new URLSearchParams();
  if (orbitalIndex !== undefined && orbitalIndex >= 0) {
    queryParams.append('orbital_index', orbitalIndex.toString());
  }

  const queryString = queryParams.toString();
  const endpoint = `/api/quantum/calculations/${calculationId}/orbitals/cube-files${
    queryString ? `?${queryString}` : ''
  }`;

  return request<CubeFilesDeleteResponseData>(endpoint, { method: 'DELETE' });
};

/**
 * Get supported quantum chemistry parameters
 */
export const getSupportedParameters =
  (): Promise<SupportedParametersResponseData> => {
    return request<SupportedParametersResponseData>(
      '/api/quantum/supported-parameters',
      { method: 'GET' }
    );
  };

/**
 * Get IR spectrum for a calculation
 */
export const getIRSpectrum = (
  calculationId: string,
  options?: {
    broadening_fwhm?: number;
    x_min?: number;
    x_max?: number;
    show_peaks?: boolean;
  }
): Promise<IRSpectrumResponseData> => {
  if (
    !calculationId ||
    calculationId === 'undefined' ||
    calculationId === 'null'
  ) {
    return Promise.reject(
      new ApiError(
        'Invalid calculation ID provided.',
        400,
        'Bad Request',
        `/api/quantum/calculations/${calculationId}/ir-spectrum`,
        null,
        false
      )
    );
  }

  // Build query parameters
  const queryParams = new URLSearchParams();
  if (options?.broadening_fwhm !== undefined) {
    queryParams.append('broadening_fwhm', options.broadening_fwhm.toString());
  }
  if (options?.x_min !== undefined) {
    queryParams.append('x_min', options.x_min.toString());
  }
  if (options?.x_max !== undefined) {
    queryParams.append('x_max', options.x_max.toString());
  }
  if (options?.show_peaks !== undefined) {
    queryParams.append('show_peaks', options.show_peaks.toString());
  }

  const queryString = queryParams.toString();
  const endpoint = `/api/quantum/calculations/${calculationId}/ir-spectrum${
    queryString ? `?${queryString}` : ''
  }`;

  return request<IRSpectrumResponseData>(endpoint, { method: 'GET' });
};

/**
 * Get current application settings
 */
export const getSettings = (): Promise<SettingsResponse['data']> => {
  return request<SettingsResponse['data']>('/api/settings', { method: 'GET' });
};

/**
 * Update application settings
 */
export const updateSettings = (
  settings: AppSettings
): Promise<SettingsResponse['data']> => {
  return request<SettingsResponse['data']>('/api/settings', {
    method: 'PUT',
    body: JSON.stringify(settings),
  });
};

/**
 * Get current system resource status
 */
export const getSystemResourceStatus = (): Promise<SystemResourceResponse> => {
  return request<SystemResourceResponse>('/api/system/resource-status', {
    method: 'GET',
  });
};
