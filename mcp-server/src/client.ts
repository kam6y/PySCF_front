import axios, { AxiosInstance, AxiosResponse, AxiosError } from 'axios';
import type { components } from './types.js';

export type ApiResponse<T> = {
  success: boolean;
  data?: T;
  error?: string;
};

export interface DetailedErrorInfo {
  message: string;
  status?: number;
  statusText?: string;
  url?: string;
  method?: string;
  requestData?: any;
  responseData?: any;
  timestamp: string;
}

export class PySCFApiError extends Error {
  public readonly details: DetailedErrorInfo;

  constructor(details: DetailedErrorInfo) {
    super(details.message);
    this.name = 'PySCFApiError';
    this.details = details;
  }
}

export type PubChemSearchRequest = components['schemas']['PubChemSearchRequest'];
export type PubChemSearchResponse = components['schemas']['PubChemSearchResponse'];
export type SMILESConvertRequest = components['schemas']['SMILESConvertRequest'];
export type SMILESConvertResponse = components['schemas']['SMILESConvertResponse'];
export type XYZValidateRequest = components['schemas']['XYZValidateRequest'];
export type XYZValidateResponse = components['schemas']['XYZValidateResponse'];
export type QuantumCalculationRequest = components['schemas']['QuantumCalculationRequest'];
export type StartCalculationResponse = components['schemas']['StartCalculationResponse'];
export type CalculationListResponse = components['schemas']['CalculationListResponse'];
export type CalculationDetailsResponse = components['schemas']['CalculationDetailsResponse'];
export type CalculationUpdateRequest = components['schemas']['CalculationUpdateRequest'];
export type CalculationUpdateResponse = components['schemas']['CalculationUpdateResponse'];
export type CalculationDeleteResponse = components['schemas']['CalculationDeleteResponse'];
export type OrbitalsResponse = components['schemas']['OrbitalsResponse'];
export type OrbitalCubeResponse = components['schemas']['OrbitalCubeResponse'];
export type CubeFilesListResponse = components['schemas']['CubeFilesListResponse'];
export type CubeFilesDeleteResponse = components['schemas']['CubeFilesDeleteResponse'];
export type IRSpectrumResponse = components['schemas']['IRSpectrumResponse'];
export type SupportedParametersResponse = components['schemas']['SupportedParametersResponse'];
export type SettingsResponse = components['schemas']['SettingsResponse'];
export type SettingsUpdateRequest = components['schemas']['SettingsUpdateRequest'];
export type SystemResourceResponse = components['schemas']['SystemResourceResponse'];
export type HealthResponse = components['schemas']['HealthResponse'];

export class PySCFApiClient {
  private client: AxiosInstance;
  private baseUrl: string;

  /**
   * Extract detailed error information from Axios error
   */
  private extractErrorDetails(error: any, requestData?: any): DetailedErrorInfo {
    const timestamp = new Date().toISOString();
    
    if (axios.isAxiosError(error)) {
      const axiosError = error as AxiosError;
      const details: DetailedErrorInfo = {
        message: error.message || 'Unknown API error',
        timestamp,
        requestData: requestData,
        responseData: axiosError.response?.data,
      };
      
      if (axiosError.response?.status !== undefined) details.status = axiosError.response.status;
      if (axiosError.response?.statusText !== undefined) details.statusText = axiosError.response.statusText;
      if (axiosError.config?.url !== undefined) details.url = axiosError.config.url;
      if (axiosError.config?.method !== undefined) details.method = axiosError.config.method.toUpperCase();
      
      return details;
    }

    return {
      message: error instanceof Error ? error.message : String(error),
      timestamp
    };
  }

  constructor(baseUrl: string = 'http://127.0.0.1:5000') {
    this.baseUrl = baseUrl;
    this.client = axios.create({
      baseURL: baseUrl,
      timeout: 30000, // 30 second timeout for most requests
      headers: {
        'Content-Type': 'application/json',
      },
    });
  }

  /**
   * Auto-detect available port by checking health endpoint
   */
  static async createWithAutoDetect(
    host: string = '127.0.0.1',
    portRange: { start: number; end: number } = { start: 5000, end: 5100 }
  ): Promise<PySCFApiClient> {
    for (let port = portRange.start; port <= portRange.end; port++) {
      try {
        const url = `http://${host}:${port}`;
        const client = new PySCFApiClient(url);
        await client.healthCheck();
        return client;
      } catch (error) {
        // Continue to next port
        continue;
      }
    }
    throw new Error(`No PySCF server found on ${host} in port range ${portRange.start}-${portRange.end}`);
  }

  /**
   * Health check
   */
  async healthCheck(): Promise<HealthResponse> {
    const response: AxiosResponse<HealthResponse> = await this.client.get('/health');
    return response.data;
  }

  // ====== PubChem and SMILES Methods ======

  async searchPubChem(request: PubChemSearchRequest): Promise<PubChemSearchResponse> {
    const response: AxiosResponse<PubChemSearchResponse> = await this.client.post('/api/pubchem/search', request);
    return response.data;
  }

  async convertSmiles(request: SMILESConvertRequest): Promise<SMILESConvertResponse> {
    const response: AxiosResponse<SMILESConvertResponse> = await this.client.post('/api/smiles/convert', request);
    return response.data;
  }

  async validateXYZ(request: XYZValidateRequest): Promise<XYZValidateResponse> {
    const response: AxiosResponse<XYZValidateResponse> = await this.client.post('/api/pubchem/validate', request);
    return response.data;
  }

  // ====== Quantum Calculation Methods ======

  async startCalculation(request: QuantumCalculationRequest): Promise<StartCalculationResponse> {
    try {
      // Use longer timeout for calculation requests
      const response: AxiosResponse<StartCalculationResponse> = await this.client.post('/api/quantum/calculate', request, {
        timeout: 60000, // 60 seconds
      });
      return response.data;
    } catch (error) {
      const errorDetails = this.extractErrorDetails(error, request);
      throw new PySCFApiError(errorDetails);
    }
  }

  async listCalculations(): Promise<CalculationListResponse> {
    const response: AxiosResponse<CalculationListResponse> = await this.client.get('/api/quantum/calculations');
    return response.data;
  }

  async getCalculationDetails(calculationId: string): Promise<CalculationDetailsResponse> {
    const response: AxiosResponse<CalculationDetailsResponse> = await this.client.get(
      `/api/quantum/calculations/${calculationId}`
    );
    return response.data;
  }

  async updateCalculation(calculationId: string, request: CalculationUpdateRequest): Promise<CalculationUpdateResponse> {
    const response: AxiosResponse<CalculationUpdateResponse> = await this.client.put(
      `/api/quantum/calculations/${calculationId}`,
      request
    );
    return response.data;
  }

  async deleteCalculation(calculationId: string): Promise<CalculationDeleteResponse> {
    const response: AxiosResponse<CalculationDeleteResponse> = await this.client.delete(
      `/api/quantum/calculations/${calculationId}`
    );
    return response.data;
  }

  // ====== Molecular Orbital Methods ======

  async getOrbitals(calculationId: string): Promise<OrbitalsResponse> {
    const response: AxiosResponse<OrbitalsResponse> = await this.client.get(
      `/api/quantum/calculations/${calculationId}/orbitals`
    );
    return response.data;
  }

  async getOrbitalCube(
    calculationId: string,
    orbitalIndex: number,
    options?: {
      gridSize?: number;
      isovaluePos?: number;
      isovalueNeg?: number;
    }
  ): Promise<OrbitalCubeResponse> {
    const params = new URLSearchParams();
    if (options?.gridSize) params.append('gridSize', options.gridSize.toString());
    if (options?.isovaluePos) params.append('isovaluePos', options.isovaluePos.toString());
    if (options?.isovalueNeg) params.append('isovalueNeg', options.isovalueNeg.toString());

    const url = `/api/quantum/calculations/${calculationId}/orbitals/${orbitalIndex}/cube`;
    const fullUrl = params.toString() ? `${url}?${params.toString()}` : url;

    const response: AxiosResponse<OrbitalCubeResponse> = await this.client.get(fullUrl, {
      timeout: 120000, // 2 minutes for CUBE generation
    });
    return response.data;
  }

  async listCubeFiles(calculationId: string): Promise<CubeFilesListResponse> {
    const response: AxiosResponse<CubeFilesListResponse> = await this.client.get(
      `/api/quantum/calculations/${calculationId}/orbitals/cube-files`
    );
    return response.data;
  }

  async deleteCubeFiles(calculationId: string, orbitalIndex?: number): Promise<CubeFilesDeleteResponse> {
    const params = new URLSearchParams();
    if (orbitalIndex !== undefined) params.append('orbital_index', orbitalIndex.toString());

    const url = `/api/quantum/calculations/${calculationId}/orbitals/cube-files`;
    const fullUrl = params.toString() ? `${url}?${params.toString()}` : url;

    const response: AxiosResponse<CubeFilesDeleteResponse> = await this.client.delete(fullUrl);
    return response.data;
  }

  // ====== IR Spectrum Methods ======

  async getIRSpectrum(
    calculationId: string,
    options?: {
      broadening_fwhm?: number;
      x_min?: number;
      x_max?: number;
      show_peaks?: boolean;
    }
  ): Promise<IRSpectrumResponse> {
    try {
      const params = new URLSearchParams();
      if (options?.broadening_fwhm) params.append('broadening_fwhm', options.broadening_fwhm.toString());
      if (options?.x_min) params.append('x_min', options.x_min.toString());
      if (options?.x_max) params.append('x_max', options.x_max.toString());
      if (options?.show_peaks !== undefined) params.append('show_peaks', options.show_peaks.toString());

      const url = `/api/quantum/calculations/${calculationId}/ir-spectrum`;
      const fullUrl = params.toString() ? `${url}?${params.toString()}` : url;

      const response: AxiosResponse<IRSpectrumResponse> = await this.client.get(fullUrl);
      return response.data;
    } catch (error) {
      const errorDetails = this.extractErrorDetails(error, { calculationId, options });
      throw new PySCFApiError(errorDetails);
    }
  }

  // ====== System Methods ======

  async getSupportedParameters(): Promise<SupportedParametersResponse> {
    const response: AxiosResponse<SupportedParametersResponse> = await this.client.get('/api/quantum/supported-parameters');
    return response.data;
  }

  async getSettings(): Promise<SettingsResponse> {
    const response: AxiosResponse<SettingsResponse> = await this.client.get('/api/settings');
    return response.data;
  }

  async updateSettings(request: SettingsUpdateRequest): Promise<SettingsResponse> {
    const response: AxiosResponse<SettingsResponse> = await this.client.put('/api/settings', request);
    return response.data;
  }

  async getResourceStatus(): Promise<SystemResourceResponse> {
    try {
      const response: AxiosResponse<SystemResourceResponse> = await this.client.get('/api/system/resource-status');
      return response.data;
    } catch (error) {
      const errorDetails = this.extractErrorDetails(error);
      throw new PySCFApiError(errorDetails);
    }
  }

  /**
   * Get the current base URL
   */
  getBaseUrl(): string {
    return this.baseUrl;
  }

  /**
   * Test connection to the API
   */
  async testConnection(): Promise<{ connected: boolean; version?: string; error?: string }> {
    try {
      const health = await this.healthCheck();
      return {
        connected: true,
        version: health.version,
      };
    } catch (error) {
      return {
        connected: false,
        error: error instanceof Error ? error.message : 'Unknown error',
      };
    }
  }
}

export default PySCFApiClient;