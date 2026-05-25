import {
  CalculationDetailsResponseData,
  CalculationListResponseData,
  QuantumCalculationRequest,
  CalculationUpdateResponseData,
  CalculationDeleteResponseData,
  StartCalculationResponseData,
  OrbitalsResponseData,
  OrbitalCubeResponseData,
  SupportedParametersResponseData,
  IRSpectrumResponseData,
  PauseCalculationResponseData,
  ResumeCalculationResponseData,
} from '../types/api-types';
import { request, ApiError } from './core';

export type StartCalculationResponse = StartCalculationResponseData;

const validateCalculationId = (id: string | null, endpoint: string) => {
  if (!id || id === 'undefined' || id === 'null') {
    throw new ApiError(
      'Invalid calculation ID provided.',
      400,
      'Bad Request',
      endpoint
    );
  }
};

export const getCalculations = (): Promise<CalculationListResponseData> => {
  return request<CalculationListResponseData>('/api/quantum/calculations', {
    method: 'GET',
  });
};

export const getCalculationDetails = (
  id: string
): Promise<CalculationDetailsResponseData> => {
  validateCalculationId(id, `/api/quantum/calculations/${id}`);
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

export const pauseCalculation = (
  id: string
): Promise<PauseCalculationResponseData> => {
  validateCalculationId(id, `/api/quantum/calculations/${id}/pause`);
  return request<PauseCalculationResponseData>(
    `/api/quantum/calculations/${id}/pause`,
    {
      method: 'POST',
    }
  );
};

export const resumeCalculation = (
  id: string
): Promise<ResumeCalculationResponseData> => {
  validateCalculationId(id, `/api/quantum/calculations/${id}/resume`);
  return request<ResumeCalculationResponseData>(
    `/api/quantum/calculations/${id}/resume`,
    {
      method: 'POST',
    }
  );
};

export const getOrbitals = (
  calculationId: string
): Promise<OrbitalsResponseData> => {
  validateCalculationId(
    calculationId,
    `/api/quantum/calculations/${calculationId}/orbitals`
  );
  return request<OrbitalsResponseData>(
    `/api/quantum/calculations/${calculationId}/orbitals`,
    { method: 'GET' }
  );
};

export const getOrbitalCube = (
  calculationId: string,
  orbitalIndex: number,
  options?: {
    gridSize?: number;
    isovaluePos?: number;
    isovalueNeg?: number;
  }
): Promise<OrbitalCubeResponseData> => {
  validateCalculationId(
    calculationId,
    `/api/quantum/calculations/${calculationId}/orbitals/${orbitalIndex}/cube`
  );

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

export const getSupportedParameters =
  (): Promise<SupportedParametersResponseData> => {
    return request<SupportedParametersResponseData>(
      '/api/quantum/supported-parameters',
      { method: 'GET' }
    );
  };

export const getIRSpectrum = (
  calculationId: string,
  options?: {
    broadening_fwhm?: number;
    x_min?: number;
    x_max?: number;
    show_peaks?: boolean;
  }
): Promise<IRSpectrumResponseData> => {
  validateCalculationId(
    calculationId,
    `/api/quantum/calculations/${calculationId}/ir-spectrum`
  );

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
