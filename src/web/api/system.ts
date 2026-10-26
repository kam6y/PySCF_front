import type { components } from '../types/generated-api';
import { request } from './core';

type Gpu4PyscfStatusResponse = components['schemas']['Gpu4PyscfStatusResponse'];
type Gpu4PyscfInstallResponse =
  components['schemas']['Gpu4PyscfInstallResponse'];
type Gpu4PyscfInstallRequest = components['schemas']['Gpu4PyscfInstallRequest'];

export const getGpu4PyscfStatus = (): Promise<
  Gpu4PyscfStatusResponse['data']
> => {
  return request<Gpu4PyscfStatusResponse['data']>(
    '/api/system/gpu4pyscf-status',
    { method: 'GET' }
  );
};

export const installGpu4Pyscf = (
  payload?: Gpu4PyscfInstallRequest
): Promise<Gpu4PyscfInstallResponse['data']> => {
  const options: RequestInit = { method: 'POST' };
  if (payload && Object.keys(payload).length > 0) {
    options.body = JSON.stringify(payload);
  }

  return request<Gpu4PyscfInstallResponse['data']>(
    '/api/system/gpu4pyscf-install',
    options
  );
};
