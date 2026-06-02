import { components } from '../types/generated-api';
import { request } from './core';

type AppSettings = components['schemas']['AppSettings'];
type AppSettingsResponse = components['schemas']['AppSettingsResponse'];
type SettingsResponse = components['schemas']['SettingsResponse'];

export type { AppSettings, AppSettingsResponse };

export const getSettings = (): Promise<SettingsResponse['data']> => {
  return request<SettingsResponse['data']>('/api/settings', { method: 'GET' });
};

export const updateSettings = (
  settings: AppSettings
): Promise<SettingsResponse['data']> => {
  return request<SettingsResponse['data']>('/api/settings', {
    method: 'PUT',
    body: JSON.stringify(settings),
  });
};
