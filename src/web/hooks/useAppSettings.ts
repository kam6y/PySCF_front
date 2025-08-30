import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { components } from '../types/generated-api';
import { getSettings, updateSettings } from '../apiClient';

// Type definitions
type AppSettings = components['schemas']['AppSettings'];
type SettingsResponse = components['schemas']['SettingsResponse'];

// API Client functions are now imported from ../apiClient

// Query keys
const settingsKeys = {
  all: ['settings'] as const,
  settings: () => [...settingsKeys.all, 'current'] as const,
};

// Custom hooks
export const useGetSettings = () => {
  return useQuery({
    queryKey: settingsKeys.settings(),
    queryFn: () => getSettings().then(response => response.settings),
    staleTime: 5 * 60 * 1000, // 5 minutes
    retry: 2,
  });
};

export const useUpdateSettings = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: (settings: AppSettings) =>
      updateSettings(settings).then(response => response.settings),
    onSuccess: updatedSettings => {
      // Update the cache with new settings
      queryClient.setQueryData(settingsKeys.settings(), updatedSettings);

      // Invalidate related queries to ensure consistency
      queryClient.invalidateQueries({ queryKey: settingsKeys.all });
    },
    onError: error => {
      console.error('Failed to update settings:', error);
    },
  });
};

// Main hook that combines get and update functionality
export const useAppSettings = () => {
  const getSettingsQuery = useGetSettings();
  const updateSettingsMutation = useUpdateSettings();

  return {
    // Settings data
    settings: getSettingsQuery.data,

    // Loading states
    isLoading: getSettingsQuery.isLoading,
    isUpdating: updateSettingsMutation.isPending,

    // Error states
    error: getSettingsQuery.error || updateSettingsMutation.error,

    // Update function
    updateSettings: updateSettingsMutation.mutate,
    updateSettingsAsync: updateSettingsMutation.mutateAsync,

    // Status flags
    isSuccess: getSettingsQuery.isSuccess && !updateSettingsMutation.isPending,
    isError: getSettingsQuery.isError || updateSettingsMutation.isError,

    // Refetch function
    refetch: getSettingsQuery.refetch,
  };
};
