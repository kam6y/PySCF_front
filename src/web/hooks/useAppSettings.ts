import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { components } from '../types/generated-api';
import { getSettings, updateSettings } from '../apiClient';
import { useCalculationStore } from '../store/calculationStore';

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
    onSuccess: (updatedSettings, variables) => {
      // Check if calculations_directory changed BEFORE updating cache
      const oldSettings = queryClient.getQueryData<AppSettings>(
        settingsKeys.settings()
      );
      const directoryChanged =
        oldSettings &&
        oldSettings.calculations_directory !==
          updatedSettings.calculations_directory;

      // Update the cache with new settings
      queryClient.setQueryData(settingsKeys.settings(), updatedSettings);

      if (directoryChanged) {
        console.log(
          '[useAppSettings] Calculations directory changed, clearing active calculation and invalidating calculations list'
        );

        // Clear active calculation ID since old calculations may not exist in new directory
        useCalculationStore.getState().setActiveCalculationId(null);

        // Clear staged calculation as well
        useCalculationStore.getState().clearStagedCalculation();

        // Invalidate calculations list to force refetch from new directory
        queryClient.invalidateQueries({ queryKey: ['calculations'] });
      }

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
