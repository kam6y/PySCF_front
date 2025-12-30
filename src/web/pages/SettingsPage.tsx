import React, { useState, useEffect, useMemo } from 'react';
import { useMutation, useQuery } from '@tanstack/react-query';
import { useAppSettings } from '../hooks';
import { TIMEZONE_GROUPS, getTimezoneLabel } from '../utils/dateFormatter';
import type { components } from '../types/generated-api';
import {
  cancelGpuInstallJob,
  getGpuInstallJob,
  getGpuStatus,
  installGpu,
  listGpuInstallJobs,
} from '../apiClient';
import styles from './SettingsPage.module.css';

// Type for timezone from generated API types
type Timezone = components['schemas']['AppSettings']['timezone'];
type GpuStatus = components['schemas']['GpuStatus'];
type GpuInstallRequest = components['schemas']['GpuInstallRequest'];
type GpuInstallJob = components['schemas']['GpuInstallJob'];

// Default values constants
const DEFAULT_MAX_PARALLEL_INSTANCES = 4;
const DEFAULT_MAX_CPU_UTILIZATION = 95.0;
const DEFAULT_MAX_MEMORY_UTILIZATION = 95.0;

interface SettingsPageProps {
  // Props will be added when integrating with the main app
}

export const SettingsPage: React.FC<SettingsPageProps> = () => {
  const [maxParallelInstances, setMaxParallelInstances] = useState<
    number | undefined
  >(undefined);
  const [maxCpuUtilization, setMaxCpuUtilization] = useState<
    number | undefined
  >(undefined);
  const [maxMemoryUtilization, setMaxMemoryUtilization] = useState<
    number | undefined
  >(undefined);
  const [geminiApiKey, setGeminiApiKey] = useState<string>('');
  const [calculationsDirectory, setCalculationsDirectory] =
    useState<string>('');
  const [timezone, setTimezone] = useState<Timezone>('UTC');
  const [gpuEnabled, setGpuEnabled] = useState<boolean>(false);
  const [gpuPreferredPackage, setGpuPreferredPackage] = useState<string>('');
  const [activeInstallJobId, setActiveInstallJobId] = useState<string | null>(
    null
  );
  const [lastInstallJob, setLastInstallJob] = useState<GpuInstallJob | null>(
    null
  );
  const [isSelectingFolder, setIsSelectingFolder] = useState(false);
  const [originalValues, setOriginalValues] = useState<{
    maxParallelInstances?: number;
    maxCpuUtilization?: number;
    maxMemoryUtilization?: number;
    geminiApiKey?: string;
    calculationsDirectory?: string;
    timezone?: Timezone;
    gpuEnabled?: boolean;
    gpuPreferredPackage?: string;
  }>({});

  const {
    settings,
    isLoading,
    isUpdating,
    error,
    updateSettingsAsync,
    refetch: refetchSettings,
  } = useAppSettings();

  const isMacOS = useMemo(() => {
    if (typeof navigator === 'undefined') return false;
    return /Macintosh|Mac OS X/.test(navigator.userAgent);
  }, []);

  const {
    data: gpuStatusData,
    isFetching: isFetchingGpuStatus,
    refetch: refetchGpuStatus,
    error: gpuStatusError,
  } = useQuery<GpuStatus>({
    queryKey: ['gpu-status'],
    queryFn: () => getGpuStatus().then(response => response.status),
    refetchOnWindowFocus: false,
  });

  const {
    data: installJobData,
  } = useQuery<GpuInstallJob>({
    queryKey: ['gpu-install-job', activeInstallJobId],
    queryFn: () => getGpuInstallJob(activeInstallJobId!),
    enabled: !!activeInstallJobId,
    refetchOnWindowFocus: false,
    refetchInterval: activeInstallJobId ? 2000 : false,
  });

  const {
    data: installJobList,
    refetch: refetchInstallJobList,
    isFetching: isFetchingInstallJobList,
    error: installJobListError,
  } = useQuery<GpuInstallJob[]>({
    queryKey: ['gpu-install-jobs'],
    queryFn: () => listGpuInstallJobs(),
    refetchOnWindowFocus: false,
    refetchInterval: activeInstallJobId ? false : 10000,
  });

  const installGpuMutation = useMutation<GpuInstallJob, unknown, GpuInstallRequest>({
    mutationFn: (payload: GpuInstallRequest) => installGpu(payload),
    onSuccess: async job => {
      setActiveInstallJobId(job.job_id);
      setLastInstallJob(job);
      await refetchInstallJobList();
      await refetchGpuStatus();
    },
    onError: error => {
      console.error('Failed to install gpu4pyscf:', error);
    },
  });

  const cancelInstallJobMutation = useMutation<GpuInstallJob, unknown, string>({
    mutationFn: (jobId: string) => cancelGpuInstallJob(jobId),
    onSuccess: async job => {
      setLastInstallJob(job);
      setActiveInstallJobId(null);
      await refetchGpuStatus();
      await refetchSettings();
      await refetchInstallJobList();
    },
    onError: error => {
      console.error('Failed to cancel gpu4pyscf install job:', error);
    },
  });

  // Update local state when settings are loaded
  useEffect(() => {
    if (settings) {
      // Explicit type conversion and validation
      const maxParallelInstancesValue = Number(settings.max_parallel_instances);
      const maxCpuUtilizationValue = Number(
        settings.max_cpu_utilization_percent
      );
      const maxMemoryUtilizationValue = Number(
        settings.max_memory_utilization_percent
      );

      const newValues = {
        maxParallelInstances: !isNaN(maxParallelInstancesValue)
          ? maxParallelInstancesValue
          : DEFAULT_MAX_PARALLEL_INSTANCES,
        maxCpuUtilization: !isNaN(maxCpuUtilizationValue)
          ? maxCpuUtilizationValue
          : DEFAULT_MAX_CPU_UTILIZATION,
        maxMemoryUtilization: !isNaN(maxMemoryUtilizationValue)
          ? maxMemoryUtilizationValue
          : DEFAULT_MAX_MEMORY_UTILIZATION,
        geminiApiKey: settings.gemini_api_key || '',
        calculationsDirectory: settings.calculations_directory || '',
        timezone: settings.timezone || 'UTC',
        gpuEnabled: settings.gpu_acceleration_enabled ?? false,
        gpuPreferredPackage: settings.gpu_preferred_package || '',
      };

      if (process.env.NODE_ENV === 'development') {
        console.log('SettingsPage: Loading settings', { settings, newValues });
      }

      setMaxParallelInstances(newValues.maxParallelInstances);
      setMaxCpuUtilization(newValues.maxCpuUtilization);
      setMaxMemoryUtilization(newValues.maxMemoryUtilization);
      setGeminiApiKey(newValues.geminiApiKey);
      setCalculationsDirectory(newValues.calculationsDirectory);
      setTimezone(newValues.timezone);
      setGpuEnabled(newValues.gpuEnabled);
      setGpuPreferredPackage(newValues.gpuPreferredPackage);
      setOriginalValues(newValues);

      if (process.env.NODE_ENV === 'development') {
        console.log('SettingsPage: State updated', {
          originalValues: newValues,
        });
      }
    }
  }, [settings]);

  useEffect(() => {
    if (installJobData) {
      setLastInstallJob(installJobData);
      const status = installJobData.status;
      const isTerminal =
        status === 'succeeded' ||
        status === 'failed' ||
        status === 'canceled';

      if (isTerminal) {
        setActiveInstallJobId(null);
        const result = installJobData.result_status;
        if (result?.status === 'ready') {
          setGpuEnabled(true);
          setGpuPreferredPackage(
            result.installed_package || result.recommended_package || ''
          );
        }
        refetchGpuStatus();
        refetchSettings();
        refetchInstallJobList();
      }
    }
  }, [installJobData, refetchGpuStatus, refetchSettings, refetchInstallJobList]);

  useEffect(() => {
    if (!installJobList || installJobList.length === 0) {
      setLastInstallJob(null);
      if (activeInstallJobId) {
        setActiveInstallJobId(null);
      }
      return;
    }

    const inProgress = installJobList.find(
      job => job.status === 'queued' || job.status === 'running'
    );
    const newest = installJobList[0];

    setLastInstallJob(inProgress || newest);

    const activeJobExists = activeInstallJobId
      ? installJobList.some(job => job.job_id === activeInstallJobId)
      : false;

    if (!activeJobExists && inProgress) {
      setActiveInstallJobId(inProgress.job_id);
    }
  }, [installJobList, activeInstallJobId]);

  const handleSave = async () => {
    try {
      const payload = {
        max_parallel_instances:
          maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES,
        max_cpu_utilization_percent:
          maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION,
        max_memory_utilization_percent:
          maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION,
        system_total_cores: settings?.system_total_cores || 0,
        system_total_memory_mb: settings?.system_total_memory_mb || 0,
        calculations_directory: calculationsDirectory,
        timezone: timezone,
        gemini_api_key: geminiApiKey || null,
        gpu_acceleration_enabled: gpuEnabled,
        gpu_preferred_package: gpuPreferredPackage || null,
      };
      await updateSettingsAsync(payload);

      const newValues = {
        maxParallelInstances,
        maxCpuUtilization,
        maxMemoryUtilization,
        geminiApiKey,
        calculationsDirectory,
        timezone,
        gpuEnabled,
        gpuPreferredPackage,
      };
      setOriginalValues(newValues);
    } catch (error) {
      console.error('Failed to save settings:', error);
      // Reset to original values on error
      setMaxParallelInstances(originalValues.maxParallelInstances);
      setMaxCpuUtilization(originalValues.maxCpuUtilization);
      setMaxMemoryUtilization(originalValues.maxMemoryUtilization);
      setGeminiApiKey(originalValues.geminiApiKey || '');
      setCalculationsDirectory(originalValues.calculationsDirectory || '');
      setTimezone(originalValues.timezone || 'UTC');
      setGpuEnabled(originalValues.gpuEnabled ?? false);
      setGpuPreferredPackage(originalValues.gpuPreferredPackage || '');
    }
  };

  const handleCancel = () => {
    setMaxParallelInstances(originalValues.maxParallelInstances);
    setMaxCpuUtilization(originalValues.maxCpuUtilization);
    setMaxMemoryUtilization(originalValues.maxMemoryUtilization);
    setGeminiApiKey(originalValues.geminiApiKey || '');
    setCalculationsDirectory(originalValues.calculationsDirectory || '');
    setTimezone(originalValues.timezone || 'UTC');
    setGpuEnabled(originalValues.gpuEnabled ?? false);
    setGpuPreferredPackage(originalValues.gpuPreferredPackage || '');
  };

  const handleSelectFolder = async () => {
    setIsSelectingFolder(true);
    try {
      const result = await window.electronAPI.selectFolder();
      if (!result.canceled && result.filePath) {
        // Append /PySCF_calculations to the selected path
        const fullPath = `${result.filePath}/PySCF_calculations`;
        setCalculationsDirectory(fullPath);
      }
    } catch (error) {
      console.error('Failed to select folder:', error);
    } finally {
      setIsSelectingFolder(false);
    }
  };

  const handleInstallGpu = () => {
    const packageName =
      gpuPreferredPackage ||
      gpuStatusData?.recommended_package ||
      'gpu4pyscf';

    installGpuMutation.mutate({
      package: packageName,
      enable_gpu: true,
    });
  };

  const handleCancelInstall = () => {
    if (!activeInstallJobId) return;
    cancelInstallJobMutation.mutate(activeInstallJobId);
  };

  const gpuStatus = gpuStatusData;
  const gpuStatusReady = gpuStatus?.status === 'ready';
  const canInstallGpu =
    !!gpuStatus &&
    gpuStatus.gpu_supported_platform &&
    gpuStatus.has_nvidia_gpu;
  const currentInstallJob = installJobData || lastInstallJob;
  const installJobStatus = currentInstallJob?.status;
  const installJobInProgress =
    installJobStatus === 'queued' || installJobStatus === 'running';
  const gpuStatusLabel = useMemo(() => {
    if (isFetchingGpuStatus) {
      return 'Checking status...';
    }
    if (!gpuStatus) {
      return 'Not fetched';
    }
    switch (gpuStatus.status) {
      case 'ready':
        return 'Ready';
      case 'not_installed':
        return 'gpu4pyscf not installed';
      case 'missing_cuda':
        return 'CUDA not detected';
      case 'missing_gpu':
        return 'GPU not detected';
      case 'unsupported_platform':
        return 'Linux only';
      default:
        return 'GPU error';
    }
  }, [gpuStatus, isFetchingGpuStatus]);

  const gpuBadgeClass = useMemo(() => {
    const base = styles.statusBadge;
    if (!gpuStatus || isFetchingGpuStatus) {
      return base;
    }
    if (gpuStatus.status === 'ready') {
      return `${base} ${styles.statusReady}`;
    }
    if (
      gpuStatus.status === 'not_installed' ||
      gpuStatus.status === 'missing_cuda' ||
      gpuStatus.status === 'missing_gpu' ||
      gpuStatus.status === 'unsupported_platform'
    ) {
      return `${base} ${styles.statusWarning}`;
    }
    return `${base} ${styles.statusError}`;
  }, [gpuStatus, isFetchingGpuStatus]);

  const installJobStatusLabel = useMemo(() => {
    if (!currentInstallJob) return null;
    switch (currentInstallJob.status) {
      case 'queued':
        return 'Queued';
      case 'running':
        return 'Installing';
      case 'succeeded':
        return 'Completed';
      case 'failed':
        return 'Failed';
      case 'canceled':
        return 'Canceled';
      default:
        return currentInstallJob.status;
    }
  }, [currentInstallJob]);

  const recommendedPackage =
    gpuStatus?.recommended_package ||
    currentInstallJob?.result_status?.recommended_package ||
    gpuPreferredPackage ||
    'gpu4pyscf';
  const detectedGpuNames =
    gpuStatus?.detected_gpus
      ?.map(gpu => gpu?.name)
      .filter(Boolean)
      .join(', ') || 'Not detected';
  const gpuMessage =
    gpuStatus?.message ||
    'On Linux with an NVIDIA GPU, install the gpu4pyscf package that matches your CUDA version to use it.';
  const gpuErrorText = useMemo(() => {
    if (gpuStatusError) {
      if (gpuStatusError instanceof Error) return gpuStatusError.message;
      if (typeof gpuStatusError === 'string') return gpuStatusError;
    }
    if (installGpuMutation.error) {
      if (installGpuMutation.error instanceof Error)
        return installGpuMutation.error.message;
      if (typeof installGpuMutation.error === 'string')
        return installGpuMutation.error;
    }
    if (installJobListError) {
      if (installJobListError instanceof Error) return installJobListError.message;
      if (typeof installJobListError === 'string') return installJobListError;
    }
    if (cancelInstallJobMutation.error) {
      if (cancelInstallJobMutation.error instanceof Error)
        return cancelInstallJobMutation.error.message;
      if (typeof cancelInstallJobMutation.error === 'string')
        return cancelInstallJobMutation.error;
    }
    if (currentInstallJob?.error) {
      return currentInstallJob.error;
    }
    return null;
  }, [
    gpuStatusError,
    installGpuMutation.error,
    installJobListError,
    cancelInstallJobMutation.error,
    currentInstallJob?.error,
  ]);

  const hasUnsavedChanges = useMemo(() => {
    // Return false if settings haven't been loaded yet
    if (!originalValues || Object.keys(originalValues).length === 0) {
      return false;
    }

    // Use default values for comparison if values are undefined
    const currentParallel =
      maxParallelInstances ?? DEFAULT_MAX_PARALLEL_INSTANCES;
    const currentCpu = maxCpuUtilization ?? DEFAULT_MAX_CPU_UTILIZATION;
    const currentMemory =
      maxMemoryUtilization ?? DEFAULT_MAX_MEMORY_UTILIZATION;
    const currentGeminiApiKey = geminiApiKey;
    const currentCalculationsDirectory = calculationsDirectory;
    const currentTimezone = timezone;
    const currentGpuEnabled = gpuEnabled;
    const currentGpuPreferred = gpuPreferredPackage;

    const originalParallel =
      originalValues.maxParallelInstances ?? DEFAULT_MAX_PARALLEL_INSTANCES;
    const originalCpu =
      originalValues.maxCpuUtilization ?? DEFAULT_MAX_CPU_UTILIZATION;
    const originalMemory =
      originalValues.maxMemoryUtilization ?? DEFAULT_MAX_MEMORY_UTILIZATION;
    const originalGeminiApiKey = originalValues.geminiApiKey || '';
    const originalCalculationsDirectory =
      originalValues.calculationsDirectory || '';
    const originalTimezone = originalValues.timezone || 'UTC';
    const originalGpuEnabled = originalValues.gpuEnabled ?? false;
    const originalGpuPreferred = originalValues.gpuPreferredPackage || '';

    const parallelChanged = currentParallel !== originalParallel;
    const cpuChanged = Math.abs(currentCpu - originalCpu) > 0.001;
    const memoryChanged = Math.abs(currentMemory - originalMemory) > 0.001;
    const geminiApiKeyChanged = currentGeminiApiKey !== originalGeminiApiKey;
    const calculationsDirectoryChanged =
      currentCalculationsDirectory !== originalCalculationsDirectory;
    const timezoneChanged = currentTimezone !== originalTimezone;
    const gpuEnabledChanged = currentGpuEnabled !== originalGpuEnabled;
    const gpuPreferredChanged = currentGpuPreferred !== originalGpuPreferred;

    const hasChanges =
      parallelChanged ||
      cpuChanged ||
      memoryChanged ||
      geminiApiKeyChanged ||
      calculationsDirectoryChanged ||
      timezoneChanged ||
      gpuEnabledChanged ||
      gpuPreferredChanged;

    // Debug logging in development
    if (process.env.NODE_ENV === 'development') {
      console.log('SettingsPage: hasUnsavedChanges check', {
        current: {
          currentParallel,
          currentCpu,
          currentMemory,
          currentGeminiApiKey,
        },
        original: {
          originalParallel,
          originalCpu,
          originalMemory,
          originalGeminiApiKey,
        },
        changes: {
          parallelChanged,
          cpuChanged,
          memoryChanged,
          geminiApiKeyChanged,
        },
        hasChanges,
      });
    }

    return hasChanges;
  }, [
    maxParallelInstances,
    maxCpuUtilization,
    maxMemoryUtilization,
    geminiApiKey,
    calculationsDirectory,
    timezone,
    gpuEnabled,
    gpuPreferredPackage,
    originalValues?.maxParallelInstances,
    originalValues?.maxCpuUtilization,
    originalValues?.maxMemoryUtilization,
    originalValues?.geminiApiKey,
    originalValues?.calculationsDirectory,
    originalValues?.timezone,
    originalValues?.gpuEnabled,
    originalValues?.gpuPreferredPackage,
  ]);

  if (isLoading) {
    return (
      <div className={styles.settingsPage}>
        <div className={styles.settingsHeader}>
          <h2>Settings</h2>
        </div>
        <div className={styles.loadingContainer}>
          <div className={styles.loadingSpinner}></div>
          <p>Loading settings...</p>
        </div>
      </div>
    );
  }

  return (
    <div className={styles.settingsPage}>
      <div className={styles.settingsHeader}>
        <div className={styles.settingsHeaderLeft}>
          <h2>Settings</h2>
          <p className={styles.settingsDescription}>
            Configure application behavior and performance settings
          </p>
        </div>
        {hasUnsavedChanges && (
          <div className={styles.settingsHeaderActions}>
            <button
              onClick={handleCancel}
              className={styles.cancelButton}
              disabled={isUpdating}
            >
              Cancel
            </button>
            <button
              onClick={handleSave}
              className={styles.saveButton}
              disabled={isUpdating}
            >
              {isUpdating ? (
                <>
                  <div className={styles.buttonSpinner}></div>
                  Saving...
                </>
              ) : (
                'Save Changes'
              )}
            </button>
          </div>
        )}
      </div>

      <div className={styles.settingsContent}>
        {/* Error Display */}
        {error && (
          <div className={styles.errorContainer}>
            <div className={styles.errorIcon}>⚠️</div>
            <div className={styles.errorContent}>
              <h4>Error Loading Settings</h4>
              <p>
                {error instanceof Error
                  ? error.message
                  : 'An unexpected error occurred'}
              </p>
            </div>
          </div>
        )}

        <div className={styles.settingsSection}>
          <h3>Performance</h3>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label htmlFor="maxParallelInstances">
                Maximum Parallel Calculations
              </label>
              <p className={styles.settingHelp}>
                Maximum number of calculations that can run simultaneously.
                Additional calculations will wait in queue.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.sliderContainer}>
                <input
                  id="maxParallelInstances"
                  type="range"
                  min="1"
                  max="16"
                  step="1"
                  value={maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES}
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    if (process.env.NODE_ENV === 'development') {
                      console.log(
                        'SettingsPage: maxParallelInstances changed',
                        newValue
                      );
                    }
                    setMaxParallelInstances(newValue);
                  }}
                  className={styles.slider}
                  disabled={isUpdating}
                />
                <div className={styles.sliderLabels}>
                  <span>1</span>
                  <span>4</span>
                  <span>8</span>
                  <span>12</span>
                  <span>16</span>
                </div>
              </div>

              <div className={styles.valueDisplay}>
                <span className={styles.currentValue}>
                  {maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES}
                </span>
                <span className={styles.valueUnit}>calculations</span>
              </div>
            </div>
          </div>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label htmlFor="maxCpuUtilization">Maximum CPU Utilization</label>
              <p className={styles.settingHelp}>
                Maximum percentage of system CPU that can be used for
                calculations. Helps prevent system overload and maintains
                responsiveness.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.sliderContainer}>
                <input
                  type="range"
                  id="maxCpuUtilization"
                  min="10"
                  max="100"
                  step="5"
                  value={maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION}
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    if (process.env.NODE_ENV === 'development') {
                      console.log(
                        'SettingsPage: maxCpuUtilization changed',
                        newValue
                      );
                    }
                    setMaxCpuUtilization(newValue);
                  }}
                  className={styles.slider}
                  disabled={isUpdating}
                />
                <div className={styles.sliderLabels}>
                  <span>10%</span>
                  <span>50%</span>
                  <span>100%</span>
                </div>
              </div>

              <div className={styles.valueDisplay}>
                <span className={styles.currentValue}>
                  {(maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION).toFixed(
                    0
                  )}
                  %
                </span>
                <span className={styles.valueUnit}>CPU</span>
              </div>
            </div>
          </div>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label htmlFor="maxMemoryUtilization">
                Maximum Memory Utilization
              </label>
              <p className={styles.settingHelp}>
                Maximum percentage of system memory that can be used for
                calculations. Prevents out-of-memory issues and maintains system
                stability.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.sliderContainer}>
                <input
                  type="range"
                  id="maxMemoryUtilization"
                  min="10"
                  max="100"
                  step="5"
                  value={maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION}
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    if (process.env.NODE_ENV === 'development') {
                      console.log(
                        'SettingsPage: maxMemoryUtilization changed',
                        newValue
                      );
                    }
                    setMaxMemoryUtilization(newValue);
                  }}
                  className={styles.slider}
                  disabled={isUpdating}
                />
                <div className={styles.sliderLabels}>
                  <span>10%</span>
                  <span>50%</span>
                  <span>100%</span>
                </div>
              </div>

              <div className={styles.valueDisplay}>
                <span className={styles.currentValue}>
                  {(
                    maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION
                  ).toFixed(0)}
                  %
                </span>
                <span className={styles.valueUnit}>Memory</span>
              </div>
            </div>
          </div>
          {settings && (
            <div className={styles.systemInfoSection}>
              <h4>System Information</h4>
              <div className={styles.systemInfoGrid}>
                <div className={styles.systemInfoItem}>
                  <span className={styles.infoLabel}>CPU Cores:</span>
                  <span className={styles.infoValue}>
                    {settings.system_total_cores}
                  </span>
                </div>
                <div className={styles.systemInfoItem}>
                  <span className={styles.infoLabel}>Total Memory:</span>
                  <span className={styles.infoValue}>
                    {(settings.system_total_memory_mb / 1024).toFixed(1)} GB
                  </span>
                </div>
              </div>
            </div>
          )}
        </div>

        <div className={styles.settingsSection}>
          <h3>GPU Acceleration (Linux + NVIDIA)</h3>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label>GPU computations with gpu4pyscf</label>
              <p className={styles.settingHelp}>
                On Linux, install the gpu4pyscf package matching your CUDA version to accelerate HF/DFT/TDDFT calculations on NVIDIA GPUs.
              </p>
            </div>

            <div className={styles.settingControl}>
              {isMacOS ? (
                <>
                  <div className={styles.gpuStatusRow}>
                    <span className={`${styles.statusBadge} ${styles.statusWarning}`}>
                      Linux only
                    </span>
                  </div>
                  <p className={styles.settingHelp}>
                    GPU acceleration is supported only on Linux.
                  </p>
                </>
              ) : (
                <>
                  <div className={styles.gpuStatusRow}>
                    <span className={gpuBadgeClass}>{gpuStatusLabel}</span>
                    <div className={styles.gpuButtons}>
                      <button
                        onClick={() => {
                          refetchGpuStatus();
                          refetchInstallJobList();
                        }}
                        className={styles.refreshButton}
                        disabled={isFetchingGpuStatus || isFetchingInstallJobList}
                      >
                        {isFetchingGpuStatus || isFetchingInstallJobList
                          ? 'Rescanning...'
                          : 'Rescan status'}
                      </button>
                      <button
                        onClick={handleInstallGpu}
                        className={styles.installGpuButton}
                        disabled={
                          installGpuMutation.isPending ||
                          !canInstallGpu ||
                          isFetchingGpuStatus ||
                          installJobInProgress
                        }
                      >
                        {installGpuMutation.isPending
                          ? 'Submitting job...'
                          : installJobInProgress
                          ? 'Installation in progress'
                          : 'Install recommended gpu4pyscf and enable'}
                      </button>
                      {installJobInProgress && (
                        <button
                          onClick={handleCancelInstall}
                          className={styles.cancelInstallButton}
                          disabled={cancelInstallJobMutation.isPending}
                        >
                          {cancelInstallJobMutation.isPending
                            ? 'Canceling...'
                            : 'Cancel installation'}
                        </button>
                      )}
                    </div>
                  </div>

                  {currentInstallJob && (
                    <p className={styles.settingHelp}>
                      Install job {currentInstallJob.job_id} :{' '}
                      {installJobStatusLabel || 'Fetching status'}
                    </p>
                  )}

                  <p className={styles.settingHelp}>{gpuMessage}</p>

                  {gpuErrorText && (
                    <div className={styles.warningBox}>
                      <span className={styles.warningIcon}>⚠</span>
                      <span className={styles.warningText}>{gpuErrorText}</span>
                    </div>
                  )}

                  <div className={styles.gpuDetailGrid}>
                    <div className={styles.gpuDetailItem}>
                      <div className={styles.infoLabel}>CUDA version</div>
                      <div className={styles.infoValue}>
                        {gpuStatus?.cuda_version || 'Not detected'}
                      </div>
                    </div>
                    <div className={styles.gpuDetailItem}>
                      <div className={styles.infoLabel}>Driver</div>
                      <div className={styles.infoValue}>
                        {gpuStatus?.driver_version || 'Not detected'}
                      </div>
                    </div>
                    <div className={styles.gpuDetailItem}>
                      <div className={styles.infoLabel}>GPU</div>
                      <div className={styles.infoValue}>{detectedGpuNames}</div>
                    </div>
                    <div className={styles.gpuDetailItem}>
                      <div className={styles.infoLabel}>Recommended package</div>
                      <div className={styles.infoValue}>{recommendedPackage}</div>
                    </div>
                    <div className={styles.gpuDetailItem}>
                      <div className={styles.infoLabel}>Installed</div>
                      <div className={styles.infoValue}>
                        {gpuStatus?.installed_package || 'Not installed'}
                      </div>
                    </div>
                    <div className={styles.gpuDetailItem}>
                      <div className={styles.infoLabel}>Configured package</div>
                      <div className={styles.infoValue}>
                        {gpuPreferredPackage || 'Not configured'}
                      </div>
                    </div>
                  </div>

                  <div className={styles.checkboxRow}>
                    <input
                      id="gpuEnabled"
                      type="checkbox"
                      checked={gpuEnabled}
                      disabled={!gpuStatusReady && !gpuEnabled}
                      onChange={e => setGpuEnabled(e.target.checked)}
                    />
                    <label htmlFor="gpuEnabled">
                      Use GPU acceleration for HF/DFT/TDDFT (Linux only)
                    </label>
                  </div>
                  <p className={styles.settingHelp}>
                    {gpuStatusReady
                      ? 'Saving will run HF/DFT/TDDFT calculations through gpu4pyscf.'
                      : 'Available when Linux + NVIDIA GPU + CUDA are detected.'}
                  </p>
                </>
              )}
            </div>
          </div>
        </div>

        <div className={styles.settingsSection}>
          <h3>AI Agent</h3>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label htmlFor="geminiApiKey">Google Gemini API Key</label>
              <p className={styles.settingHelp}>
                API key for Google Gemini AI to enable intelligent molecular
                analysis and assistance. Leave empty to use fallback responses
                without AI features.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.textInputContainer}>
                <input
                  id="geminiApiKey"
                  type="password"
                  placeholder="Enter your Gemini API key..."
                  value={geminiApiKey}
                  onChange={e => {
                    const newValue = e.target.value;
                    if (process.env.NODE_ENV === 'development') {
                      console.log(
                        'SettingsPage: geminiApiKey changed (length)',
                        newValue.length
                      );
                    }
                    setGeminiApiKey(newValue);
                  }}
                  className={styles.textInput}
                  disabled={isUpdating}
                />
                <div className={styles.inputStatus}>
                  {geminiApiKey ? (
                    <span className={styles.statusConfigured}>
                      ✓ API Key Configured
                    </span>
                  ) : (
                    <span className={styles.statusNotConfigured}>
                      ⚠ API Key Not Set
                    </span>
                  )}
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className={styles.settingsSection}>
          <h3>Display</h3>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label htmlFor="timezone">Timezone for Date & Time Display</label>
              <p className={styles.settingHelp}>
                Timezone used for displaying timestamps in chat history and
                calculation history. Changes apply immediately to all displayed
                dates.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.selectContainer}>
                <select
                  id="timezone"
                  value={timezone}
                  onChange={e => {
                    const newValue = e.target.value as Timezone;
                    if (process.env.NODE_ENV === 'development') {
                      console.log('SettingsPage: timezone changed', newValue);
                    }
                    setTimezone(newValue);
                  }}
                  className={styles.selectInput}
                  disabled={isUpdating}
                >
                  <optgroup label="Standard">
                    {TIMEZONE_GROUPS.standard.map(tz => (
                      <option key={tz} value={tz}>
                        {getTimezoneLabel(tz)}
                      </option>
                    ))}
                  </optgroup>
                  <optgroup label="Asia / Pacific">
                    {TIMEZONE_GROUPS.asiaPacific.map(tz => (
                      <option key={tz} value={tz}>
                        {getTimezoneLabel(tz)}
                      </option>
                    ))}
                  </optgroup>
                  <optgroup label="Europe">
                    {TIMEZONE_GROUPS.europe.map(tz => (
                      <option key={tz} value={tz}>
                        {getTimezoneLabel(tz)}
                      </option>
                    ))}
                  </optgroup>
                  <optgroup label="Americas">
                    {TIMEZONE_GROUPS.americas.map(tz => (
                      <option key={tz} value={tz}>
                        {getTimezoneLabel(tz)}
                      </option>
                    ))}
                  </optgroup>
                </select>
                <div className={styles.selectedTimezone}>
                  <span className={styles.timezoneLabel}>Selected:</span>
                  <span className={styles.timezoneValue}>
                    {getTimezoneLabel(timezone)}
                  </span>
                </div>
              </div>
            </div>
          </div>
        </div>

        <div className={styles.settingsSection}>
          <h3>Storage</h3>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label htmlFor="calculationsDirectory">
                Calculation Data Folder
              </label>
              <p className={styles.settingHelp}>
                Location where all calculation data is stored. A
                'PySCF_calculations' subfolder will be created in the selected
                directory. Existing data will be moved when changed.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.textInputContainer}>
                <input
                  id="calculationsDirectory"
                  type="text"
                  placeholder="Calculations directory path..."
                  value={calculationsDirectory}
                  readOnly
                  className={styles.textInput}
                  style={{
                    fontFamily: 'Monaco, Menlo, Courier New, monospace',
                  }}
                />
                <button
                  onClick={handleSelectFolder}
                  className={styles.selectButton}
                  disabled={isUpdating || isSelectingFolder}
                >
                  {isSelectingFolder ? 'Selecting...' : 'Change Folder...'}
                </button>
              </div>
              {calculationsDirectory !==
                originalValues.calculationsDirectory && (
                <div className={styles.warningBox}>
                  <span className={styles.warningIcon}>⚠</span>
                  <span className={styles.warningText}>
                    Changing this folder will move all existing calculation data
                    to the new location. This operation may take a while.
                  </span>
                </div>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
