import React, { useState, useEffect, useMemo, useRef } from 'react';
import { useAppSettings, useGpu4Pyscf } from '../hooks';
import {
  TIMEZONE_LABELS,
  TIMEZONE_GROUPS,
  getTimezoneLabel,
} from '../utils/dateFormatter';
import type { components } from '../types/generated-api';
import styles from './SettingsPage.module.css';

// Type for timezone from generated API types
type Timezone = components['schemas']['AppSettings']['timezone'];

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
  const [isSelectingFolder, setIsSelectingFolder] = useState(false);
  const [gpuAccelerationEnabled, setGpuAccelerationEnabled] =
    useState<boolean>(false);
  const [originalValues, setOriginalValues] = useState<{
    maxParallelInstances?: number;
    maxCpuUtilization?: number;
    maxMemoryUtilization?: number;
    geminiApiKey?: string;
    calculationsDirectory?: string;
    timezone?: Timezone;
    gpuAccelerationEnabled?: boolean;
  }>({});

  const { settings, isLoading, isUpdating, error, updateSettings } =
    useAppSettings();
  const {
    status: gpuStatus,
    isLoading: isGpuStatusLoading,
    isFetching: isGpuStatusFetching,
    isInstalling: isGpuInstalling,
    error: gpuError,
    installGpu4Pyscf,
    refetch: refetchGpuStatus,
  } = useGpu4Pyscf();

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
        gpuAccelerationEnabled: settings.gpu_acceleration_enabled ?? false,
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
      setGpuAccelerationEnabled(newValues.gpuAccelerationEnabled ?? false);
      setOriginalValues(newValues);

      if (process.env.NODE_ENV === 'development') {
        console.log('SettingsPage: State updated', {
          originalValues: newValues,
        });
      }
    }
  }, [settings]);

  const handleSave = async () => {
    try {
      updateSettings({
        max_parallel_instances:
          maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES,
        max_cpu_utilization_percent:
          maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION,
        max_memory_utilization_percent:
          maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION,
        gpu_acceleration_enabled: gpuAccelerationEnabled,
        system_total_cores: settings?.system_total_cores || 0,
        system_total_memory_mb: settings?.system_total_memory_mb || 0,
        calculations_directory: calculationsDirectory,
        timezone: timezone,
        gemini_api_key: geminiApiKey || null,
      });

      const newValues = {
        maxParallelInstances,
        maxCpuUtilization,
        maxMemoryUtilization,
        geminiApiKey,
        calculationsDirectory,
        timezone,
        gpuAccelerationEnabled,
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
      setGpuAccelerationEnabled(originalValues.gpuAccelerationEnabled ?? false);
    }
  };

  const handleCancel = () => {
    setMaxParallelInstances(originalValues.maxParallelInstances);
    setMaxCpuUtilization(originalValues.maxCpuUtilization);
    setMaxMemoryUtilization(originalValues.maxMemoryUtilization);
    setGeminiApiKey(originalValues.geminiApiKey || '');
    setCalculationsDirectory(originalValues.calculationsDirectory || '');
    setTimezone(originalValues.timezone || 'UTC');
    setGpuAccelerationEnabled(originalValues.gpuAccelerationEnabled ?? false);
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
    const currentGpuEnabled = gpuAccelerationEnabled ?? false;

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
    const originalGpuEnabled = originalValues.gpuAccelerationEnabled ?? false;

    const parallelChanged = currentParallel !== originalParallel;
    const cpuChanged = Math.abs(currentCpu - originalCpu) > 0.001;
    const memoryChanged = Math.abs(currentMemory - originalMemory) > 0.001;
    const geminiApiKeyChanged = currentGeminiApiKey !== originalGeminiApiKey;
    const calculationsDirectoryChanged =
      currentCalculationsDirectory !== originalCalculationsDirectory;
    const timezoneChanged = currentTimezone !== originalTimezone;
    const gpuEnabledChanged = currentGpuEnabled !== originalGpuEnabled;

    const hasChanges =
      parallelChanged ||
      cpuChanged ||
      memoryChanged ||
      geminiApiKeyChanged ||
      calculationsDirectoryChanged ||
      timezoneChanged ||
      gpuEnabledChanged;

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
          gpuEnabledChanged,
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
    gpuAccelerationEnabled,
    originalValues?.maxParallelInstances,
    originalValues?.maxCpuUtilization,
    originalValues?.maxMemoryUtilization,
    originalValues?.geminiApiKey,
    originalValues?.calculationsDirectory,
    originalValues?.timezone,
    originalValues?.gpuAccelerationEnabled,
  ]);

  const gpuInstallLabel = useMemo(() => {
    if (!gpuStatus) {
      return 'Install GPU4PySCF';
    }
    return gpuStatus.gpu4pyscf_installed
      ? 'Reinstall GPU4PySCF'
      : 'Install GPU4PySCF';
  }, [gpuStatus]);

  const canEnableGpuAcceleration = Boolean(
    gpuStatus?.gpu4pyscf_installed &&
      gpuStatus?.is_linux &&
      gpuStatus?.cuda_detected &&
      gpuStatus?.cuda_supported
  );

  const gpuStatusMessage = useMemo(() => {
    if (!gpuStatus) {
      return null;
    }

    if (!gpuStatus.is_linux) {
      return 'GPU4PySCF is supported on Linux only.';
    }

    if (!gpuStatus.cuda_detected) {
      const detail = gpuStatus.cuda_detection_message
        ? ` (${gpuStatus.cuda_detection_message})`
        : '';
      return `CUDA Toolkit not detected. Install CUDA 11/12/13 and ensure nvcc is in PATH.${detail}`;
    }

    if (!gpuStatus.cuda_supported) {
      const detail = gpuStatus.cuda_detection_message
        ? ` (${gpuStatus.cuda_detection_message})`
        : '';
      return `Detected CUDA ${gpuStatus.cuda_version || 'unknown'} is not supported. Supported versions: 11.x, 12.x, 13.x.${detail}`;
    }

    if (gpuStatus.gpu4pyscf_installed) {
      return 'GPU4PySCF is installed. You can reinstall if needed.';
    }

    return 'CUDA detected. Install GPU4PySCF to enable GPU acceleration.';
  }, [gpuStatus]);

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
          <div className={styles.sectionHeader}>
            <h3>GPU Acceleration (Linux)</h3>
            <label className={styles.gpuToggle}>
              <input
                type="checkbox"
                checked={gpuAccelerationEnabled}
                onChange={event =>
                  setGpuAccelerationEnabled(event.target.checked)
                }
                disabled={
                  isUpdating ||
                  isGpuStatusLoading ||
                  (!canEnableGpuAcceleration && !gpuAccelerationEnabled)
                }
              />
              <span className={styles.gpuToggleTrack}></span>
              <span className={styles.gpuToggleLabel}>Enabled</span>
            </label>
          </div>

          <div className={styles.settingItem}>
            <div className={styles.settingLabel}>
              <label>GPU4PySCF Installation</label>
              <p className={styles.settingHelp}>
                Detects CUDA Toolkit via <code>nvcc --version</code> and
                installs the matching GPU4PySCF package. cuTENSOR is strongly
                recommended and will be installed together.
              </p>
            </div>

            <div className={styles.settingControl}>
              {isGpuStatusLoading ? (
                <div className={styles.inlineStatus}>
                  <div className={styles.inlineSpinner}></div>
                  <span>Detecting CUDA environment...</span>
                </div>
              ) : (
                <>
                  <div className={styles.gpuStatusGrid}>
                    <div className={styles.gpuStatusItem}>
                      <span className={styles.infoLabel}>Platform</span>
                      <span className={styles.infoValue}>
                        {gpuStatus?.is_linux ? 'Linux' : 'Unsupported'}
                      </span>
                    </div>
                    <div className={styles.gpuStatusItem}>
                      <span className={styles.infoLabel}>CUDA Toolkit</span>
                      <span className={styles.infoValue}>
                        {gpuStatus?.cuda_version
                          ? `CUDA ${gpuStatus.cuda_version}`
                          : 'Not detected'}
                      </span>
                    </div>
                    <div className={styles.gpuStatusItem}>
                      <span className={styles.infoLabel}>CUDA Support</span>
                      <span
                        className={`${styles.gpuStatusBadge} ${
                          gpuStatus?.cuda_supported
                            ? styles.gpuStatusBadgeSuccess
                            : styles.gpuStatusBadgeWarning
                        }`}
                      >
                        {gpuStatus?.cuda_supported
                          ? 'Supported'
                          : 'Unsupported'}
                      </span>
                    </div>
                    <div className={styles.gpuStatusItem}>
                      <span className={styles.infoLabel}>GPU4PySCF</span>
                      <span
                        className={`${styles.gpuStatusBadge} ${
                          gpuStatus?.gpu4pyscf_installed
                            ? styles.gpuStatusBadgeSuccess
                            : styles.gpuStatusBadgeWarning
                        }`}
                      >
                        {gpuStatus?.gpu4pyscf_installed
                          ? `Installed${gpuStatus.gpu4pyscf_version ? ` v${gpuStatus.gpu4pyscf_version}` : ''}`
                          : 'Not installed'}
                      </span>
                    </div>
                    <div className={styles.gpuStatusItem}>
                      <span className={styles.infoLabel}>cuTENSOR</span>
                      <span
                        className={`${styles.gpuStatusBadge} ${
                          gpuStatus?.cutensor_installed
                            ? styles.gpuStatusBadgeSuccess
                            : styles.gpuStatusBadgeWarning
                        }`}
                      >
                        {gpuStatus?.cutensor_installed
                          ? `Installed${gpuStatus.cutensor_version ? ` v${gpuStatus.cutensor_version}` : ''}`
                          : 'Not installed'}
                      </span>
                    </div>
                    <div className={styles.gpuStatusItem}>
                      <span className={styles.infoLabel}>Recommended</span>
                      <span className={styles.gpuPackageValue}>
                        {gpuStatus?.recommended_gpu4pyscf_package ? (
                          <>
                            <code>
                              {gpuStatus.recommended_gpu4pyscf_package}
                            </code>
                            {gpuStatus.recommended_cutensor_package && (
                              <>
                                <span className={styles.gpuPackageDivider}>
                                  +
                                </span>
                                <code>
                                  {gpuStatus.recommended_cutensor_package}
                                </code>
                              </>
                            )}
                          </>
                        ) : (
                          'Not available'
                        )}
                      </span>
                    </div>
                  </div>

                  {gpuStatusMessage && (
                    <div className={styles.gpuNotice}>{gpuStatusMessage}</div>
                  )}

                  {gpuError && (
                    <div className={styles.gpuErrorBox}>
                      {gpuError instanceof Error
                        ? gpuError.message
                        : 'Failed to retrieve GPU status.'}
                    </div>
                  )}

                  <div className={styles.gpuActionRow}>
                    <button
                      className={styles.gpuActionButton}
                      onClick={() =>
                        installGpu4Pyscf({
                          include_cutensor: true,
                          force_reinstall: Boolean(
                            gpuStatus?.gpu4pyscf_installed
                          ),
                        })
                      }
                      disabled={
                        isGpuInstalling ||
                        isGpuStatusFetching ||
                        !gpuStatus?.is_linux ||
                        !gpuStatus?.cuda_supported
                      }
                    >
                      {isGpuInstalling ? (
                        <>
                          <div className={styles.buttonSpinner}></div>
                          Installing...
                        </>
                      ) : (
                        gpuInstallLabel
                      )}
                    </button>
                    <button
                      className={styles.gpuSecondaryButton}
                      onClick={() => refetchGpuStatus()}
                      disabled={isGpuStatusFetching || isGpuInstalling}
                    >
                      Refresh
                    </button>
                  </div>
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
