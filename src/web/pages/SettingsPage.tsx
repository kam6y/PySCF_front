import React, { useState, useEffect, useMemo } from 'react';
import { useAppSettings, useGetCalculations, useGpu4Pyscf } from '../hooks';
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
const DIRECTORY_CHANGE_BLOCKING_STATUSES = [
  'pending',
  'running',
  'waiting',
  'pausing',
] as const;

interface SettingsPageProps {
  // Props will be added when integrating with the main app
}

interface SettingsFormValues {
  maxParallelInstances?: number;
  maxCpuUtilization?: number;
  maxMemoryUtilization?: number;
  geminiApiKey: string;
  calculationsDirectory: string;
  timezone: Timezone;
  gpuAccelerationEnabled: boolean;
}

export const SettingsPage: React.FC<SettingsPageProps> = () => {
  const [formValues, setFormValues] = useState<SettingsFormValues>({
    geminiApiKey: '',
    calculationsDirectory: '',
    timezone: 'UTC',
    gpuAccelerationEnabled: false,
  });
  const [isSelectingFolder, setIsSelectingFolder] = useState(false);
  const [originalValues, setOriginalValues] = useState<SettingsFormValues>({
    geminiApiKey: '',
    calculationsDirectory: '',
    timezone: 'UTC',
    gpuAccelerationEnabled: false,
  });

  const { settings, isLoading, isUpdating, error, updateSettingsAsync } =
    useAppSettings();
  const { data: calculationsData } = useGetCalculations();
  const {
    status: gpuStatus,
    isLoading: isGpuStatusLoading,
    isFetching: isGpuStatusFetching,
    isInstalling: isGpuInstalling,
    error: gpuError,
    installGpu4Pyscf,
    refetch: refetchGpuStatus,
  } = useGpu4Pyscf();

  const updateFormValue = <K extends keyof SettingsFormValues>(
    key: K,
    value: SettingsFormValues[K]
  ) => {
    setFormValues(prev => ({ ...prev, [key]: value }));
  };

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

      const newValues: SettingsFormValues = {
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

      setFormValues(newValues);
      setOriginalValues(newValues);
    }
  }, [settings]);

  const hasBlockingCalculations = useMemo(
    () =>
      calculationsData?.calculations?.some(calculation =>
        DIRECTORY_CHANGE_BLOCKING_STATUSES.includes(
          calculation.status as (typeof DIRECTORY_CHANGE_BLOCKING_STATUSES)[number]
        )
      ) ?? false,
    [calculationsData]
  );

  const calculationsDirectoryChanged =
    formValues.calculationsDirectory !==
    (originalValues.calculationsDirectory || '');
  const blockSaveForDirectoryChange =
    hasBlockingCalculations && calculationsDirectoryChanged;

  const handleSave = async () => {
    if (blockSaveForDirectoryChange) {
      return;
    }

    try {
      await updateSettingsAsync({
        max_parallel_instances:
          formValues.maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES,
        max_cpu_utilization_percent:
          formValues.maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION,
        max_memory_utilization_percent:
          formValues.maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION,
        gpu_acceleration_enabled: formValues.gpuAccelerationEnabled,
        system_total_cores: settings?.system_total_cores || 0,
        system_total_memory_mb: settings?.system_total_memory_mb || 0,
        calculations_directory: formValues.calculationsDirectory,
        timezone: formValues.timezone,
        gemini_api_key: formValues.geminiApiKey || null,
      });

      const newValues = { ...formValues };
      setOriginalValues(newValues);
    } catch (error) {
      console.error('Failed to save settings:', error);
      setFormValues({ ...originalValues });
    }
  };

  const handleCancel = () => {
    setFormValues({ ...originalValues });
  };

  const handleSelectFolder = async () => {
    if (hasBlockingCalculations) {
      return;
    }

    setIsSelectingFolder(true);
    try {
      const result = await window.electronAPI.selectFolder();
      if (!result.canceled && result.filePath) {
        // Append /PySCF_calculations to the selected path
        const fullPath = `${result.filePath}/PySCF_calculations`;
        updateFormValue('calculationsDirectory', fullPath);
      }
    } catch (error) {
      console.error('Failed to select folder:', error);
    } finally {
      setIsSelectingFolder(false);
    }
  };

  const hasUnsavedChanges = useMemo(() => {
    // settings 未ロード時は変更なしとみなす
    if (!settings) return false;

    const floatClose = (
      a: number | undefined,
      b: number | undefined,
      defaultValue: number
    ) => Math.abs((a ?? defaultValue) - (b ?? defaultValue)) <= 0.001;

    if (
      (formValues.maxParallelInstances ?? DEFAULT_MAX_PARALLEL_INSTANCES) !==
      (originalValues.maxParallelInstances ?? DEFAULT_MAX_PARALLEL_INSTANCES)
    )
      return true;
    if (
      !floatClose(
        formValues.maxCpuUtilization,
        originalValues.maxCpuUtilization,
        DEFAULT_MAX_CPU_UTILIZATION
      )
    )
      return true;
    if (
      !floatClose(
        formValues.maxMemoryUtilization,
        originalValues.maxMemoryUtilization,
        DEFAULT_MAX_MEMORY_UTILIZATION
      )
    )
      return true;
    if (formValues.geminiApiKey !== (originalValues.geminiApiKey || ''))
      return true;
    if (
      formValues.calculationsDirectory !==
      (originalValues.calculationsDirectory || '')
    )
      return true;
    if (formValues.timezone !== (originalValues.timezone || 'UTC')) return true;
    if (
      (formValues.gpuAccelerationEnabled ?? false) !==
      (originalValues.gpuAccelerationEnabled ?? false)
    )
      return true;

    return false;
  }, [formValues, originalValues, settings]);

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
              disabled={isUpdating || blockSaveForDirectoryChange}
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
                  value={
                    formValues.maxParallelInstances ||
                    DEFAULT_MAX_PARALLEL_INSTANCES
                  }
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    updateFormValue('maxParallelInstances', newValue);
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
                  {formValues.maxParallelInstances ||
                    DEFAULT_MAX_PARALLEL_INSTANCES}
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
                  value={
                    formValues.maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION
                  }
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    updateFormValue('maxCpuUtilization', newValue);
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
                    formValues.maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION
                  ).toFixed(0)}
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
                  value={
                    formValues.maxMemoryUtilization ||
                    DEFAULT_MAX_MEMORY_UTILIZATION
                  }
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    updateFormValue('maxMemoryUtilization', newValue);
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
                    formValues.maxMemoryUtilization ||
                    DEFAULT_MAX_MEMORY_UTILIZATION
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
                checked={formValues.gpuAccelerationEnabled}
                onChange={event =>
                  updateFormValue(
                    'gpuAccelerationEnabled',
                    event.target.checked
                  )
                }
                disabled={
                  isUpdating ||
                  isGpuStatusLoading ||
                  (!canEnableGpuAcceleration &&
                    !formValues.gpuAccelerationEnabled)
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
                    {gpuStatus?.is_linux && (
                      <>
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
                      </>
                    )}
                  </div>

                  {gpuStatus?.is_linux && gpuStatusMessage && (
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
                Required API key for the Gemini-based AI agent. Chat assistance
                is unavailable until a valid key is configured.
              </p>
            </div>

            <div className={styles.settingControl}>
              <div className={styles.textInputContainer}>
                <input
                  id="geminiApiKey"
                  type="password"
                  placeholder="Enter your Gemini API key..."
                  value={formValues.geminiApiKey}
                  onChange={e => {
                    const newValue = e.target.value;
                    updateFormValue('geminiApiKey', newValue);
                  }}
                  className={styles.textInput}
                  disabled={isUpdating}
                />
                <div className={styles.inputStatus}>
                  {formValues.geminiApiKey ? (
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
                  value={formValues.timezone}
                  onChange={e => {
                    const newValue = e.target.value as Timezone;
                    updateFormValue('timezone', newValue);
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
                    {getTimezoneLabel(formValues.timezone)}
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
                  value={formValues.calculationsDirectory}
                  readOnly
                  className={styles.textInput}
                  style={{
                    fontFamily: 'Monaco, Menlo, Courier New, monospace',
                  }}
                />
                <button
                  onClick={handleSelectFolder}
                  className={styles.selectButton}
                  disabled={
                    isUpdating || isSelectingFolder || hasBlockingCalculations
                  }
                >
                  {isSelectingFolder ? 'Selecting...' : 'Change Folder...'}
                </button>
              </div>
              {hasBlockingCalculations && (
                <div className={styles.warningBox}>
                  <span className={styles.warningIcon}>⚠</span>
                  <span className={styles.warningText}>
                    Calculation data folder cannot be changed while calculations
                    are running or queued.
                  </span>
                </div>
              )}
              {formValues.calculationsDirectory !==
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
