import React, { useState, useEffect, useMemo, useRef } from 'react';
import { useAppSettings } from '../hooks';
import styles from './SettingsPage.module.css';

// Default values constants
const DEFAULT_MAX_PARALLEL_INSTANCES = 4;
const DEFAULT_MAX_CPU_UTILIZATION = 80.0;
const DEFAULT_MAX_MEMORY_UTILIZATION = 80.0;

interface SettingsPageProps {
  // Props will be added when integrating with the main app
}

export const SettingsPage: React.FC<SettingsPageProps> = () => {
  const [maxParallelInstances, setMaxParallelInstances] = useState<number | undefined>(undefined);
  const [maxCpuUtilization, setMaxCpuUtilization] = useState<number | undefined>(undefined);
  const [maxMemoryUtilization, setMaxMemoryUtilization] =
    useState<number | undefined>(undefined);
  const [originalValues, setOriginalValues] = useState<{
    maxParallelInstances?: number;
    maxCpuUtilization?: number;
    maxMemoryUtilization?: number;
  }>({});

  const { settings, isLoading, isUpdating, error, updateSettings } =
    useAppSettings();

  // Update local state when settings are loaded
  useEffect(() => {
    if (settings) {
      // Explicit type conversion and validation
      const maxParallelInstancesValue = Number(settings.max_parallel_instances);
      const maxCpuUtilizationValue = Number(settings.max_cpu_utilization_percent);
      const maxMemoryUtilizationValue = Number(settings.max_memory_utilization_percent);
      
      const newValues = {
        maxParallelInstances: !isNaN(maxParallelInstancesValue) ? maxParallelInstancesValue : DEFAULT_MAX_PARALLEL_INSTANCES,
        maxCpuUtilization: !isNaN(maxCpuUtilizationValue) ? maxCpuUtilizationValue : DEFAULT_MAX_CPU_UTILIZATION,
        maxMemoryUtilization: !isNaN(maxMemoryUtilizationValue) ? maxMemoryUtilizationValue : DEFAULT_MAX_MEMORY_UTILIZATION,
      };

      if (process.env.NODE_ENV === 'development') {
        console.log('SettingsPage: Loading settings', { settings, newValues });
      }

      setMaxParallelInstances(newValues.maxParallelInstances);
      setMaxCpuUtilization(newValues.maxCpuUtilization);
      setMaxMemoryUtilization(newValues.maxMemoryUtilization);
      setOriginalValues(newValues);
      
      if (process.env.NODE_ENV === 'development') {
        console.log('SettingsPage: State updated', { originalValues: newValues });
      }
    }
  }, [settings]);

  const handleSave = async () => {
    try {
      updateSettings({
        max_parallel_instances: maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES,
        max_cpu_utilization_percent: maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION,
        max_memory_utilization_percent: maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION,
        system_total_cores: settings?.system_total_cores || 0,
        system_total_memory_mb: settings?.system_total_memory_mb || 0,
      });

      const newValues = {
        maxParallelInstances,
        maxCpuUtilization,
        maxMemoryUtilization,
      };
      setOriginalValues(newValues);
    } catch (error) {
      console.error('Failed to save settings:', error);
      // Reset to original values on error
      setMaxParallelInstances(originalValues.maxParallelInstances);
      setMaxCpuUtilization(originalValues.maxCpuUtilization);
      setMaxMemoryUtilization(originalValues.maxMemoryUtilization);
    }
  };

  const handleCancel = () => {
    setMaxParallelInstances(originalValues.maxParallelInstances);
    setMaxCpuUtilization(originalValues.maxCpuUtilization);
    setMaxMemoryUtilization(originalValues.maxMemoryUtilization);
  };

  const hasUnsavedChanges = useMemo(() => {
    // Return false if settings haven't been loaded yet
    if (!originalValues || Object.keys(originalValues).length === 0) {
      return false;
    }

    // Use default values for comparison if values are undefined
    const currentParallel = maxParallelInstances ?? DEFAULT_MAX_PARALLEL_INSTANCES;
    const currentCpu = maxCpuUtilization ?? DEFAULT_MAX_CPU_UTILIZATION;
    const currentMemory = maxMemoryUtilization ?? DEFAULT_MAX_MEMORY_UTILIZATION;
    
    const originalParallel = originalValues.maxParallelInstances ?? DEFAULT_MAX_PARALLEL_INSTANCES;
    const originalCpu = originalValues.maxCpuUtilization ?? DEFAULT_MAX_CPU_UTILIZATION;
    const originalMemory = originalValues.maxMemoryUtilization ?? DEFAULT_MAX_MEMORY_UTILIZATION;

    const parallelChanged = currentParallel !== originalParallel;
    const cpuChanged = Math.abs(currentCpu - originalCpu) > 0.001;
    const memoryChanged = Math.abs(currentMemory - originalMemory) > 0.001;

    const hasChanges = parallelChanged || cpuChanged || memoryChanged;

    // Debug logging in development
    if (process.env.NODE_ENV === 'development') {
      console.log('SettingsPage: hasUnsavedChanges check', {
        current: { currentParallel, currentCpu, currentMemory },
        original: { originalParallel, originalCpu, originalMemory },
        changes: { parallelChanged, cpuChanged, memoryChanged },
        hasChanges
      });
    }

    return hasChanges;
  }, [
    maxParallelInstances,
    maxCpuUtilization,
    maxMemoryUtilization,
    originalValues?.maxParallelInstances,
    originalValues?.maxCpuUtilization,
    originalValues?.maxMemoryUtilization
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
                  value={maxParallelInstances || DEFAULT_MAX_PARALLEL_INSTANCES}
                  onChange={e => {
                    const newValue = Number(e.target.value);
                    if (process.env.NODE_ENV === 'development') {
                      console.log('SettingsPage: maxParallelInstances changed', newValue);
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
                      console.log('SettingsPage: maxCpuUtilization changed', newValue);
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
                  {(maxCpuUtilization || DEFAULT_MAX_CPU_UTILIZATION).toFixed(0)}%
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
                      console.log('SettingsPage: maxMemoryUtilization changed', newValue);
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
                  {(maxMemoryUtilization || DEFAULT_MAX_MEMORY_UTILIZATION).toFixed(0)}%
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

      </div>
    </div>
  );
};
