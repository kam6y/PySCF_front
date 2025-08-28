import React, { useState, useEffect } from 'react';
import { useAppSettings } from '../hooks';
import styles from './SettingsPage.module.css';

interface SettingsPageProps {
  // Props will be added when integrating with the main app
}

export const SettingsPage: React.FC<SettingsPageProps> = () => {
  const [maxParallelInstances, setMaxParallelInstances] = useState<number>(4);
  const [originalValue, setOriginalValue] = useState<number>(4);

  const {
    settings,
    isLoading,
    isUpdating,
    error,
    updateSettings,
  } = useAppSettings();

  // Update local state when settings are loaded
  useEffect(() => {
    if (settings && settings.max_parallel_instances) {
      setMaxParallelInstances(settings.max_parallel_instances);
      setOriginalValue(settings.max_parallel_instances);
    }
  }, [settings]);

  const handleSave = async () => {
    try {
      await updateSettings({
        max_parallel_instances: maxParallelInstances,
      });
      setOriginalValue(maxParallelInstances);
    } catch (error) {
      console.error('Failed to save settings:', error);
      // Reset to original value on error
      setMaxParallelInstances(originalValue);
    }
  };

  const handleCancel = () => {
    setMaxParallelInstances(originalValue);
  };

  const hasUnsavedChanges = maxParallelInstances !== originalValue;

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
        <h2>Settings</h2>
        <p className={styles.settingsDescription}>
          Configure application behavior and performance settings
        </p>
      </div>

      <div className={styles.settingsContent}>
        {/* Error Display */}
        {error && (
          <div className={styles.errorContainer}>
            <div className={styles.errorIcon}>⚠️</div>
            <div className={styles.errorContent}>
              <h4>Error Loading Settings</h4>
              <p>{error instanceof Error ? error.message : 'An unexpected error occurred'}</p>
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
                  value={maxParallelInstances}
                  onChange={(e) => setMaxParallelInstances(parseInt(e.target.value))}
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
                  {maxParallelInstances}
                </span>
                <span className={styles.valueUnit}>calculations</span>
              </div>
            </div>
          </div>
        </div>

        {/* Action buttons */}
        {hasUnsavedChanges && (
          <div className={styles.settingsActions}>
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
    </div>
  );
};