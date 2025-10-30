import React, { useState } from 'react';
import styles from './InitialSetupDialog.module.css';

interface InitialSetupDialogProps {
  onComplete: (calculationsDirectory: string) => void;
  defaultDirectory: string;
}

export const InitialSetupDialog: React.FC<InitialSetupDialogProps> = ({
  onComplete,
  defaultDirectory,
}) => {
  const [selectedPath, setSelectedPath] = useState<string>(defaultDirectory);
  const [isSelecting, setIsSelecting] = useState(false);

  const handleSelectFolder = async () => {
    setIsSelecting(true);
    try {
      const result = await window.electronAPI.selectFolder();
      if (!result.canceled && result.filePath) {
        // Append /PySCF_calculations to the selected path
        const fullPath = `${result.filePath}/PySCF_calculations`;
        setSelectedPath(fullPath);
      }
    } catch (error) {
      console.error('Failed to select folder:', error);
    } finally {
      setIsSelecting(false);
    }
  };

  const handleUseDefault = () => {
    onComplete(defaultDirectory);
  };

  const handleContinue = () => {
    onComplete(selectedPath);
  };

  return (
    <div className={styles.overlay}>
      <div className={styles.dialog}>
        <div className={styles.header}>
          <h2>Welcome to PySCF Native App</h2>
        </div>

        <div className={styles.content}>
          <p className={styles.description}>
            Choose where to store calculation data. A 'PySCF_calculations'
            subfolder will be created in the selected directory.
          </p>

          <div className={styles.pathSection}>
            <div className={styles.pathDisplay}>
              <span className={styles.pathIcon}>📁</span>
              <span className={styles.pathText}>{selectedPath}</span>
            </div>

            <button
              onClick={handleSelectFolder}
              className={styles.selectButton}
              disabled={isSelecting}
            >
              {isSelecting ? 'Selecting...' : 'Choose Folder...'}
            </button>
          </div>

          <p className={styles.note}>
            You can change this location later in Settings.
          </p>
        </div>

        <div className={styles.actions}>
          <button onClick={handleUseDefault} className={styles.defaultButton}>
            Use Default
          </button>
          <button onClick={handleContinue} className={styles.continueButton}>
            Continue
          </button>
        </div>
      </div>
    </div>
  );
};
