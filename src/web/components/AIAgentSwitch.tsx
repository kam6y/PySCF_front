import React from 'react';
import styles from './AIAgentSwitch.module.css';

interface AIAgentSwitchProps {
  isEnabled: boolean;
  onChange: (enabled: boolean) => void;
  disabled?: boolean;
}

export const AIAgentSwitch: React.FC<AIAgentSwitchProps> = ({
  isEnabled,
  onChange,
  disabled = false,
}) => {
  const handleClick = () => {
    if (!disabled) {
      onChange(!isEnabled);
    }
  };

  const handleKeyDown = (event: React.KeyboardEvent) => {
    if ((event.key === 'Enter' || event.key === ' ') && !disabled) {
      event.preventDefault();
      onChange(!isEnabled);
    }
  };

  return (
    <button
      className={`${styles.switchContainer} ${isEnabled ? styles.enabled : ''} ${disabled ? styles.disabled : ''}`}
      onClick={handleClick}
      onKeyDown={handleKeyDown}
      aria-pressed={isEnabled}
      aria-label="Toggle AI Agent"
      disabled={disabled}
      type="button"
    >
      <div className={styles.switchBackground}>
        <div className={styles.switchToggle}>
          <span className={styles.aiText}>AI</span>
        </div>
      </div>
    </button>
  );
};
