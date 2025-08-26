import React from 'react';
import { CalculationInstance } from '../types/api-types';
import styles from './Sidebar.module.css';

// Status display configuration
const STATUS_CONFIG = {
  pending: {
    icon: '⏳',
    label: 'Pending',
    className: styles.statusPending,
    color: '#ffa500',
  },
  running: {
    icon: '⚛️',
    label: 'Running',
    className: styles.statusRunning,
    color: '#2196f3',
  },
  completed: {
    icon: '✅',
    label: 'Completed',
    className: styles.statusCompleted,
    color: '#4caf50',
  },
  error: {
    icon: '❌',
    label: 'Error',
    className: styles.statusError,
    color: '#f44336',
  },
} as const;

type StatusType = keyof typeof STATUS_CONFIG;

interface SidebarProps {
  isOpen: boolean;
  onClose: () => void;
  calculations: CalculationInstance[];
  activeCalculationId: string | null;
  calculationsLoading: boolean;
  calculationsError: string | null;
  onCalculationSelect: (calculationId: string) => void;
  onCalculationDelete: (calculationId: string) => Promise<void>;
}

export const Sidebar: React.FC<SidebarProps> = ({
  isOpen,
  onClose,
  calculations,
  activeCalculationId,
  calculationsLoading,
  calculationsError,
  onCalculationSelect,
  onCalculationDelete,
}) => {
  return (
    <>
      {/* Backdrop/Overlay */}
      {isOpen && (
        <div
          className={styles.sidebarBackdrop}
          onClick={onClose}
          aria-label="Close sidebar"
        />
      )}

      {/* Sidebar Panel */}
      <aside className={`${styles.sidebar} ${isOpen ? styles.open : ''}`}>
        <div className={styles.sidebarContent}>
          <div className={styles.sidebarHeader}>
            <h3>Calculations</h3>
          </div>

          <div className={styles.sidebarCalculations}>
            {calculationsLoading ? (
              <div className={styles.sidebarLoading}>
                <p>Loading calculations...</p>
              </div>
            ) : calculationsError ? (
              <div className={styles.sidebarError}>
                <p>Error: {calculationsError}</p>
              </div>
            ) : calculations.length === 0 ? (
              <div className={styles.sidebarEmpty}>
                <p>No calculations yet</p>
                <p>Click the + button to create one</p>
              </div>
            ) : (
              calculations.map(calculation => (
                <div
                  key={calculation.id}
                  className={`${styles.sidebarCalculationItem} ${
                    calculation.id === activeCalculationId ? styles.active : ''
                  }`}
                  onClick={() => onCalculationSelect(calculation.id)}
                >
                  <div className={styles.calculationInfo}>
                    <div className={styles.calculationName}>{calculation.name}</div>
                    <div className={styles.calculationMeta}>
                      <div className={styles.calculationDate}>
                        {new Date(calculation.updatedAt).toLocaleDateString()}
                      </div>
                      <div className={styles.calculationStatus}>
                        <span
                          className={`${styles.statusBadge} ${STATUS_CONFIG[calculation.status as StatusType]?.className || styles.statusUnknown}`}
                          style={{
                            color:
                              STATUS_CONFIG[calculation.status as StatusType]
                                ?.color || '#666',
                            animation:
                              calculation.status === 'running'
                                ? 'pulse 2s infinite'
                                : 'none',
                          }}
                        >
                          <span className={styles.statusIcon}>
                            {STATUS_CONFIG[calculation.status as StatusType]
                              ?.icon || '❓'}
                          </span>
                          <span className={styles.statusLabel}>
                            {STATUS_CONFIG[calculation.status as StatusType]
                              ?.label || calculation.status}
                          </span>
                        </span>
                      </div>
                    </div>
                  </div>
                  <div className={styles.calculationActions}>
                    <button
                      className={styles.deleteBtn}
                      onClick={e => {
                        e.stopPropagation();
                        if (
                          confirm(`Delete calculation "${calculation.name}"?`)
                        ) {
                          onCalculationDelete(calculation.id);
                        }
                      }}
                    >
                      🗑️
                    </button>
                  </div>
                </div>
              ))
            )}
          </div>
        </div>
      </aside>
    </>
  );
};
