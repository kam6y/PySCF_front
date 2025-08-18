import React from 'react';
import { CalculationInstance } from '../types/api-types';

// Status display configuration
const STATUS_CONFIG = {
  pending: {
    icon: '⏳',
    label: 'Pending',
    className: 'status-pending',
    color: '#ffa500',
  },
  running: {
    icon: '⚛️',
    label: 'Running',
    className: 'status-running',
    color: '#2196f3',
  },
  completed: {
    icon: '✅',
    label: 'Completed',
    className: 'status-completed',
    color: '#4caf50',
  },
  error: {
    icon: '❌',
    label: 'Error',
    className: 'status-error',
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
          className="sidebar-backdrop"
          onClick={onClose}
          aria-label="Close sidebar"
        />
      )}

      {/* Sidebar Panel */}
      <aside className={`sidebar ${isOpen ? 'open' : ''}`}>
        <div className="sidebar-content">
          <div className="sidebar-header">
            <h3>Calculations</h3>
          </div>

          <div className="sidebar-calculations">
            {calculationsLoading ? (
              <div className="sidebar-loading">
                <p>Loading calculations...</p>
              </div>
            ) : calculationsError ? (
              <div className="sidebar-error">
                <p>Error: {calculationsError}</p>
              </div>
            ) : calculations.length === 0 ? (
              <div className="sidebar-empty">
                <p>No calculations yet</p>
                <p>Click the + button to create one</p>
              </div>
            ) : (
              calculations.map(calculation => (
                <div
                  key={calculation.id}
                  className={`sidebar-calculation-item ${
                    calculation.id === activeCalculationId ? 'active' : ''
                  }`}
                  onClick={() => onCalculationSelect(calculation.id)}
                >
                  <div className="calculation-info">
                    <div className="calculation-name">{calculation.name}</div>
                    <div className="calculation-meta">
                      <div className="calculation-date">
                        {new Date(calculation.updatedAt).toLocaleDateString()}
                      </div>
                      <div className="calculation-status">
                        <span
                          className={`status-badge ${STATUS_CONFIG[calculation.status as StatusType]?.className || 'status-unknown'}`}
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
                          <span className="status-icon">
                            {STATUS_CONFIG[calculation.status as StatusType]
                              ?.icon || '❓'}
                          </span>
                          <span className="status-label">
                            {STATUS_CONFIG[calculation.status as StatusType]
                              ?.label || calculation.status}
                          </span>
                        </span>
                      </div>
                    </div>
                  </div>
                  <div className="calculation-actions">
                    <button
                      className="delete-btn"
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
