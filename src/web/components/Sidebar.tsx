import React, { useState } from 'react';
import { CalculationSummary } from '../types/api-types';
import { useGetCalculationDetails } from '../hooks/useCalculationQueries';
import styles from './Sidebar.module.css';

// Status display configuration
const STATUS_CONFIG = {
  pending: {
    icon: (
      <svg
        width="20"
        height="20"
        viewBox="0 0 16 16"
        fill="none"
        xmlns="http://www.w3.org/2000/svg"
      >
        <path
          d="M8 1C4.13401 1 1 4.13401 1 8C1 11.866 4.13401 15 8 15C11.866 15 15 11.866 15 8C15 4.13401 11.866 1 8 1Z"
          stroke="currentColor"
          strokeWidth="1.5"
          fill="none"
        />
        <path
          d="M8 4V8L10.5 10.5"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
      </svg>
    ),
    label: 'Pending',
    className: styles.statusPending,
    color: '#ffa500',
  },
  running: {
    icon: (
      <svg
        width="20"
        height="20"
        viewBox="0 0 16 16"
        fill="none"
        xmlns="http://www.w3.org/2000/svg"
        className={styles.animateSpin}
      >
        <path
          d="M8 1V3M8 13V15M3.05 3.05L4.46 4.46M11.54 11.54L12.95 12.95M1 8H3M13 8H15M3.05 12.95L4.46 11.54M11.54 4.46L12.95 3.05"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
      </svg>
    ),
    label: 'Running',
    className: styles.statusRunning,
    color: '#2196f3',
  },
  completed: {
    icon: (
      <svg
        width="20"
        height="20"
        viewBox="0 0 16 16"
        fill="none"
        xmlns="http://www.w3.org/2000/svg"
      >
        <path
          d="M8 1C4.13401 1 1 4.13401 1 8C1 11.866 4.13401 15 8 15C11.866 15 15 11.866 15 8C15 4.13401 11.866 1 8 1Z"
          stroke="currentColor"
          strokeWidth="1.5"
          fill="none"
        />
        <path
          d="M11 6L7 10L5 8"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
      </svg>
    ),
    label: 'Completed',
    className: styles.statusCompleted,
    color: '#4caf50',
  },
  error: {
    icon: (
      <svg
        width="20"
        height="20"
        viewBox="0 0 16 16"
        fill="none"
        xmlns="http://www.w3.org/2000/svg"
      >
        <path
          d="M8 1C4.13401 1 1 4.13401 1 8C1 11.866 4.13401 15 8 15C11.866 15 15 11.866 15 8C15 4.13401 11.866 1 8 1Z"
          stroke="currentColor"
          strokeWidth="1.5"
          fill="none"
        />
        <path
          d="M10 6L6 10M6 6L10 10"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
      </svg>
    ),
    label: 'Error',
    className: styles.statusError,
    color: '#f44336',
  },
  waiting: {
    icon: (
      <svg
        width="20"
        height="20"
        viewBox="0 0 16 16"
        fill="none"
        xmlns="http://www.w3.org/2000/svg"
      >
        <path
          d="M8 1C4.13401 1 1 4.13401 1 8C1 11.866 4.13401 15 8 15C11.866 15 15 11.866 15 8C15 4.13401 11.866 1 8 1Z"
          stroke="currentColor"
          strokeWidth="1.5"
          fill="none"
        />
        <path
          d="M8 4V8L10.5 10.5"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
        <path
          d="M12 3L13 2M13 2L14 3M13 2V1"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
      </svg>
    ),
    label: 'Waiting',
    className: styles.statusWaiting,
    color: '#ff9800',
  },
} as const;

// Status order for grouping
const STATUS_ORDER: (keyof typeof STATUS_CONFIG)[] = [
  'error',
  'pending',
  'waiting',
  'running',
  'completed',
];

// Group calculations by status
const groupCalculationsByStatus = (calculations: CalculationSummary[]) => {
  const grouped = calculations.reduce(
    (acc, calculation) => {
      const status = calculation.status as keyof typeof STATUS_CONFIG;
      if (!acc[status]) {
        acc[status] = [];
      }
      acc[status].push(calculation);
      return acc;
    },
    {} as Record<keyof typeof STATUS_CONFIG, CalculationSummary[]>
  );

  // Return groups in the specified order
  return STATUS_ORDER.map(status => ({
    status,
    config: STATUS_CONFIG[status],
    calculations: grouped[status] || [],
  })).filter(group => group.calculations.length > 0);
};

// Individual calculation card component
interface CalculationCardProps {
  calculation: CalculationSummary;
  isActive: boolean;
  onSelect: (calculationId: string) => void;
  onDelete: (calculationId: string) => Promise<void>;
}

const CalculationCard: React.FC<CalculationCardProps> = ({
  calculation,
  isActive,
  onSelect,
  onDelete,
}) => {
  const { data: detailsData, isLoading: detailsLoading } =
    useGetCalculationDetails(calculation.id);

  const handleCardClick = () => {
    onSelect(calculation.id);
  };

  const calculation_details = detailsData?.calculation?.parameters;

  return (
    <div
      className={`${styles.sidebarCalculationItem} ${
        isActive ? styles.active : ''
      }`}
      onClick={handleCardClick}
    >
      <div className={styles.calculationInfo}>
        <div className={styles.calculationName}>{calculation.name}</div>
        <div className={styles.calculationMeta}>
          <div className={styles.calculationDate}>
            {new Date(calculation.date).toLocaleString('en-US', {
              year: 'numeric',
              month: '2-digit',
              day: '2-digit',
              hour: '2-digit',
              minute: '2-digit',
              second: '2-digit',
            })}
          </div>
        </div>

        {/* Calculation Details - Always visible */}
        {calculation_details && !detailsLoading ? (
          <div className={styles.calculationTags}>
            <span className={`${styles.calculationTag} ${styles.tagMethod}`}>
              {calculation_details.calculation_method}
            </span>
            <span className={`${styles.calculationTag} ${styles.tagBasis}`}>
              {calculation_details.basis_function}
            </span>
            {calculation_details.exchange_correlation && (
              <span className={`${styles.calculationTag} ${styles.tagXC}`}>
                {calculation_details.exchange_correlation}
              </span>
            )}
          </div>
        ) : detailsLoading ? (
          <div className={styles.calculationTagsLoading}>
            <div className={styles.tagSkeleton}></div>
            <div className={styles.tagSkeleton}></div>
            <div className={styles.tagSkeleton}></div>
          </div>
        ) : null}
      </div>
      <div className={styles.calculationActions}>
        <button
          className={styles.deleteBtn}
          onClick={e => {
            e.stopPropagation();
            if (confirm(`Delete calculation "${calculation.name}"?`)) {
              onDelete(calculation.id);
            }
          }}
        >
          <svg
            width="14"
            height="14"
            viewBox="0 0 24 24"
            fill="none"
            stroke="currentColor"
            strokeWidth="2"
            strokeLinecap="round"
            strokeLinejoin="round"
          >
            <polyline points="3,6 5,6 21,6"></polyline>
            <path d="m5,6 1,14 c0,1 1,2 2,2 h8 c1,0 2,-1 2,-2 l1,-14"></path>
            <path d="m10,11 v6"></path>
            <path d="m14,11 v6"></path>
            <path d="m7,6 V4 c0,-1 1,-2 2,-2 h6 c1,0 2,1 2,2 v2"></path>
          </svg>
        </button>
      </div>
    </div>
  );
};

interface SidebarProps {
  isOpen: boolean;
  onClose: () => void;
  calculations: CalculationSummary[];
  activeCalculationId: string | null;
  calculationsLoading: boolean;
  calculationsError: string | null;
  onCalculationSelect: (calculationId: string) => void;
  onCalculationDelete: (calculationId: string) => Promise<void>;
  onCreateNew: () => void;
  searchQuery: string;
  onSearchChange: (query: string) => void;
  onUserMenuToggle: () => void;
  isUserMenuOpen: boolean;
  onSettingsOpen: () => void;
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
  onCreateNew,
  searchQuery,
  onSearchChange,
  onUserMenuToggle,
  isUserMenuOpen,
  onSettingsOpen,
}) => {
  // Bulk delete handler for error instances
  const handleBulkDeleteError = async (
    errorCalculations: CalculationSummary[]
  ) => {
    if (errorCalculations.length === 0) return;

    const count = errorCalculations.length;
    const confirmed = confirm(
      `Do you want to delete ${count} error instances? This operation cannot be undone.`
    );

    if (!confirmed) return;

    try {
      // Delete all error calculations in parallel
      await Promise.all(
        errorCalculations.map(calculation =>
          onCalculationDelete(calculation.id)
        )
      );
    } catch (error) {
      console.error('一括削除中にエラーが発生しました:', error);
      alert('Failed to delete some instances.');
    }
  };

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
          {/* Top Section - Create button and Search */}
          <div className={styles.sidebarTopSection}>
            <button className={styles.createNewButton} onClick={onCreateNew}>
              + Make a instance
            </button>
            <div className={styles.searchContainer}>
              <svg
                className={styles.searchIcon}
                width="16"
                height="16"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              >
                <circle cx="11" cy="11" r="8"></circle>
                <path d="m21 21-4.35-4.35"></path>
              </svg>
              <input
                type="text"
                placeholder="Value"
                value={searchQuery}
                onChange={e => onSearchChange(e.target.value)}
                className={styles.searchInput}
              />
              <button
                className={styles.searchClearButton}
                onClick={e => {
                  e.preventDefault();
                  onSearchChange('');
                }}
                style={{ display: searchQuery ? 'block' : 'none' }}
              >
                ×
              </button>
            </div>
          </div>

          {/* Main Content Area */}
          <div className={styles.sidebarMainContent}>
            <div className={styles.sidebarHeader}>
              <h3>Instances</h3>
            </div>

            <div className={styles.sidebarCalculations}>
              {calculationsLoading ? (
                <div className={styles.sidebarLoading}>
                  <p>Loading instances...</p>
                </div>
              ) : calculationsError ? (
                <div className={styles.sidebarError}>
                  <p>Error: {calculationsError}</p>
                </div>
              ) : calculations.length === 0 ? (
                <div className={styles.sidebarEmpty}>
                  <p>No instances yet</p>
                  <p>Click the + button to create one</p>
                </div>
              ) : (
                groupCalculationsByStatus(calculations).map(group => (
                  <div key={group.status} className={styles.statusSection}>
                    <div className={styles.statusSectionHeader}>
                      <div className={styles.statusSectionInfo}>
                        <span className={styles.statusSectionIcon}>
                          {group.config.icon}
                        </span>
                        <h4 className={styles.statusSectionTitle}>
                          {group.config.label}
                        </h4>
                      </div>
                      {group.status === 'error' &&
                        group.calculations.length > 0 && (
                          <button
                            className={styles.bulkDeleteButton}
                            onClick={e => {
                              e.stopPropagation();
                              handleBulkDeleteError(group.calculations);
                            }}
                            title={`Bulk delete ${group.calculations.length} error instances`}
                          >
                            <svg
                              width="14"
                              height="14"
                              viewBox="0 0 24 24"
                              fill="none"
                              stroke="currentColor"
                              strokeWidth="2"
                              strokeLinecap="round"
                              strokeLinejoin="round"
                            >
                              <polyline points="3,6 5,6 21,6"></polyline>
                              <path d="m5,6 1,14 c0,1 1,2 2,2 h8 c1,0 2,-1 2,-2 l1,-14"></path>
                              <path d="m10,11 v6"></path>
                              <path d="m14,11 v6"></path>
                              <path d="m7,6 V4 c0,-1 1,-2 2,-2 h6 c1,0 2,1 2,2 v2"></path>
                            </svg>
                          </button>
                        )}
                    </div>
                    <div className={styles.statusSectionCalculations}>
                      {group.calculations.map(calculation => (
                        <CalculationCard
                          key={calculation.id}
                          calculation={calculation}
                          isActive={calculation.id === activeCalculationId}
                          onSelect={onCalculationSelect}
                          onDelete={onCalculationDelete}
                        />
                      ))}
                    </div>
                  </div>
                ))
              )}
            </div>
          </div>

          {/* Bottom Section - User Info */}
          <div className={styles.sidebarBottomSection}>
            <div
              className={`${styles.userInfoSection} ${isUserMenuOpen ? styles.userMenuOpen : ''}`}
              onClick={onUserMenuToggle}
            >
              {isUserMenuOpen && (
                <div className={styles.userMenu}>
                  <button
                    className={styles.userMenuItem}
                    onClick={onSettingsOpen}
                  >
                    Settings
                  </button>
                  <button className={styles.userMenuItem}>Profile</button>
                  <button className={styles.userMenuItem}>About</button>
                </div>
              )}
              <div className={styles.userInfoMain}>
                <span className={styles.userName}>User name</span>
                <span
                  className={`${styles.userMenuToggle} ${isUserMenuOpen ? styles.rotated : ''}`}
                >
                  <svg
                    width="24"
                    height="24"
                    viewBox="0 0 16 16"
                    fill="none"
                    xmlns="http://www.w3.org/2000/svg"
                  >
                    <path
                      d="M4 6L8 10L12 6"
                      stroke="currentColor"
                      strokeWidth="1.5"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    />
                  </svg>
                </span>
              </div>
            </div>
          </div>
        </div>
      </aside>
    </>
  );
};
