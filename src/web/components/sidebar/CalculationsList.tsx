import React from 'react';
import { CalculationSummary } from '../../types/api-types';
import styles from '../Sidebar.module.css';
import { CalculationCard } from './CalculationCard';
import { groupCalculationsByStatus } from './StatusConfig';

interface CalculationsListProps {
  calculations: CalculationSummary[];
  activeCalculationId: string | null;
  calculationsLoading: boolean;
  calculationsError: string | null;
  onCalculationSelect: (calculationId: string) => void;
  onRequestDelete: (calculationId: string, calculationName: string) => void;
  onBulkDeleteError: (errorCalculations: CalculationSummary[]) => void;
}

export const CalculationsList: React.FC<CalculationsListProps> = ({
  calculations,
  activeCalculationId,
  calculationsLoading,
  calculationsError,
  onCalculationSelect,
  onRequestDelete,
  onBulkDeleteError,
}) => {
  if (calculationsLoading) {
    return (
      <div className={styles.sidebarLoading}>
        <p>Loading calculations...</p>
      </div>
    );
  }

  if (calculationsError) {
    return (
      <div className={styles.sidebarError}>
        <p>Error: {calculationsError}</p>
      </div>
    );
  }

  if (calculations.length === 0) {
    return (
      <div className={styles.sidebarEmpty}>
        <p>No calculations yet</p>
        <p>Click the + button to create one</p>
      </div>
    );
  }

  return (
    <>
      {groupCalculationsByStatus(calculations).map(group => (
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
            {group.status === 'error' && group.calculations.length > 0 && (
              <button
                className={styles.bulkDeleteButton}
                onClick={e => {
                  e.stopPropagation();
                  onBulkDeleteError(group.calculations);
                }}
                title={`Bulk delete ${group.calculations.length} error calculations`}
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
                onRequestDelete={onRequestDelete}
              />
            ))}
          </div>
        </div>
      ))}
    </>
  );
};
