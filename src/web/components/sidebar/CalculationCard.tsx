import React from 'react';
import { useAppSettings } from '../../hooks';
import { CalculationSummary } from '../../types/api-types';
import { formatDateTimeWithSeconds } from '../../utils/dateFormatter';
import styles from '../Sidebar.module.css';

export interface CalculationCardProps {
  calculation: CalculationSummary;
  isActive: boolean;
  onSelect: (calculationId: string) => void;
  onRequestDelete: (calculationId: string, calculationName: string) => void;
}

export const CalculationCard = React.memo<CalculationCardProps>(({
  calculation,
  isActive,
  onSelect,
  onRequestDelete,
}) => {
  const { settings } = useAppSettings();

  const handleCardClick = () => {
    onSelect(calculation.id);
  };

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
            {formatDateTimeWithSeconds(
              calculation.created_at || calculation.date,
              settings?.timezone || 'UTC'
            )}
          </div>
        </div>

        {(calculation.calculation_method || calculation.basis_function || calculation.exchange_correlation) && (
          <div className={styles.calculationTags}>
            {calculation.calculation_method && (
              <span className={`${styles.calculationTag} ${styles.tagMethod}`}>
                {calculation.calculation_method}
              </span>
            )}
            {calculation.basis_function && (
              <span className={`${styles.calculationTag} ${styles.tagBasis}`}>
                {calculation.basis_function}
              </span>
            )}
            {calculation.exchange_correlation && (
              <span className={`${styles.calculationTag} ${styles.tagXC}`}>
                {calculation.exchange_correlation}
              </span>
            )}
          </div>
        )}
      </div>
      <div className={styles.calculationActions}>
        <button
          className={styles.deleteBtn}
          onClick={e => {
            e.stopPropagation();
            onRequestDelete(calculation.id, calculation.name);
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
});
