import { CalculationSummary } from '../../types/api-types';
import styles from '../Sidebar.module.css';

export const STATUS_CONFIG = {
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
  },
  pausing: {
    icon: (
      <span
        style={{
          position: 'relative',
          display: 'inline-block',
          width: '20px',
          height: '20px',
        }}
      >
        <svg
          width="20"
          height="20"
          viewBox="0 0 16 16"
          fill="none"
          xmlns="http://www.w3.org/2000/svg"
          className={styles.animateSpin}
          style={{ position: 'absolute', top: 0, left: 0 }}
        >
          <circle
            cx="8"
            cy="8"
            r="6"
            stroke="currentColor"
            strokeWidth="1.5"
            strokeDasharray="30 8"
            strokeLinecap="round"
            fill="none"
          />
        </svg>
        <svg
          width="20"
          height="20"
          viewBox="0 0 16 16"
          fill="none"
          xmlns="http://www.w3.org/2000/svg"
          style={{ position: 'absolute', top: 0, left: 0 }}
        >
          <rect x="5" y="5" width="1.5" height="6" fill="currentColor" />
          <rect x="9.5" y="5" width="1.5" height="6" fill="currentColor" />
        </svg>
      </span>
    ),
    label: 'Pausing',
  },
  paused: {
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
        <rect x="5" y="5" width="1.5" height="6" fill="currentColor" />
        <rect x="9.5" y="5" width="1.5" height="6" fill="currentColor" />
      </svg>
    ),
    label: 'Paused',
  },
} as const;

export const STATUS_ORDER: (keyof typeof STATUS_CONFIG)[] = [
  'error',
  'pending',
  'pausing',
  'paused',
  'waiting',
  'running',
  'completed',
];

export const groupCalculationsByStatus = (
  calculations: CalculationSummary[]
) => {
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

  return STATUS_ORDER.map(status => ({
    status,
    config: STATUS_CONFIG[status],
    calculations: grouped[status] || [],
  })).filter(group => group.calculations.length > 0);
};
