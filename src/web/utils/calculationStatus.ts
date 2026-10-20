import type { CalculationStatus } from '../types/api-types';

export function isCalculationEditable(
  status: CalculationStatus | undefined
): boolean {
  return (
    status === undefined ||
    status === 'pending' ||
    status === 'completed' ||
    status === 'error'
  );
}
