export const ERROR_CODES = {
  CPU_INSUFFICIENT_SYSTEM: 'CPU_INSUFFICIENT_SYSTEM',
  CPU_INSUFFICIENT_LIMIT: 'CPU_INSUFFICIENT_LIMIT',
  MEMORY_INSUFFICIENT_SYSTEM: 'MEMORY_INSUFFICIENT_SYSTEM',
  MEMORY_INSUFFICIENT_LIMIT: 'MEMORY_INSUFFICIENT_LIMIT',
  RESOURCE_INSUFFICIENT: 'RESOURCE_INSUFFICIENT',
} as const;

export type ResourceErrorCode = (typeof ERROR_CODES)[keyof typeof ERROR_CODES];

/**
 * リソース不足エラーかどうかを判定する。
 * 統一パターン: cpu usage, memory usage, system cpu usage, system memory usage, no active calculations
 */
export function isResourceInsufficientError(message: string): boolean {
  const lower = message.toLowerCase();
  return (
    lower.includes('cpu usage') ||
    lower.includes('memory usage') ||
    lower.includes('system cpu usage') ||
    lower.includes('system memory usage') ||
    lower.includes('no active calculations')
  );
}

/**
 * リソース不足エラーメッセージを詳細な errorCode に分類する。
 */
export function classifyResourceError(message: string): ResourceErrorCode {
  const lower = message.toLowerCase();
  const hasNoActiveCalc = lower.includes('no active calculations');

  if (lower.includes('cpu')) {
    return hasNoActiveCalc
      ? ERROR_CODES.CPU_INSUFFICIENT_SYSTEM
      : ERROR_CODES.CPU_INSUFFICIENT_LIMIT;
  }

  if (lower.includes('memory')) {
    return hasNoActiveCalc
      ? ERROR_CODES.MEMORY_INSUFFICIENT_SYSTEM
      : ERROR_CODES.MEMORY_INSUFFICIENT_LIMIT;
  }

  return ERROR_CODES.RESOURCE_INSUFFICIENT;
}
