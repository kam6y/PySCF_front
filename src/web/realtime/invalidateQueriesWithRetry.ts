import { QueryClient } from '@tanstack/react-query';

export interface InvalidateQueriesWithRetryOptions {
  queryClient: QueryClient;
  queryKey: unknown[];
  maxRetries?: number;
  baseDelayMs?: number;
  maxJitterMs?: number;
}

export async function invalidateQueriesWithRetry(
  options: InvalidateQueriesWithRetryOptions
): Promise<void> {
  const {
    queryClient,
    queryKey,
    maxRetries = 3,
    baseDelayMs = 500,
    maxJitterMs = 200,
  } = options;

  for (let attempt = 0; attempt < maxRetries; attempt++) {
    try {
      if (attempt > 0) {
        const baseDelay = baseDelayMs * Math.pow(2, attempt - 1);
        const jitter = Math.random() * maxJitterMs;
        const delay = baseDelay + jitter;

        console.log(
          `[CalculationUpdates] Retrying query invalidation for ${JSON.stringify(queryKey)} (attempt ${attempt + 1}/${maxRetries}) after ${Math.round(delay)}ms`
        );
        await new Promise(resolve => setTimeout(resolve, delay));
      }

      await queryClient.invalidateQueries({
        queryKey,
        refetchType: 'active',
      });

      console.log(
        `[CalculationUpdates] Successfully invalidated queries for ${JSON.stringify(queryKey)}`
      );
      return;
    } catch (error) {
      const isLastAttempt = attempt === maxRetries - 1;
      const errorMessage =
        error instanceof Error ? error.message : String(error);
      const isNetworkError =
        errorMessage.includes('ERR_NETWORK_CHANGED') ||
        errorMessage.includes('NetworkError') ||
        errorMessage.includes('Failed to fetch');

      console.error(
        `[CalculationUpdates] Query invalidation failed for ${JSON.stringify(queryKey)} (attempt ${attempt + 1}/${maxRetries}):`,
        errorMessage
      );

      if (isLastAttempt) {
        throw new Error(
          `Failed to invalidate queries after ${maxRetries} attempts: ${errorMessage}`
        );
      }

      if (!isNetworkError) {
        console.warn(
          '[CalculationUpdates] Non-network error detected, skipping retry'
        );
        throw error;
      }
    }
  }
}
