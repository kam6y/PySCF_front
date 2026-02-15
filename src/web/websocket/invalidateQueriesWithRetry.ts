import { QueryClient } from '@tanstack/react-query';

export interface InvalidateQueriesWithRetryOptions {
  queryClient: QueryClient;
  queryKey: unknown[];
  maxRetries?: number; // default: 3
  baseDelayMs?: number; // default: 500
  maxJitterMs?: number; // default: 200
}

/**
 * リトライ機能付きクエリ無効化ヘルパー
 * 指数バックオフとジッターを使用してネットワークエラーから回復
 */
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
      // 初回以外は指数バックオフで待機
      if (attempt > 0) {
        const baseDelay = baseDelayMs * Math.pow(2, attempt - 1);
        const jitter = Math.random() * maxJitterMs;
        const delay = baseDelay + jitter;

        console.log(
          `[UnifiedWebSocket] Retrying query invalidation for ${JSON.stringify(queryKey)} (attempt ${attempt + 1}/${maxRetries}) after ${Math.round(delay)}ms`
        );
        await new Promise(resolve => setTimeout(resolve, delay));
      }

      // クエリを無効化して再フェッチ
      await queryClient.invalidateQueries({
        queryKey,
        refetchType: 'active',
      });

      console.log(
        `[UnifiedWebSocket] Successfully invalidated queries for ${JSON.stringify(queryKey)}`
      );
      return; // 成功したら終了
    } catch (error) {
      const isLastAttempt = attempt === maxRetries - 1;

      // エラーの種類を判定
      const errorMessage = error instanceof Error ? error.message : String(error);
      const isNetworkError =
        errorMessage.includes('ERR_NETWORK_CHANGED') ||
        errorMessage.includes('NetworkError') ||
        errorMessage.includes('Failed to fetch');

      console.error(
        `[UnifiedWebSocket] Query invalidation failed for ${JSON.stringify(queryKey)} (attempt ${attempt + 1}/${maxRetries}):`,
        errorMessage
      );

      // 最後の試行でエラーの場合は例外をスロー
      if (isLastAttempt) {
        throw new Error(
          `Failed to invalidate queries after ${maxRetries} attempts: ${errorMessage}`
        );
      }

      // ネットワークエラー以外の場合は即座に失敗
      if (!isNetworkError) {
        console.warn(
          '[UnifiedWebSocket] Non-network error detected, skipping retry'
        );
        throw error;
      }
    }
  }
}

