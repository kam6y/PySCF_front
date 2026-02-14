// src/web/index.tsx

import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { App } from './App';
import { ApiError } from './api/core';

const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      // ネットワーク復帰時の動作
      networkMode: 'online', // オフライン時は自動的に一時停止
      refetchOnReconnect: true, // ネットワーク復帰時に再フェッチ
      refetchOnWindowFocus: false, // フォーカス時の再フェッチは無効化（パフォーマンス向上）

      // キャッシュとデータの鮮度
      staleTime: 30 * 1000, // 30秒（デフォルト）- 頻繁な再フェッチを防ぐ
      gcTime: 5 * 60 * 1000, // 5分（デフォルト）- メモリに保持する時間

      // リトライ戦略
      retry: (failureCount, error) => {
        // ApiError の場合
        if (error instanceof ApiError) {
          // ネットワークエラー: 3回リトライ
          if (error.isNetworkError) {
            return failureCount < 3;
          }
          // 4xx エラー: リトライしない（クライアント側の問題）
          if (error.status >= 400 && error.status < 500) {
            return false;
          }
          // 5xx エラー: 2回リトライ（サーバー側の問題）
          if (error.status >= 500) {
            return failureCount < 2;
          }
        }
        // その他のエラー: 1回リトライ
        return failureCount < 1;
      },

      // 指数バックオフによるリトライ遅延
      retryDelay: attemptIndex => {
        // 指数バックオフ: 1秒 * 2^attemptIndex + ランダムジッター
        // 1回目: ~1秒, 2回目: ~2秒, 3回目: ~4秒
        const baseDelay = Math.min(1000 * Math.pow(2, attemptIndex), 10000);
        const jitter = Math.random() * 1000; // 0-1秒のランダムジッター
        return baseDelay + jitter;
      },
    },
    mutations: {
      // ミューテーションはデフォルトではリトライしない
      retry: false,
    },
  },
});
const root = createRoot(document.getElementById('root') as Element);

// URLパラメータ経由で取得したポート番号（preloadで設定済み）
const flaskPort = window.electronAPI?.flaskPort;

if (!flaskPort) {
  console.error('[index.tsx] CRITICAL: Flask port not available from preload.');
  console.error('[index.tsx] This indicates a serious initialization failure.');
  // 本番環境では、ここでユーザーにエラーダイアログを表示することを推奨
}

// グローバル変数に保存（WebSocket接続用）
window.flaskPort = flaskPort || 5000;

console.log('[index.tsx] Rendering app with Flask port:', window.flaskPort);

// 即座にレンダリング開始（ポートは既に確定済み）
root.render(
  <StrictMode>
    <QueryClientProvider client={queryClient}>
      <App />
    </QueryClientProvider>
  </StrictMode>
);
