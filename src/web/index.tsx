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
const backendPort = window.electronAPI?.flaskPort;

if (!backendPort) {
  console.error(
    '[index.tsx] CRITICAL: Backend port not available from preload.'
  );

  // ポート未取得時はエラーUIをレンダリングし、アプリを起動しない（fail-fast）
  root.render(
    <StrictMode>
      <div
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          justifyContent: 'center',
          height: '100vh',
          fontFamily: 'system-ui, sans-serif',
          color: '#c0392b',
          textAlign: 'center',
          padding: '2rem',
        }}
      >
        <h1 style={{ fontSize: '1.5rem', marginBottom: '1rem' }}>
          アプリケーションの起動に失敗しました
        </h1>
        <p>
          バックエンドサーバーのポート番号を取得できませんでした。
          <br />
          アプリを再起動してください。
        </p>
      </div>
    </StrictMode>
  );
} else {
  // グローバル変数に保存（WebSocket接続用）
  window.flaskPort = backendPort;

  console.log('[index.tsx] Rendering app with backend port:', window.flaskPort);

  root.render(
    <StrictMode>
      <QueryClientProvider client={queryClient}>
        <App />
      </QueryClientProvider>
    </StrictMode>
  );
}
