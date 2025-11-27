// src/web/index.tsx

import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { App } from './App';

const queryClient = new QueryClient();
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
