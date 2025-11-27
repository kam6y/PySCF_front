// src/preload.ts

import { contextBridge, ipcRenderer } from 'electron';

// URLパラメータからポート番号を取得（IPC不要）
const urlParams = new URLSearchParams(window.location.search);
const flaskPortParam = urlParams.get('flask_port');
const flaskPort = flaskPortParam ? parseInt(flaskPortParam, 10) : null;

console.log(`[Preload] Flask port from URL: ${flaskPort}`);

// 検証: ポート番号が有効な範囲かチェック
const isValidPort = flaskPort && flaskPort > 0 && flaskPort < 65536;
if (!isValidPort) {
  console.error(`[Preload] Invalid Flask port: ${flaskPort}`);
}

// Expose protected methods that allow the renderer process to use
// the ipcRenderer without exposing the entire object
contextBridge.exposeInMainWorld('electronAPI', {
  // URLパラメータから取得したポート番号を公開
  flaskPort: isValidPort ? flaskPort : null,

  // Electron API methods
  openExternalUrl: (url: string) =>
    ipcRenderer.invoke('open-external-url', url),
  showAboutDialog: () => ipcRenderer.invoke('show-about-dialog'),
  selectFolder: () =>
    ipcRenderer.invoke('dialog:select-folder') as Promise<{
      canceled: boolean;
      filePath: string | null;
      error?: string;
    }>,
  getPlatform: () => ipcRenderer.invoke('get-platform'),
  isFullScreen: () => ipcRenderer.invoke('get-fullscreen'),
  onFullScreenChange: (callback: (isFullScreen: boolean) => void) => {
    const handler = (_event: any, isFullScreen: boolean) =>
      callback(isFullScreen);
    ipcRenderer.on('fullscreen-changed', handler);
    // Return a cleanup function
    return () => {
      ipcRenderer.removeListener('fullscreen-changed', handler);
    };
  },
});
