import { contextBridge, ipcRenderer } from 'electron';
import type { SplashAPI, SplashStatusUpdate } from '../types/splash';

/**
 * スプラッシュウィンドウ用のPreloadスクリプト
 * Main ProcessからのIPC通信を受信してRenderer Processに公開する
 */

const splashAPI: SplashAPI = {
  /**
   * 進捗状態の更新を受信
   */
  onUpdateStatus: (callback: (update: SplashStatusUpdate) => void) => {
    const handler = (
      _event: Electron.IpcRendererEvent,
      update: SplashStatusUpdate
    ) => {
      callback(update);
    };
    ipcRenderer.on('splash:update-status', handler);

    // クリーンアップ関数を返す
    return () => {
      ipcRenderer.removeListener('splash:update-status', handler);
    };
  },

  /**
   * エラーメッセージの受信
   */
  onShowError: (callback: (message: string) => void) => {
    const handler = (_event: Electron.IpcRendererEvent, message: string) => {
      callback(message);
    };
    ipcRenderer.on('splash:show-error', handler);

    // クリーンアップ関数を返す
    return () => {
      ipcRenderer.removeListener('splash:show-error', handler);
    };
  },

  /**
   * スプラッシュクローズ指示の受信
   */
  onClose: (callback: () => void) => {
    const handler = () => {
      callback();
    };
    ipcRenderer.on('splash:close', handler);

    // クリーンアップ関数を返す
    return () => {
      ipcRenderer.removeListener('splash:close', handler);
    };
  },
};

// contextBridge経由でRenderer ProcessにsplashAPIを公開
contextBridge.exposeInMainWorld('splashAPI', splashAPI);

console.log('Splash preload script loaded successfully');
