import { BrowserWindow, app } from 'electron';
import path from 'path';
import type { SplashStage, SplashStatusUpdate } from '../types/splash';
import { getSplashRendererEntry } from './renderer-entry';

let splashWindow: BrowserWindow | null = null;

/**
 * スプラッシュウィンドウを作成
 * 600x400pxのフレームなしウィンドウ
 */
export const createSplashWindow = (): void => {
  // 既に存在する場合は何もしない
  if (splashWindow && !splashWindow.isDestroyed()) {
    console.log('Splash window already exists');
    return;
  }

  splashWindow = new BrowserWindow({
    width: 600,
    height: 400,
    frame: false,
    transparent: false,
    resizable: false,
    center: true,
    alwaysOnTop: true,
    backgroundColor: '#ffffff',
    webPreferences: {
      preload: path.join(__dirname, 'splashPreload.js'),
      nodeIntegration: false,
      contextIsolation: true,
    },
  });

  const splashPath = path.join(__dirname, 'splash.html');
  const rendererEntry = getSplashRendererEntry({
    htmlPath: splashPath,
    isPackaged: app.isPackaged,
    rendererUrl: process.env.ELECTRON_RENDERER_URL,
  });

  if (rendererEntry.type === 'url') {
    splashWindow.loadURL(rendererEntry.url);
  } else {
    splashWindow.loadFile(rendererEntry.path);
  }

  // 開発環境でのみDevToolsを自動で開く（オプション）
  if (!app.isPackaged && process.env.DEBUG_SPLASH === 'true') {
    splashWindow.webContents.openDevTools({ mode: 'detach' });
  }

  console.log('Splash window created successfully');
};

/**
 * スプラッシュウィンドウの参照を取得
 */
export const getSplashWindow = (): BrowserWindow | null => {
  return splashWindow;
};

/**
 * 進捗状態を更新
 * @param stage 現在のステージ
 * @param message 表示するメッセージ
 * @param retryCount リトライカウント（オプション）
 */
export const updateSplashStatus = (
  stage: SplashStage,
  message: string,
  retryCount?: number
): void => {
  if (
    !splashWindow ||
    splashWindow.isDestroyed() ||
    !splashWindow.webContents
  ) {
    console.warn('Splash window is not available for status update');
    return;
  }

  const update: SplashStatusUpdate = {
    stage,
    message,
    retryCount,
  };

  splashWindow.webContents.send('splash:update-status', update);
  console.log('Splash status updated:', update);
};

/**
 * エラーメッセージを表示
 * @param message エラーメッセージ
 */
export const showSplashError = (message: string): void => {
  if (
    !splashWindow ||
    splashWindow.isDestroyed() ||
    !splashWindow.webContents
  ) {
    console.warn('Splash window is not available for error display');
    return;
  }

  splashWindow.webContents.send('splash:show-error', message);
  console.log('Splash error displayed:', message);

  // 5秒後にアプリを終了
  setTimeout(() => {
    console.log('Closing application after error...');
    app.quit();
  }, 5000);
};

/**
 * スプラッシュウィンドウをクローズ
 * フェードアウトアニメーション後にクローズ
 */
export const closeSplashWindow = (): void => {
  if (
    !splashWindow ||
    splashWindow.isDestroyed() ||
    !splashWindow.webContents
  ) {
    console.warn('Splash window is already closed');
    return;
  }

  console.log('Closing splash window...');

  // フェードアウトアニメーション指示を送信
  splashWindow.webContents.send('splash:close');

  // アニメーション完了後にウィンドウをクローズ（300ms + 余裕）
  setTimeout(() => {
    if (splashWindow && !splashWindow.isDestroyed()) {
      splashWindow.close();
      splashWindow = null;
      console.log('Splash window closed');
    }
  }, 500);
};

/**
 * メインウィンドウの準備完了を待ってスプラッシュをクローズ
 * @param mainWindow メインウィンドウのインスタンス
 */
export const closeSplashWindowWhenReady = (mainWindow: BrowserWindow): void => {
  mainWindow.once('ready-to-show', () => {
    setTimeout(() => {
      closeSplashWindow();
    }, 300);
  });
};
