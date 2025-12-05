import { app, BrowserWindow } from 'electron';
import crypto from 'crypto';
import { loadServerConfig } from './main/config';
import { startPythonServer, stopPythonServer } from './main/python-server';
import { createWindow, getMainWindow } from './main/window-manager';
import { createApplicationMenu } from './main/menu';
import { registerIpcHandlers } from './main/ipc';
import {
  createSplashWindow,
  updateSplashStatus,
  showSplashError,
  closeSplashWindow,
  closeSplashWindowWhenReady,
} from './main/splash-window-manager';

let flaskPort: number | null = null;
let authToken: string = '';
let isCreatingWindow = false;
let isQuitting = false;

// グローバル設定を読み込み
let serverConfig: any = null;

const initializeApp = async () => {
  // スプラッシュウィンドウを作成
  createSplashWindow();
  updateSplashStatus('initializing', 'Initializing PySCF Native App...');

  if (!serverConfig) {
    serverConfig = loadServerConfig();
  }

  // ガード条件: 既に作成中または既存ウィンドウがある場合は何もしない
  if (isCreatingWindow || getMainWindow()) {
    console.log('Window creation already in progress or window already exists');
    const mainWindow = getMainWindow();
    if (mainWindow && !mainWindow.isDestroyed()) {
      mainWindow.focus();
    }
    closeSplashWindow();
    return;
  }

  isCreatingWindow = true;

  // Generate a random authentication token if not already generated
  if (!authToken) {
    authToken = crypto.randomBytes(32).toString('hex');
    console.log('Generated authentication token for Python backend');
  }
  console.log('Starting window creation...');

  try {
    // Python環境検出開始を通知
    updateSplashStatus('detecting-env', 'Detecting Python environment...');

    flaskPort = await startPythonServer(serverConfig, authToken);
    if (!flaskPort) {
      throw new Error('Could not determine Flask server port.');
    }
    console.log('Python server started successfully.');

    // メインウィンドウ作成開始を通知
    updateSplashStatus('creating-window', 'Creating main window...');
  } catch (error) {
    console.error('Failed to start Python server:', error);
    const errorMessage = error instanceof Error ? error.message : String(error);

    // スプラッシュにエラーを表示（5秒後に自動終了）
    showSplashError(
      `Could not start the Python backend.\n\nError details: ${errorMessage}`
    );
    return;
  } finally {
    isCreatingWindow = false;
  }

  const mainWindow = createWindow(flaskPort, authToken);

  // メインウィンドウの準備完了を待ってスプラッシュをクローズ
  closeSplashWindowWhenReady(mainWindow);
};

// シングルインスタンスロックを取得
const gotTheLock = app.requestSingleInstanceLock();

if (!gotTheLock) {
  // 既に別のインスタンスが起動している場合は即座に終了
  console.log('Another instance is already running. Exiting...');
  app.quit();
} else {
  // 2つ目のインスタンスが起動しようとした時の処理
  app.on('second-instance', (event, commandLine, workingDirectory) => {
    console.log('Second instance detected. Focusing existing window...');
    // 既存のウィンドウをフォーカス
    const mainWindow = getMainWindow();
    if (mainWindow && !mainWindow.isDestroyed()) {
      if (mainWindow.isMinimized()) {
        mainWindow.restore();
      }
      mainWindow.focus();
      mainWindow.show();
    } else {
      console.log('Main window not available, creating new window...');
      initializeApp();
    }
  });

  // アプリが準備できたら通常の起動処理
  app.whenReady().then(() => {
    createApplicationMenu();
    registerIpcHandlers(getMainWindow);
    initializeApp();
  });
}

app.on('window-all-closed', () => {
  stopPythonServer();
  app.quit();
});

app.on('activate', () => {
  // ウィンドウがない場合のみ新しいウィンドウを作成
  // シングルインスタンスロックにより、複数プロセスの起動は防止される
  if (BrowserWindow.getAllWindows().length === 0) {
    initializeApp();
  }
});

app.on('before-quit', () => {
  isQuitting = true;
  stopPythonServer();
});
