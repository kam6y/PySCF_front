import { BrowserWindow, app } from 'electron';
import path from 'node:path';

let mainWindow: BrowserWindow | null = null;

export const getMainWindow = (): BrowserWindow | null => {
  return mainWindow;
};

export const setMainWindow = (window: BrowserWindow | null): void => {
  mainWindow = window;
};

export const createWindow = (
  backendPort: number,
  authToken: string
): BrowserWindow => {
  // Create the browser window.
  const newWindow = new BrowserWindow({
    width: 1400,
    height: 800,
    minWidth: 1200,
    minHeight: 800,
    show: false, // 初期状態で非表示
    titleBarStyle: 'hidden',
    titleBarOverlay: {
      color: 'rgba(0, 0, 0, 0)',
      symbolColor: '#000000',
      height: 40,
    },
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'), // dist/preload.js
      nodeIntegration: false,
      contextIsolation: true,
      devTools: true,
    },
  });

  if (!app.isPackaged) {
    newWindow.webContents.openDevTools({
      mode: 'detach',
      activate: false,
    });

    // DevToolsでAutofillエラーを抑制
    newWindow.webContents.once('devtools-opened', () => {
      newWindow?.webContents.devToolsWebContents
        ?.executeJavaScript(
          `
        // Autofill機能を無効化
        try {
          if (window.DevToolsAPI) {
            window.DevToolsAPI.disableAutofill = true;
          }
        } catch (e) {
          // エラーを無視
        }
      `
        )
        .catch(() => {
          // executeJavaScriptのエラーを無視
        });
    });
  }

  // ウィンドウの準備が完了したら表示
  newWindow.once('ready-to-show', () => {
    newWindow.show();
  });

  // 全画面状態の変更をレンダラープロセスに通知
  newWindow.on('enter-full-screen', () => {
    newWindow?.webContents.send('fullscreen-changed', true);
  });

  newWindow.on('leave-full-screen', () => {
    newWindow?.webContents.send('fullscreen-changed', false);
  });

  // ウィンドウが閉じられた時の処理
  newWindow.on('closed', () => {
    console.log('Main window closed');
    mainWindow = null;
  });

  // ポート番号をURLパラメータとして渡す
  // 認証トークンはセキュリティのためIPC経由で送信
  // dist/index.html
  const htmlPath = path.join(__dirname, 'index.html');
  newWindow.loadFile(htmlPath, {
    query: {
      backend_port: String(backendPort),
    },
  });

  // ウィンドウのロード完了後にIPCで認証トークンを送信
  // リロード時にもトークンを送信するため 'on' を使用
  newWindow.webContents.on('did-finish-load', () => {
    newWindow.webContents.send('auth-token', authToken);
    console.log('[Main] Auth token sent via IPC');
  });

  console.log(`[Main] Loading window with backend port: ${backendPort}`);

  mainWindow = newWindow;
  return newWindow;
};
