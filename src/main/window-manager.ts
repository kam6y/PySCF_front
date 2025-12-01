import { BrowserWindow, app, screen } from 'electron';
import path from 'node:path';

let mainWindow: BrowserWindow | null = null;

export const getMainWindow = (): BrowserWindow | null => {
  return mainWindow;
};

export const setMainWindow = (window: BrowserWindow | null): void => {
  mainWindow = window;
};

export const createWindow = (
  flaskPort: number,
  authToken: string
): BrowserWindow => {
  // Create the browser window.
  const newWindow = new BrowserWindow({
    width: 1400,
    height: 800,
    minWidth: 1200,
    minHeight: 800,
    titleBarStyle: 'hidden',
    titleBarOverlay: {
      color: 'rgba(0, 0, 0, 0)',
      symbolColor: '#000000',
      height: 36,
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

  // ポート番号をURLパラメータとして渡す（IPC不要の堅牢な方式）
  // dist/index.html
  const htmlPath = path.join(__dirname, 'index.html');
  newWindow.loadFile(htmlPath, {
    query: {
      flask_port: String(flaskPort),
      auth_token: authToken,
    },
  });

  console.log(`[Main] Loading window with Flask port: ${flaskPort}`);

  mainWindow = newWindow;
  return newWindow;
};
