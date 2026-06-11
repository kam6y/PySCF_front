import { BrowserWindow, app } from 'electron';
import path from 'node:path';
import {
  getMainRendererEntry,
  installNavigationGuards,
} from './renderer-entry';
import { hardenSession } from './session-hardening';
import { registerBackendAuthInjection } from './backend-auth-injection';

let mainWindow: BrowserWindow | null = null;

export const getMainWindow = (): BrowserWindow | null => {
  return mainWindow;
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
      devTools: !app.isPackaged,
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

  // SEC-006: Harden the window session with CSP response headers and
  // permission restrictions before any content is loaded.
  hardenSession(newWindow.webContents.session, app.isPackaged);

  // SEC-002: inject the auth token into renderer->backend requests at the
  // network layer so the token never enters the renderer/DOM world. Registered
  // before any content loads; backendPort and authToken are valid here (resolved
  // in main.ts before createWindow is called).
  registerBackendAuthInjection(newWindow.webContents.session, backendPort, authToken);

  const htmlPath = path.join(__dirname, 'index.html');
  const rendererEntry = getMainRendererEntry({
    backendPort,
    htmlPath,
    isPackaged: app.isPackaged,
    rendererUrl: process.env.ELECTRON_RENDERER_URL,
  });

  // --- SEC-001: Lock down renderer navigation ---
  installNavigationGuards(newWindow.webContents, rendererEntry);

  if (rendererEntry.type === 'url') {
    newWindow.loadURL(rendererEntry.url);
  } else {
    newWindow.loadFile(rendererEntry.path, {
      query: rendererEntry.query,
    });
  }

  console.log(`[Main] Loading window with backend port: ${backendPort}`);

  mainWindow = newWindow;
  return newWindow;
};
