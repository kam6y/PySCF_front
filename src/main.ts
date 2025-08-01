import path from "node:path";
import { BrowserWindow, app, ipcMain } from "electron";

let mainWindow: BrowserWindow;

app.whenReady().then(() => {
  // アプリの起動イベント発火で BrowserWindow インスタンスを作成
  mainWindow = new BrowserWindow({
    width: 1200,
    height: 800,
    minWidth: 800,
    minHeight: 600,
    titleBarStyle: 'hiddenInset', // macOS style with hidden title bar
    webPreferences: {
      // webpack が出力したプリロードスクリプトを読み込み
      preload: path.join(__dirname, "preload.js"),
      nodeIntegration: false,
      contextIsolation: true,
    },
  });

  // レンダラープロセスをロード
  mainWindow.loadFile("dist/index.html");

  // Window state change events
  mainWindow.on('maximize', () => {
    mainWindow.webContents.send('window-state-changed', true);
  });

  mainWindow.on('unmaximize', () => {
    mainWindow.webContents.send('window-state-changed', false);
  });
});

// IPC handlers for window controls
ipcMain.handle('window-close', () => {
  if (mainWindow) {
    mainWindow.close();
  }
});

ipcMain.handle('window-minimize', () => {
  if (mainWindow) {
    mainWindow.minimize();
  }
});

ipcMain.handle('window-maximize', () => {
  if (mainWindow) {
    if (mainWindow.isMaximized()) {
      mainWindow.unmaximize();
    } else {
      mainWindow.maximize();
    }
  }
});

ipcMain.handle('window-is-maximized', () => {
  return mainWindow ? mainWindow.isMaximized() : false;
});

// すべてのウィンドウが閉じられたらアプリを終了する
app.once("window-all-closed", () => app.quit());