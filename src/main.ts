// src/main.ts

import path from "node:path";
import { BrowserWindow, app, ipcMain, dialog } from "electron";
import { spawn, ChildProcess } from "child_process";
import http from "node:http";

let mainWindow: BrowserWindow;
let pythonProcess: ChildProcess | null = null;
let flaskPort: number | null = null;
let isQuitting = false;

/**
 * Pythonサーバーのヘルスチェックを行い、起動完了を待つ関数
 * @param port チェック対象のポート番号
 * @param retries リトライ回数
 * @param delay リトライ間隔（ミリ秒）
 * @returns サーバーが正常に応答すれば解決されるPromise
 */
const checkServerHealth = (port: number, retries = 20, delay = 500): Promise<void> => {
  return new Promise((resolve, reject) => {
    let attempts = 0;
    const url = `http://127.0.0.1:${port}/health`;
    const interval = setInterval(() => {
      http.get(url, (res) => {
        if (res.statusCode === 200) {
          clearInterval(interval);
          console.log("Python server is healthy.");
          resolve();
        } else {
          console.log(`Health check failed with status: ${res.statusCode}`);
        }
      }).on('error', (err) => {
        attempts++;
        console.log(`Health check attempt ${attempts} failed...`);
        if (attempts >= retries) {
          clearInterval(interval);
          reject(new Error('Python server health check failed after multiple retries.'));
        }
      });
    }, delay);
  });
};

/**
 * Python/Flaskサーバーを起動する関数
 * app.isPackaged の値に応じて、開発モードとパッケージ後で起動方法を切り替える
 */
const startPythonServer = (): Promise<void> => {
  return new Promise((resolve, reject) => {
    console.log('Starting Python Flask server...');
    
    if (app.isPackaged) {
      const executableName = process.platform === 'win32' ? 'pyscf_front_api.exe' : 'pyscf_front_api';
      const pythonExecutablePath = path.join(process.resourcesPath, 'python_dist', executableName);
      
      console.log(`Executing packaged Python server at: ${pythonExecutablePath}`);
      pythonProcess = spawn(pythonExecutablePath, [], {
        stdio: ['pipe', 'pipe', 'pipe']
      });

    } else {
      const pythonPath = path.join(__dirname, '..', 'src', 'python');
      const appPath = path.join(pythonPath, 'app.py');
      
      console.log('Executing Python server in development mode...');
      pythonProcess = spawn('uv', ['run', 'python', appPath], {
        cwd: pythonPath,
        stdio: ['pipe', 'pipe', 'pipe']
      });
    }

    // stdout/stderrのログ出力
    pythonProcess.stdout?.on('data', (data) => {
      const output = data.toString().trim();
      console.log(`Python server: ${output}`);
      // Pythonサーバーからポート番号を受け取る
      const match = output.match(/FLASK_SERVER_PORT:(\d+)/);
      if (match && match[1]) {
        flaskPort = parseInt(match[1], 10);
        console.log(`Detected Flask server port: ${flaskPort}`);
        // ヘルスチェックを開始してサーバー起動を待つ
        checkServerHealth(flaskPort)
          .then(resolve)
          .catch(reject);
      }
    });

    pythonProcess.stderr?.on('data', (data) => {
      console.error(`Python server error: ${data.toString().trim()}`);
    });

    pythonProcess.on('error', (error) => {
      console.error(`Failed to start Python server process: ${error.message}`);
      reject(error);
    });

    pythonProcess.on('close', (code) => {
      console.log(`Python server exited with code ${code}`);
      pythonProcess = null;
      // ユーザーがアプリを終了しようとしている場合以外での終了は異常とみなす
      if (!isQuitting) {
        dialog.showErrorBox(
          "Backend Process Error",
          `The Python backend process has unexpectedly stopped (exit code: ${code}). Please restart the application.`
        );
        app.quit();
      }
    });
  });
};

/**
 * Python/Flaskサーバーを停止する関数
 */
const stopPythonServer = (): void => {
  if (pythonProcess && !pythonProcess.killed) {
    console.log('Stopping Python Flask server...');
    isQuitting = true;
    pythonProcess.kill('SIGTERM'); 
    
    setTimeout(() => {
      if (pythonProcess && !pythonProcess.killed) {
        console.log('Force killing Python server...');
        pythonProcess.kill('SIGKILL');
      }
    }, 5000);
  }
};

const createWindow = async () => {
  try {
    await startPythonServer();
    if (!flaskPort) {
        throw new Error("Could not determine Flask server port.");
    }
    console.log('Python server started successfully.');
  } catch (error) {
    console.error('Failed to start Python server:', error);
    dialog.showErrorBox("Fatal Error", "Could not start the Python backend. The application will now close.");
    app.quit();
    return;
  }
  
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 800,
    minWidth: 1400,
    minHeight: 600,
    titleBarStyle: 'hidden',
    titleBarOverlay: {
      color: 'rgba(0, 0, 0, 0)',
      symbolColor: '#000000',
      height: 40
    },
    webPreferences: {
      preload: path.join(__dirname, "preload.js"),
      nodeIntegration: false,
      contextIsolation: true,
    },
  });

  if (!app.isPackaged) {
    mainWindow.webContents.openDevTools({ mode: 'detach' });
  }
  
  // レンダラープロセスにFlaskサーバーのポートを通知
  mainWindow.webContents.on('did-finish-load', () => {
    mainWindow.webContents.send('set-flask-port', flaskPort);
  });

  mainWindow.loadFile("dist/index.html");
};

app.whenReady().then(createWindow);

app.on("window-all-closed", () => {
  stopPythonServer();
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow();
  }
});

app.on('before-quit', () => {
  isQuitting = true;
  stopPythonServer();
});

// IPC handler for renderer to get port (fallback)
ipcMain.handle('get-flask-port', () => flaskPort);