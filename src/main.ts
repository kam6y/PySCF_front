// src/main.ts

import path from 'node:path';
import { BrowserWindow, app, ipcMain, dialog } from 'electron';
import { spawn, ChildProcess } from 'child_process';
import http from 'node:http';
import { execFile } from 'child_process';
import { promisify } from 'util';
import fs from 'fs';

let mainWindow: BrowserWindow;
let pythonProcess: ChildProcess | null = null;
let flaskPort: number | null = null;
let isQuitting = false;

const execFilePromise = promisify(execFile);

/**
 * condaコマンドを実行し、成功すれば出力を返す
 * @param command condaコマンドの引数配列
 * @returns 成功時は出力文字列、失敗時はnull
 */
const tryCondaCommand = async (command: string[]): Promise<string | null> => {
  try {
    const { stdout } = await execFilePromise('conda', command, {
      timeout: 5000,
      encoding: 'utf8',
    });
    return stdout.trim();
  } catch (error) {
    console.log(`conda command failed: conda ${command.join(' ')}`);
    return null;
  }
};

/**
 * conda環境のPythonパスを段階的に検出する
 * 1. 環境変数 CONDA_ENV_PATH をチェック
 * 2. conda info --base でベースパスを取得
 * 3. conda info --envs で pyscf-env 環境のパスを取得
 * 4. fallback候補パスを試行
 * @returns conda環境のPythonパス、見つからなければnull
 */
const detectCondaEnvironmentPath = async (): Promise<string | null> => {
  const envName = 'pyscf-env';

  // 1. 環境変数でのパス指定をチェック
  const envPath = process.env.CONDA_ENV_PATH;
  if (envPath) {
    const pythonPath = path.join(envPath, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Using conda environment from CONDA_ENV_PATH: ${pythonPath}`);
      return pythonPath;
    } else {
      console.log(
        `CONDA_ENV_PATH specified but python not found at: ${pythonPath}`
      );
    }
  }

  // 2. conda info --base でベースパスを取得
  const baseOutput = await tryCondaCommand(['info', '--base']);
  if (baseOutput) {
    const basePath = baseOutput.trim();
    const pythonPath = path.join(basePath, 'envs', envName, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Found conda environment using base path: ${pythonPath}`);
      return pythonPath;
    }
  }

  // 3. conda info --envs で環境一覧から pyscf-env を検索
  const envsOutput = await tryCondaCommand(['info', '--envs']);
  if (envsOutput) {
    const lines = envsOutput.split('\n');
    for (const line of lines) {
      if (line.includes(envName)) {
        const parts = line.split(/\s+/);
        if (parts.length >= 2) {
          const envPath = parts[parts.length - 1];
          const pythonPath = path.join(envPath, 'bin', 'python');
          if (fs.existsSync(pythonPath)) {
            console.log(
              `Found conda environment from envs list: ${pythonPath}`
            );
            return pythonPath;
          }
        }
      }
    }
  }

  // 4. fallback: 一般的なcondaインストール場所を試行
  const homeDir = process.env.HOME || process.env.USERPROFILE || '';
  const fallbackPaths = ['miniforge3', 'miniconda3', 'anaconda3', 'mambaforge'];

  for (const condaDir of fallbackPaths) {
    const pythonPath = path.join(
      homeDir,
      condaDir,
      'envs',
      envName,
      'bin',
      'python'
    );
    if (fs.existsSync(pythonPath)) {
      console.log(`Found conda environment using fallback path: ${pythonPath}`);
      return pythonPath;
    }
  }

  console.log(
    `conda environment '${envName}' not found in any expected locations`
  );
  return null;
};

/**
 * Python環境のパスを階層的に検出する統合関数
 * 1. 同梱conda環境 (パッケージ時優先)
 * 2. 開発時conda環境検出
 * 3. PyInstaller実行可能ファイル (フォールバック)
 * @returns Python環境のパス、見つからなければnull
 */
const detectPythonEnvironmentPath = async (): Promise<string | null> => {
  // 1. パッケージ時：同梱conda環境を優先チェック
  if (app.isPackaged) {
    const bundledCondaPath = path.join(
      process.resourcesPath,
      'conda_env',
      'bin',
      'python'
    );
    if (fs.existsSync(bundledCondaPath)) {
      console.log(`Using bundled conda environment: ${bundledCondaPath}`);
      return bundledCondaPath;
    } else {
      console.log(
        `Bundled conda environment not found at: ${bundledCondaPath}`
      );
    }
  }

  // 2. 開発時またはパッケージ時フォールバック：conda環境検出
  console.log('Attempting to detect conda environment...');
  const condaPath = await detectCondaEnvironmentPath();
  if (condaPath) {
    return condaPath;
  }

  // 3. パッケージ時のフォールバック：PyInstaller実行可能ファイル
  if (app.isPackaged) {
    const executableName =
      process.platform === 'win32'
        ? 'pyscf_front_api.exe'
        : 'pyscf_front_api';
    const pythonExecutablePath = path.join(
      process.resourcesPath,
      'python_dist',
      'pyscf_front_api',
      executableName
    );
    if (fs.existsSync(pythonExecutablePath)) {
      console.log(
        `Using PyInstaller executable as fallback: ${pythonExecutablePath}`
      );
      return pythonExecutablePath;
    }
  }

  console.log('No Python environment found in any expected locations');
  return null;
};

/**
 * Pythonサーバーのヘルスチェックを行い、起動完了を待つ関数
 * @param port チェック対象のポート番号
 * @param retries リトライ回数
 * @param delay リトライ間隔（ミリ秒）
 * @returns サーバーが正常に応答すれば解決されるPromise
 */
const checkServerHealth = (
  port: number,
  retries = 20,
  delay = 500
): Promise<void> => {
  return new Promise((resolve, reject) => {
    let attempts = 0;
    const url = `http://127.0.0.1:${port}/health`;
    const interval = setInterval(() => {
      http
        .get(url, res => {
          if (res.statusCode === 200) {
            clearInterval(interval);
            console.log('Python server is healthy.');
            resolve();
          } else {
            console.log(`Health check failed with status: ${res.statusCode}`);
          }
        })
        .on('error', _err => {
          attempts++;
          console.log(
            `Health check attempt ${attempts}/${retries} failed for ${url}`
          );
          if (attempts >= retries) {
            clearInterval(interval);
            const diagnosticMessage = app.isPackaged
              ? `Python backend failed to start after ${retries} attempts.\n\nDiagnostic information:\n- Port: ${port}\n- Health endpoint: ${url}\n\nThis may indicate:\n1. Bundled Python environment is corrupted\n2. Port ${port} is blocked by firewall\n3. Python dependencies are missing\n\nPlease report this issue with the console output.`
              : `Python backend failed to start after ${retries} attempts.\n\nDiagnostic information:\n- Port: ${port}\n- Health endpoint: ${url}\n- Environment: Development mode\n\nTroubleshooting steps:\n1. Check if conda environment 'pyscf-env' is activated\n2. Verify all dependencies are installed: conda env create -f .github/environment.yml\n3. Test the Flask server manually: cd src/python && python app.py\n4. Check if port ${port} is available\n\nFor more details, see CLAUDE.md`;
            reject(new Error(diagnosticMessage));
          }
        });
    }, delay);
  });
};

/**
 * Python/Flaskサーバーを起動する関数
 * app.isPackaged の値に応じて、開発モードとパッケージ後で起動方法を切り替える
 */
const startPythonServer = async (): Promise<void> => {
  return new Promise(async (resolve, reject) => {
    console.log('Starting Python Flask server...');

    // 統一されたPython環境検出
    const pythonExecutablePath = await detectPythonEnvironmentPath();
    if (!pythonExecutablePath) {
      const errorMessage = app.isPackaged
        ? `Python environment not found.\n\nThis appears to be a packaging issue. Please report this as a bug.\nThe application requires both:\n1. Bundled conda environment at: ${path.join(
            process.resourcesPath,
            'conda_env'
          )}\n2. PyInstaller executable at: ${path.join(
            process.resourcesPath,
            'python_dist',
            'pyscf_front_api'
          )}`
        : `Python environment not found.\n\nSetup options:\n1. Set environment variable: export CONDA_ENV_PATH=/path/to/your/pyscf-env\n2. Or ensure conda is available and run:\n   conda env create -f .github/environment.yml\n   conda activate pyscf-env\n\nFor more details, see CLAUDE.md`;
      console.error(errorMessage);
      reject(new Error(errorMessage));
      return;
    }

    console.log(`Using Python environment: ${pythonExecutablePath}`);

    // パッケージ時：実行ファイルを直接実行
    if (app.isPackaged) {
      pythonProcess = spawn(pythonExecutablePath, [], {
        stdio: ['pipe', 'pipe', 'pipe'],
      });
    } else {
      // 開発時：Gunicornを使用してFlask-SocketIOアプリケーションを実行
      const pythonPath = path.join(__dirname, '..', 'src', 'python');
      
      // Find available port for Gunicorn
      let gunicornPort = 5000; // Default fallback
      try {
        const net = require('net');
        const server = net.createServer();
        await new Promise((resolve, reject) => {
          server.listen(0, '127.0.0.1', () => {
            gunicornPort = (server.address() as any).port;
            server.close(resolve);
          });
          server.on('error', reject);
        });
      } catch (error) {
        console.log('Could not detect free port, using default 5000');
      }

      // Set the port for Electron to use
      flaskPort = gunicornPort;
      console.log(`Using Gunicorn port: ${gunicornPort}`);

      // Start Gunicorn with Flask-SocketIO app for long-running quantum calculations
      pythonProcess = spawn(pythonExecutablePath, [
        '-m', 'gunicorn',
        '--workers', '1',
        '--worker-class', 'eventlet',
        '--bind', `127.0.0.1:${gunicornPort}`,
        '--timeout', '0',        // 0 = disable timeout for long quantum calculations
        '--keep-alive', '30',
        '--access-logfile', '-',
        '--log-level', 'info',
        '--preload',
        'app:app'  // Flask application with SocketIO integration
      ], {
        cwd: pythonPath,
        stdio: ['pipe', 'pipe', 'pipe'],
        env: { ...process.env, CONDA_DEFAULT_ENV: 'pyscf-env' },
      });

      // For Gunicorn, start health check immediately since port is predetermined
      setTimeout(() => {
        if (flaskPort) {
          console.log(`Starting health check for Gunicorn server on port ${flaskPort}`);
          checkServerHealth(flaskPort).then(resolve).catch(reject);
        }
      }, 2000); // Give Gunicorn a moment to start up
    }

    // stdout/stderrのログ出力
    pythonProcess.stdout?.on('data', data => {
      const output = data.toString().trim();
      console.log(`Python server: ${output}`);
      
      // パッケージ時のみポート検出を実行（Gunicorn使用時は事前に決定済み）
      if (app.isPackaged) {
        // Pythonサーバーからポート番号を受け取る
        const match = output.match(/FLASK_SERVER_PORT:(\d+)/);
        if (match && match[1]) {
          flaskPort = parseInt(match[1], 10);
          console.log(`Detected Flask server port: ${flaskPort}`);
          // ヘルスチェックを開始してサーバー起動を待つ
          checkServerHealth(flaskPort).then(resolve).catch(reject);
        }
      }
    });

    pythonProcess.stderr?.on('data', data => {
      console.error(`Python server error: ${data.toString().trim()}`);
    });

    pythonProcess.on('error', error => {
      console.error(`Failed to start Python server process: ${error.message}`);
      reject(error);
    });

    pythonProcess.on('close', code => {
      console.log(`Python server exited with code ${code}`);
      pythonProcess = null;
      // ユーザーがアプリを終了しようとしている場合以外での終了は異常とみなす
      if (!isQuitting) {
        dialog.showErrorBox(
          'Backend Process Error',
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
      throw new Error('Could not determine Flask server port.');
    }
    console.log('Python server started successfully.');
  } catch (error) {
    console.error('Failed to start Python server:', error);
    dialog.showErrorBox(
      'Fatal Error',
      'Could not start the Python backend. The application will now close.'
    );
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
      height: 40,
    },
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      nodeIntegration: false,
      contextIsolation: true,
      devTools: true,
    },
  });

  if (!app.isPackaged) {
    mainWindow.webContents.openDevTools({
      mode: 'detach',
      activate: false,
    });

    // DevToolsでAutofillエラーを抑制
    mainWindow.webContents.once('devtools-opened', () => {
      mainWindow.webContents.devToolsWebContents
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

  // レンダラープロセスにFlaskサーバーのポートを通知
  mainWindow.webContents.on('did-finish-load', () => {
    mainWindow.webContents.send('set-flask-port', flaskPort);
  });

  mainWindow.loadFile(path.join(__dirname, 'index.html'));
};

app.whenReady().then(createWindow);

app.on('window-all-closed', () => {
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
