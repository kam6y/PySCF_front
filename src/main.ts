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
            reject(
              new Error(
                `Python server health check failed after ${retries} attempts. Check if the conda environment is properly set up and the Python server can start.`
              )
            );
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

      console.log(
        `Executing packaged Python server at: ${pythonExecutablePath}`
      );
      pythonProcess = spawn(pythonExecutablePath, [], {
        stdio: ['pipe', 'pipe', 'pipe'],
      });
    } else {
      const pythonPath = path.join(__dirname, '..', 'src', 'python');
      const appPath = path.join(pythonPath, 'app.py');

      console.log('Executing Python server in development mode...');
      console.log('Detecting conda environment...');

      // Use conda environment (required for development)
      const condaPath = await detectCondaEnvironmentPath();
      if (!condaPath) {
        const errorMessage = `Conda environment 'pyscf-env' not found.\n\nSetup options:\n1. Set environment variable: export CONDA_ENV_PATH=/path/to/your/pyscf-env\n2. Or ensure conda is available and run:\n   conda create -n pyscf-env python=3.12\n   conda activate pyscf-env\n   conda install -c conda-forge pyscf rdkit flask geometric requests flask-cors pydantic gevent threadpoolctl\n   pip install flask-sock flask-pydantic datamodel-code-generator pyinstaller gevent-websocket\n\nFor more details, see README.md`;
        console.error(errorMessage);
        reject(new Error(errorMessage));
        return;
      }

      console.log(`Using conda environment: ${condaPath}`);
      pythonProcess = spawn(condaPath, [appPath], {
        cwd: pythonPath,
        stdio: ['pipe', 'pipe', 'pipe'],
        env: { ...process.env, CONDA_DEFAULT_ENV: 'pyscf-env' },
      });
    }

    // stdout/stderrのログ出力
    pythonProcess.stdout?.on('data', data => {
      const output = data.toString().trim();
      console.log(`Python server: ${output}`);
      // Pythonサーバーからポート番号を受け取る
      const match = output.match(/FLASK_SERVER_PORT:(\d+)/);
      if (match && match[1]) {
        flaskPort = parseInt(match[1], 10);
        console.log(`Detected Flask server port: ${flaskPort}`);
        // ヘルスチェックを開始してサーバー起動を待つ
        checkServerHealth(flaskPort).then(resolve).catch(reject);
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

  mainWindow.loadFile('dist/index.html');
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
