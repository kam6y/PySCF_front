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
let serverConfig: any = null;

const execFilePromise = promisify(execFile);

/**
 * サーバー設定を読み込む
 */
const loadServerConfig = (): any => {
  try {
    // 開発環境では config/ ディレクトリから読み込み
    let configPath = path.join(__dirname, '..', 'config', 'server-config.json');
    
    // パッケージ環境では同梱された設定ファイルを使用
    if (app.isPackaged) {
      configPath = path.join(process.resourcesPath, 'config', 'server-config.json');
      // フォールバック: python_distディレクトリ内の設定
      if (!fs.existsSync(configPath)) {
        configPath = path.join(process.resourcesPath, 'python_dist', 'pyscf_front_api', 'config', 'server-config.json');
      }
    }
    
    if (fs.existsSync(configPath)) {
      const configContent = fs.readFileSync(configPath, 'utf8');
      const config = JSON.parse(configContent);
      console.log(`Loaded server configuration from: ${configPath}`);
      return config;
    } else {
      console.log(`Configuration file not found at: ${configPath}. Using defaults.`);
    }
  } catch (error) {
    console.log(`Failed to load server configuration: ${error}. Using defaults.`);
  }
  
  // デフォルト設定を返す
  return {
    server: { host: '127.0.0.1', port: { default: 5000, auto_detect: true, range: { start: 5000, end: 5100 } } },
    gunicorn: { workers: 1, threads: 4, timeout: 0, worker_class: 'sync', keep_alive: 30, log_level: 'info' },
    production: { use_gunicorn: true },
    development: { debug: false }
  };
};

// グローバル設定を読み込み
serverConfig = loadServerConfig();

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
 * ポート検出を行う統一関数
 */
const findAvailablePort = async (startPort: number, endPort: number): Promise<number> => {
  const net = require('net');
  
  for (let port = startPort; port <= endPort; port++) {
    try {
      await new Promise((resolve, reject) => {
        const server = net.createServer();
        server.listen(port, '127.0.0.1', () => {
          server.close(resolve);
        });
        server.on('error', reject);
      });
      return port;
    } catch (error) {
      // Port is in use, try next
      continue;
    }
  }
  throw new Error(`No available port found in range ${startPort}-${endPort}`);
};

/**
 * Python/Flaskサーバーを起動する統一関数
 * 設定ファイルに基づいて開発・本番環境で同一の起動方法を使用
 */
const startPythonServer = async (): Promise<void> => {
  return new Promise(async (resolve, reject) => {
    console.log('Starting Python Flask server...');
    
    const serverSettings = serverConfig.server;
    const gunicornSettings = serverConfig.gunicorn;
    const productionSettings = serverConfig.production;
    
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

    // 統一されたポート決定ロジック
    let serverPort: number;
    if (serverSettings.port.auto_detect) {
      try {
        serverPort = await findAvailablePort(
          serverSettings.port.default,
          serverSettings.port.range.end
        );
      } catch (error) {
        console.log(`Failed to find available port: ${error}. Using default: ${serverSettings.port.default}`);
        serverPort = serverSettings.port.default;
      }
    } else {
      serverPort = serverSettings.port.default;
    }
    
    flaskPort = serverPort;
    console.log(`Using server port: ${serverPort}`);

    // Python作業ディレクトリの決定
    const pythonPath = app.isPackaged 
      ? path.dirname(pythonExecutablePath)  // パッケージ時はPython実行ファイルと同じディレクトリ
      : path.join(__dirname, '..', 'src', 'python');  // 開発時
    
    // 本番環境でもGunicornを使用する統一ロジック
    if (productionSettings.use_gunicorn) {
      console.log('Starting server with Gunicorn (unified mode)');
      
      const gunicornArgs = [
        '-m', 'gunicorn',
        '--workers', String(gunicornSettings.workers),
        '--threads', String(gunicornSettings.threads),
        '--worker-class', gunicornSettings.worker_class,
        '--bind', `${serverSettings.host}:${serverPort}`,
        '--timeout', String(gunicornSettings.timeout),
        '--keep-alive', String(gunicornSettings.keep_alive),
        '--access-logfile', '-',
        '--log-level', gunicornSettings.log_level,
        '--preload',
        'app:app'  // Flask application (not socketio object)
      ];
      
      console.log(`Gunicorn command: ${pythonExecutablePath} ${gunicornArgs.join(' ')}`);
      
      pythonProcess = spawn(pythonExecutablePath, gunicornArgs, {
        cwd: pythonPath,
        stdio: ['pipe', 'pipe', 'pipe'],
        env: { ...process.env, CONDA_DEFAULT_ENV: 'pyscf-env' },
      });
      
      // Gunicorn使用時は事前にポートが決まっているので、少し待ってからヘルスチェック開始
      setTimeout(() => {
        console.log(`Starting health check for Gunicorn server on port ${flaskPort}`);
        checkServerHealth(flaskPort!).then(resolve).catch(reject);
      }, 2000);
      
    } else {
      // フォールバック: 直接実行（設定でGunicorn無効時のみ）
      console.log('Starting server with direct execution (fallback mode)');
      pythonProcess = spawn(pythonExecutablePath, [], {
        stdio: ['pipe', 'pipe', 'pipe'],
      });
      
      // 直接実行時は実行時ポート検出を使用
      // （startPythonServer関数の残りの部分で処理）
    }

    // stdout/stderrのログ出力
    pythonProcess.stdout?.on('data', data => {
      const output = data.toString().trim();
      console.log(`Python server: ${output}`);
      
      // 直接実行モード（フォールバック）でのポート検出
      if (!productionSettings.use_gunicorn) {
        const match = output.match(/FLASK_SERVER_PORT:(\d+)/);
        if (match && match[1]) {
          flaskPort = parseInt(match[1], 10);
          console.log(`Detected Flask server port from direct execution: ${flaskPort}`);
          // ヘルスチェックを開始してサーバー起動を待つ
          checkServerHealth(flaskPort).then(resolve).catch(reject);
        }
      }
    });

    pythonProcess.stderr?.on('data', data => {
      console.log(`Python server: ${data.toString().trim()}`);
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
