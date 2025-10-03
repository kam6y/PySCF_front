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
      configPath = path.join(
        process.resourcesPath,
        'config',
        'server-config.json'
      );
    }

    if (fs.existsSync(configPath)) {
      const configContent = fs.readFileSync(configPath, 'utf8');
      const config = JSON.parse(configContent);
      console.log(`Loaded server configuration from: ${configPath}`);
      return config;
    } else {
      console.log(
        `Configuration file not found at: ${configPath}. Using defaults.`
      );
    }
  } catch (error) {
    console.log(
      `Failed to load server configuration: ${error}. Using defaults.`
    );
  }

  // デフォルト設定を返す
  return {
    server: {
      host: '127.0.0.1',
      port: {
        default: 5000,
        auto_detect: true,
        range: { start: 5000, end: 5100 },
      },
    },
    gunicorn: {
      workers: 1,
      threads: 4,
      timeout: 0,
      worker_class: 'sync',
      keep_alive: 30,
      log_level: 'info',
    },
    production: { use_gunicorn: true },
    development: { debug: false },
  };
};

// グローバル設定を読み込み
serverConfig = loadServerConfig();

/**
 * 開発時conda環境のPythonパスを簡素化して検出する
 * @returns conda環境のPythonパス、見つからなければnull
 */
const detectCondaEnvironmentPath = async (): Promise<string | null> => {
  const envName = 'pyscf-env';

  // 1. 環境変数での指定をチェック
  const envPath = process.env.CONDA_ENV_PATH;
  if (envPath) {
    const pythonPath = path.join(envPath, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Using conda environment from CONDA_ENV_PATH: ${pythonPath}`);
      return pythonPath;
    }
  }

  // 2. conda info --base で環境パスを取得（1回のみ試行）
  try {
    const { stdout } = await execFilePromise('conda', ['info', '--base'], {
      timeout: 3000,
      encoding: 'utf8',
    });
    const basePath = stdout.trim();
    const pythonPath = path.join(basePath, 'envs', envName, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Found conda environment: ${pythonPath}`);
      return pythonPath;
    }
  } catch (error) {
    console.log('conda command unavailable or failed');
  }

  console.log(`conda environment '${envName}' not found`);
  return null;
};

/**
 * Python環境のパスを簡素化して検出する統合関数
 * 1. パッケージ時: 同梱conda環境のみ
 * 2. 開発時: pyscf-env環境のみ
 * @returns Python環境のパス、見つからなければnull
 */
const detectPythonEnvironmentPath = async (): Promise<string | null> => {
  console.log('=== Python Environment Detection (Simplified) ===');
  console.log(`Running in packaged mode: ${app.isPackaged}`);

  // 1. パッケージ時：同梱conda環境のみ
  if (app.isPackaged) {
    const bundledCondaPath = path.join(
      process.resourcesPath,
      'conda_env',
      'bin',
      'python'
    );
    const bundledGunicornPath = path.join(
      process.resourcesPath,
      'conda_env',
      'bin',
      'gunicorn'
    );

    console.log(`Checking bundled conda environment: ${bundledCondaPath}`);

    if (fs.existsSync(bundledCondaPath) && fs.existsSync(bundledGunicornPath)) {
      console.log(`✓ Using bundled conda environment: ${bundledCondaPath}`);
      return bundledCondaPath;
    } else {
      console.log(`✗ Bundled conda environment incomplete or missing`);
      console.log(`  - Python exists: ${fs.existsSync(bundledCondaPath)}`);
      console.log(`  - Gunicorn exists: ${fs.existsSync(bundledGunicornPath)}`);
      return null;
    }
  }

  // 2. 開発時：pyscf-env環境のみ
  console.log('Detecting development conda environment...');
  const condaPath = await detectCondaEnvironmentPath();
  if (condaPath) {
    console.log(`✓ Using conda environment: ${condaPath}`);
    return condaPath;
  }

  console.log('✗ No Python environment found');
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
const findAvailablePort = async (
  startPort: number,
  endPort: number
): Promise<number> => {
  const net = require('net');
  const attemptedPorts: number[] = [];

  for (let port = startPort; port <= endPort; port++) {
    attemptedPorts.push(port);
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
      if (endPort - startPort <= 5) {
        // 少ない範囲の場合は詳細ログを出力
        console.log(`  Port ${port} is in use, trying next...`);
      }
      continue;
    }
  }

  const rangeSize = endPort - startPort + 1;
  const errorDetails =
    rangeSize <= 10
      ? ` (tried: ${attemptedPorts.join(', ')})`
      : ` (checked ${rangeSize} ports)`;

  throw new Error(
    `No available port found in range ${startPort}-${endPort}${errorDetails}`
  );
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
        ? `Python environment not found.\n\nThe bundled conda environment is missing or incomplete.\nThis appears to be a packaging issue. Please report this as a bug.\n\nRequired location: ${path.join(
            process.resourcesPath,
            'conda_env'
          )}\nRequired components: python, gunicorn, and all dependencies`
        : `Python environment not found.\n\nSetup instructions:\n1. Run automated setup: npm run setup-env\n2. Or set environment variable: export CONDA_ENV_PATH=/path/to/your/pyscf-env\n3. Verify setup: npm run verify-env\n\nFor detailed setup instructions, see CLAUDE.md`;
      console.error(errorMessage);
      reject(new Error(errorMessage));
      return;
    }

    console.log(`Using Python environment: ${pythonExecutablePath}`);

    // 統一されたポート決定ロジック
    let serverPort: number;
    if (serverSettings.port.auto_detect) {
      try {
        console.log(
          `Auto-detecting available port in range ${serverSettings.port.default}-${serverSettings.port.range.end}...`
        );
        serverPort = await findAvailablePort(
          serverSettings.port.default,
          serverSettings.port.range.end
        );
        console.log(`✓ Found available port: ${serverPort}`);
      } catch (error) {
        console.log(`⚠️  Auto-detection failed: ${error}`);
        console.log(
          `Attempting to use fallback port: ${serverSettings.port.default}`
        );

        // フォールバック時もポートの利用可能性をチェック
        try {
          await findAvailablePort(
            serverSettings.port.default,
            serverSettings.port.default
          );
          serverPort = serverSettings.port.default;
          console.log(
            `✓ Fallback port ${serverSettings.port.default} is available`
          );
        } catch (fallbackError) {
          console.log(
            `✗ CRITICAL: Fallback port ${serverSettings.port.default} is also unavailable`
          );
          // 最後の手段として、さらに広い範囲で検索
          try {
            console.log(
              `Searching in extended range ${serverSettings.port.range.end + 1}-${serverSettings.port.range.end + 100}...`
            );
            serverPort = await findAvailablePort(
              serverSettings.port.range.end + 1,
              serverSettings.port.range.end + 100
            );
            console.log(`✓ Found port in extended range: ${serverPort}`);
          } catch (extendedError) {
            const errorMessage = `CRITICAL: No available ports found in any range. This may indicate:\n1. Too many services running on localhost\n2. Firewall blocking port access\n3. System resource limitations\n\nTried ranges: ${serverSettings.port.default}-${serverSettings.port.range.end}, ${serverSettings.port.range.end + 1}-${serverSettings.port.range.end + 100}`;
            console.error(errorMessage);
            reject(new Error(errorMessage));
            return;
          }
        }
      }
    } else {
      serverPort = serverSettings.port.default;
      console.log(`Using configured fixed port: ${serverPort}`);
    }

    flaskPort = serverPort;
    console.log(`Using server port: ${serverPort}`);

    // Python作業ディレクトリの決定
    const pythonPath = app.isPackaged
      ? path.join(process.resourcesPath, 'src', 'python') // パッケージ時は同梱されたPythonソースコード
      : path.join(__dirname, '..', 'src', 'python'); // 開発時

    // 本番環境でもGunicornを使用する統一ロジック
    if (productionSettings.use_gunicorn) {
      console.log('Starting server with Gunicorn (unified mode)');

      const gunicornArgs = [
        '-m',
        'gunicorn',
        '--workers',
        String(gunicornSettings.workers),
        '--threads',
        String(gunicornSettings.threads),
        '--worker-class',
        gunicornSettings.worker_class,
        '--bind',
        `${serverSettings.host}:${serverPort}`,
        '--timeout',
        String(gunicornSettings.timeout),
        '--keep-alive',
        String(gunicornSettings.keep_alive),
        // access_logfileがnullまたは存在しない場合は引数を追加しない
        ...(gunicornSettings.access_logfile !== null &&
        gunicornSettings.access_logfile !== undefined
          ? ['--access-logfile', gunicornSettings.access_logfile]
          : []),
        '--log-level',
        gunicornSettings.log_level,
        ...(gunicornSettings.preload_app === true ? ['--preload'] : []),
        'app:app', // Flask application (not socketio object)
      ];

      console.log(`=== Starting Gunicorn Server ===`);
      console.log(`Python executable: ${pythonExecutablePath}`);
      console.log(`Working directory: ${pythonPath}`);
      console.log(
        `Gunicorn command: ${pythonExecutablePath} ${gunicornArgs.join(' ')}`
      );
      console.log(`Environment variables:`, {
        CONDA_DEFAULT_ENV: 'pyscf-env',
        PATH: process.env.PATH?.split(':')
          .filter(p => p.includes('conda'))
          .slice(0, 3),
      });

      // Verify files exist before starting
      console.log(`File checks:`);
      console.log(
        `• Python executable exists: ${fs.existsSync(pythonExecutablePath)}`
      );
      console.log(`• Working directory exists: ${fs.existsSync(pythonPath)}`);

      const appPyPath = path.join(pythonPath, 'app.py');
      console.log(`• app.py exists: ${fs.existsSync(appPyPath)}`);

      pythonProcess = spawn(pythonExecutablePath, gunicornArgs, {
        cwd: pythonPath,
        stdio: ['pipe', 'pipe', 'pipe'],
        env: {
          ...process.env,
          CONDA_DEFAULT_ENV: 'pyscf-env',
          PYSCF_SERVER_PORT: String(serverPort),
        },
      });

      // Gunicorn使用時は事前にポートが決まっているので、少し待ってからヘルスチェック開始
      // パッケージ環境では初期化に時間がかかるため、遅延とタイムアウトを延長
      const initialDelay = app.isPackaged ? 5000 : 2000;
      const healthCheckRetries = app.isPackaged ? 60 : 20; // パッケージ環境: 30秒, 開発環境: 10秒

      setTimeout(() => {
        console.log(
          `Starting health check for Gunicorn server on port ${flaskPort}`
        );
        checkServerHealth(flaskPort!, healthCheckRetries)
          .then(resolve)
          .catch(reject);
      }, initialDelay);
    } else {
      // フォールバック: 直接実行（設定でGunicorn無効時のみ）
      console.log('Starting server with direct execution (fallback mode)');
      pythonProcess = spawn(pythonExecutablePath, [], {
        stdio: ['pipe', 'pipe', 'pipe'],
        env: {
          ...process.env,
          CONDA_DEFAULT_ENV: 'pyscf-env',
          PYSCF_SERVER_PORT: String(serverPort),
        },
      });
    }

    // stdout/stderrのログ出力
    pythonProcess.stdout?.on('data', data => {
      const output = data.toString().trim();
      console.log(`[PYTHON STDOUT] ${output}`);

      // 直接実行モード（フォールバック）でのポート検出
      if (!productionSettings.use_gunicorn) {
        const match = output.match(/FLASK_SERVER_PORT:(\d+)/);
        if (match && match[1]) {
          flaskPort = parseInt(match[1], 10);
          console.log(
            `✓ Detected Flask server port from direct execution: ${flaskPort}`
          );
          // ヘルスチェックを開始してサーバー起動を待つ
          checkServerHealth(flaskPort).then(resolve).catch(reject);
        }
      }

      // Look for important startup messages
      if (output.includes('Starting gunicorn')) {
        console.log('✓ Gunicorn is starting up...');
      }
      if (output.includes('Listening at:')) {
        console.log('✓ Server is listening for connections');
      }
      if (output.includes('Booting worker')) {
        console.log('✓ Gunicorn worker is starting...');
      }
      if (output.includes('Application object must be callable')) {
        console.log('✗ CRITICAL: Flask application object error detected');
      }
      if (
        output.includes('ModuleNotFoundError') ||
        output.includes('ImportError')
      ) {
        console.log(`✗ CRITICAL: Python import error detected - ${output}`);
      }
    });

    pythonProcess.stderr?.on('data', data => {
      const errorOutput = data.toString().trim();
      console.log(`[PYTHON STDERR] ${errorOutput}`);

      // Analyze stderr for specific error patterns
      if (errorOutput.includes('ModuleNotFoundError')) {
        console.log(`✗ CRITICAL: Missing Python module - ${errorOutput}`);
      }
      if (errorOutput.includes('ImportError')) {
        console.log(`✗ CRITICAL: Python import error - ${errorOutput}`);
      }
      if (errorOutput.includes('gunicorn')) {
        console.log(`⚠️  Gunicorn-related error - ${errorOutput}`);
      }
      if (errorOutput.includes('flask')) {
        console.log(`⚠️  Flask-related error - ${errorOutput}`);
      }
      if (errorOutput.includes('Address already in use')) {
        console.log(`✗ CRITICAL: Port ${serverPort} is already in use`);
      }
      if (
        errorOutput.includes('[CRITICAL]') ||
        errorOutput.includes('CRITICAL')
      ) {
        console.log(`✗ CRITICAL ERROR FROM PYTHON: ${errorOutput}`);
      }
    });

    pythonProcess.on('error', error => {
      console.error(`✗ CRITICAL: Failed to start Python server process`);
      console.error(`Error details: ${error.message}`);
      console.error(`Error code: ${(error as any).code || 'N/A'}`);
      console.error(`Error errno: ${(error as any).errno || 'N/A'}`);
      console.error(`Error syscall: ${(error as any).syscall || 'N/A'}`);
      console.error(`Python executable path: ${pythonExecutablePath}`);
      console.error(`Working directory: ${pythonPath}`);
      console.error(`Environment variables:`, {
        CONDA_DEFAULT_ENV: process.env.CONDA_DEFAULT_ENV,
        PATH: process.env.PATH?.split(':')
          .filter(p => p.includes('conda'))
          .slice(0, 3),
        PYTHONPATH: process.env.PYTHONPATH || 'Not set',
      });
      reject(error);
    });

    pythonProcess.on('close', (code, signal) => {
      console.log(`=== Python Server Process Terminated ===`);
      console.log(`Exit code: ${code}`);
      console.log(`Signal: ${signal || 'None'}`);
      console.log(`Was quitting: ${isQuitting}`);
      console.log(`Server port: ${serverPort}`);
      console.log(`Python executable: ${pythonExecutablePath}`);
      console.log(`Working directory: ${pythonPath}`);

      pythonProcess = null;

      // ユーザーがアプリを終了しようとしている場合以外での終了は異常とみなす
      if (!isQuitting) {
        const errorMessage = `The Python backend process has unexpectedly stopped (exit code: ${code})${signal ? `, signal: ${signal}` : ''}.\n\nDebugging Information:\n• Python executable: ${pythonExecutablePath}\n• Working directory: ${pythonPath}\n• Server port: ${serverPort}\n• Packaged mode: ${app.isPackaged}\n\nPlease check the console output for detailed error messages and restart the application.`;

        console.log(`✗ CRITICAL: Showing error dialog to user`);
        dialog.showErrorBox('Backend Process Error', errorMessage);
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

/**
 * 既存のウィンドウをフォーカス・表示する関数
 */
const focusExistingWindow = (): boolean => {
  const existingWindows = BrowserWindow.getAllWindows();
  if (existingWindows.length > 0) {
    const window = existingWindows[0];
    if (window.isMinimized()) {
      window.restore();
    }
    window.focus();
    window.show();
    console.log('Focused existing window instead of creating new one');
    return true;
  }
  return false;
};

const createWindow = async () => {
  // パッケージ環境でのシングルウィンドウ制限
  if (app.isPackaged) {
    if (focusExistingWindow()) {
      console.log('Single window restriction active - using existing window');
      return;
    }
  }

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
    minWidth: 1200,
    minHeight: 800,
    titleBarStyle: 'hidden',
    titleBarOverlay: {
      color: 'rgba(0, 0, 0, 0)',
      symbolColor: '#000000',
      height: 36,
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
  if (app.isPackaged) {
    // パッケージ環境ではシングルウィンドウ制限により、既存ウィンドウをフォーカスするか何もしない
    if (!focusExistingWindow()) {
      // 既存ウィンドウがない場合のみ新しいウィンドウを作成
      createWindow();
    }
  } else {
    // 開発環境では従来通りの動作
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  }
});

app.on('before-quit', () => {
  isQuitting = true;
  stopPythonServer();
});

// IPC handler for renderer to get port (fallback)
ipcMain.handle('get-flask-port', () => flaskPort);
