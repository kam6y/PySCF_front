// src/main.ts

import path from 'node:path';
import { BrowserWindow, app, ipcMain, dialog, shell, Menu } from 'electron';
import { spawn, ChildProcess } from 'child_process';
import http from 'node:http';
import { execFile } from 'child_process';
import { promisify } from 'util';
import fs from 'fs';

let mainWindow: BrowserWindow | null = null;
let pythonProcess: ChildProcess | null = null;
let flaskPort: number | null = null;
let isQuitting = false;
let isCreatingWindow = false;
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
 * Create a clean environment for the Python process using a whitelist approach.
 * This prevents user's local Python environment (pyenv, conda, venv) from interfering.
 *
 * @param condaBinDir - The bin directory of the conda environment to use
 * @param serverPort - The port number for the Flask server
 * @returns A clean environment object
 */
const createCleanEnvironment = (
  condaBinDir: string,
  serverPort: number
): Record<string, string> => {
  // Whitelist of environment variables to pass through
  const ALLOWED_ENV_VARS = [
    // System basics
    'HOME',
    'USER',
    'TMPDIR',
    'SHELL',
    'TERM',
    'LANG',
    'LC_ALL',
    // Display / GUI
    'DISPLAY',
    'XAUTHORITY',
    // MacOS specific
    '__CF_USER_TEXT_ENCODING',
  ];

  const cleanEnv: Record<string, string> = {};

  // 1. Pass through whitelisted variables
  for (const key of ALLOWED_ENV_VARS) {
    if (process.env[key] !== undefined) {
      cleanEnv[key] = process.env[key]!;
    }
  }

  // 2. Construct a clean PATH
  // We want to keep system paths but exclude any user-land Python paths
  const originalPath = process.env.PATH || '';
  const pathEntries = originalPath.split(':').filter(p => {
    // Exclude common Python environment paths
    if (p.includes('/.pyenv/')) return false;
    if (p.includes('/anaconda')) return false;
    if (p.includes('/miniconda')) return false;
    if (p.includes('virtualenvs')) return false;
    return true;
  });

  // Prepend our conda bin directory to ensure it takes precedence
  cleanEnv.PATH = `${condaBinDir}:${pathEntries.join(':')}`;

  // 3. Set application-specific variables
  cleanEnv.CONDA_DEFAULT_ENV = 'pyscf-env';
  cleanEnv.PYSCF_SERVER_PORT = String(serverPort);

  // Explicitly unset potentially conflicting Python variables
  // (They are already not in the whitelist, but this documents intent)
  // PYTHONPATH, PYTHONHOME, VIRTUAL_ENV are NOT copied.

  // 4. Set environment-specific variables
  if (app.isPackaged) {
    cleanEnv.PYSCF_RESOURCES_PATH = process.resourcesPath;
    cleanEnv.PYSCF_ENV = 'production';
    console.log(`Setting PYSCF_RESOURCES_PATH=${process.resourcesPath}`);
  } else {
    cleanEnv.PYSCF_ENV = 'development';
    console.log('Setting PYSCF_ENV=development');
  }

  // 5. Platform-specific fixes
  if (process.platform === 'darwin' && process.arch === 'arm64') {
    cleanEnv.OPENBLAS_CORETYPE = 'ARMV8';
    if (!cleanEnv.LC_ALL) cleanEnv.LC_ALL = 'C';
  }

  return cleanEnv;
};

/**
 * Python/Flaskサーバーを起動する統一関数
 * 設定ファイルに基づいて開発・本番環境で同一の起動方法を使用
 */
const startPythonServer = async (): Promise<void> => {
  return new Promise(async (resolve, reject) => {
    // 重複実行防止: 既にサーバーが起動している場合はスキップ
    if (pythonProcess && !pythonProcess.killed && flaskPort) {
      console.log('Python server already running, skipping startup');
      resolve();
      return;
    }

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

    // パッケージ環境では初回起動時の初期化（大規模ライブラリロード、バイトコードコンパイル等）に時間がかかるため、
    // 遅延とタイムアウトを大幅に延長。2回目以降はキャッシュにより高速化されるため、実際の待機時間は短くなる
    const initialDelay = app.isPackaged ? 10000 : 2000;
    const healthCheckRetries = app.isPackaged ? 120 : 20; // パッケージ環境: 初回最大130秒, 開発環境: 10秒
    const healthCheckInterval = app.isPackaged ? 1000 : 500; // パッケージ環境: 1秒間隔, 開発環境: 0.5秒間隔

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

      // pyenv環境変数を除外した、conda環境専用の環境変数を作成
      const condaBinDir = path.dirname(pythonExecutablePath);
      const envVars = createCleanEnvironment(condaBinDir, serverPort);

      pythonProcess = spawn(pythonExecutablePath, gunicornArgs, {
        cwd: pythonPath,
        stdio: ['pipe', 'pipe', 'pipe'],
        env: envVars,
      });

      // Gunicorn使用時は事前にポートが決まっているので、少し待ってからヘルスチェック開始

      setTimeout(() => {
        console.log(
          `Starting health check for Gunicorn server on port ${flaskPort}`
        );
        checkServerHealth(flaskPort!, healthCheckRetries, healthCheckInterval)
          .then(resolve)
          .catch(reject);
      }, initialDelay);
    } else {
      // フォールバック: 直接実行（設定でGunicorn無効時のみ）
      console.log('Starting server with direct execution (fallback mode)');

      // pyenv環境変数を除外した、conda環境専用の環境変数を作成
      const condaBinDir = path.dirname(pythonExecutablePath);
      const envVars = createCleanEnvironment(condaBinDir, serverPort);

      pythonProcess = spawn(pythonExecutablePath, [], {
        stdio: ['pipe', 'pipe', 'pipe'],
        env: envVars,
      });

      // フォールバックモードでもヘルスチェックで起動を待つ
      // ポートは既に決定済み(flaskPort)
      setTimeout(() => {
        console.log(
          `Starting health check for direct execution server on port ${flaskPort}`
        );
        checkServerHealth(flaskPort!, healthCheckRetries, healthCheckInterval)
          .then(resolve)
          .catch(reject);
      }, 1000);
    }

    // stdout/stderrのログ出力
    pythonProcess.stdout?.on('data', data => {
      const output = data.toString().trim();
      console.log(`[PYTHON STDOUT] ${output}`);

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
 * Aboutダイアログを表示する関数
 */
const showAboutDialog = (): void => {
  const packageJson = require('../package.json');

  const aboutMessage = `${packageJson.description}

Version: ${packageJson.version}

Key Features:
• 3D molecular structure visualization
• Quantum chemistry calculations (DFT, HF, MP2, CCSD, TDDFT)
• Geometry optimization & vibrational frequency analysis
• AI-powered multi-agent system
• PubChem/SMILES structure retrieval

Tech Stack:
Frontend: Electron, React, TypeScript
Backend: Python, Flask, PySCF, RDKit

License: ${packageJson.license}
Author: ${packageJson.author}

GitHub: ${packageJson.homepage}`;

  dialog.showMessageBox({
    type: 'info',
    title: 'About PySCF_front',
    message: 'PySCF_front',
    detail: aboutMessage,
    buttons: ['OK'],
    icon: app.isPackaged
      ? path.join(
        process.resourcesPath,
        'src',
        'assets',
        'icon',
        'mac',
        'Pyscf_front.icns'
      )
      : path.join(
        __dirname,
        '..',
        'src',
        'assets',
        'icon',
        'mac',
        'Pyscf_front.icns'
      ),
  });
};

/**
 * アプリケーションメニューを作成する関数
 */
const createApplicationMenu = (): void => {
  const isMac = process.platform === 'darwin';

  const template: Electron.MenuItemConstructorOptions[] = [
    // macOSの場合、最初にアプリケーションメニューを追加
    ...(isMac
      ? [
        {
          label: app.name,
          submenu: [
            {
              label: 'About PySCF_front',
              click: showAboutDialog,
            },
            { type: 'separator' as const },
            { role: 'services' as const },
            { type: 'separator' as const },
            { role: 'hide' as const },
            { role: 'hideOthers' as const },
            { role: 'unhide' as const },
            { type: 'separator' as const },
            { role: 'quit' as const },
          ],
        },
      ]
      : []),
    // Editメニュー
    {
      label: 'Edit',
      submenu: [
        { role: 'undo' as const },
        { role: 'redo' as const },
        { type: 'separator' as const },
        { role: 'cut' as const },
        { role: 'copy' as const },
        { role: 'paste' as const },
        ...(isMac
          ? [
            { role: 'pasteAndMatchStyle' as const },
            { role: 'delete' as const },
            { role: 'selectAll' as const },
          ]
          : [
            { role: 'delete' as const },
            { type: 'separator' as const },
            { role: 'selectAll' as const },
          ]),
      ],
    },
    // Viewメニュー
    {
      label: 'View',
      submenu: [
        { role: 'reload' as const },
        { role: 'forceReload' as const },
        { role: 'toggleDevTools' as const },
        { type: 'separator' as const },
        { role: 'resetZoom' as const },
        { role: 'zoomIn' as const },
        { role: 'zoomOut' as const },
        { type: 'separator' as const },
        { role: 'togglefullscreen' as const },
      ],
    },
    // Windowメニュー
    {
      label: 'Window',
      submenu: [
        { role: 'minimize' as const },
        { role: 'zoom' as const },
        ...(isMac
          ? [
            { type: 'separator' as const },
            { role: 'front' as const },
            { type: 'separator' as const },
            { role: 'window' as const },
          ]
          : [{ role: 'close' as const }]),
      ],
    },
    // Helpメニュー
    {
      role: 'help' as const,
      submenu: [
        {
          label: 'Learn More',
          click: async () => {
            const packageJson = require('../package.json');
            await shell.openExternal(packageJson.homepage);
          },
        },
        {
          label: 'Report Issue',
          click: async () => {
            const packageJson = require('../package.json');
            await shell.openExternal(packageJson.bugs.url);
          },
        },
      ],
    },
  ];

  const menu = Menu.buildFromTemplate(template);
  Menu.setApplicationMenu(menu);
};

const createWindow = async () => {
  // ガード条件: 既に作成中または既存ウィンドウがある場合は何もしない
  if (isCreatingWindow || mainWindow) {
    console.log('Window creation already in progress or window already exists');
    if (mainWindow && !mainWindow.isDestroyed()) {
      mainWindow.focus();
    }
    return;
  }

  isCreatingWindow = true;
  console.log('Starting window creation...');

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
  } finally {
    isCreatingWindow = false;
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
      mainWindow?.webContents.devToolsWebContents
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
  mainWindow.on('enter-full-screen', () => {
    mainWindow?.webContents.send('fullscreen-changed', true);
  });

  mainWindow.on('leave-full-screen', () => {
    mainWindow?.webContents.send('fullscreen-changed', false);
  });

  // ウィンドウが閉じられた時の処理
  mainWindow.on('closed', () => {
    console.log('Main window closed');
    mainWindow = null;
  });

  // ポート番号をURLパラメータとして渡す（IPC不要の堅牢な方式）
  const htmlPath = path.join(__dirname, 'index.html');
  mainWindow.loadFile(htmlPath, {
    query: { flask_port: String(flaskPort) },
  });

  console.log(`[Main] Loading window with Flask port: ${flaskPort}`);
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
    if (mainWindow && !mainWindow.isDestroyed()) {
      if (mainWindow.isMinimized()) {
        mainWindow.restore();
      }
      mainWindow.focus();
      mainWindow.show();
    } else {
      console.log('Main window not available, creating new window...');
      createWindow();
    }
  });

  // アプリが準備できたら通常の起動処理
  app.whenReady().then(() => {
    createApplicationMenu();
    createWindow();
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
    createWindow();
  }
});

app.on('before-quit', () => {
  isQuitting = true;
  stopPythonServer();
});

// IPC handler for getting platform information
ipcMain.handle('get-platform', () => process.platform);

// IPC handler for getting fullscreen state
ipcMain.handle('get-fullscreen', () => {
  if (mainWindow && !mainWindow.isDestroyed()) {
    return mainWindow.isFullScreen();
  }
  return false;
});

// IPC handler for opening external URLs in default browser
ipcMain.handle('open-external-url', async (_event, url: string) => {
  try {
    // Security: Only allow http and https protocols
    const parsedUrl = new URL(url);
    if (parsedUrl.protocol !== 'http:' && parsedUrl.protocol !== 'https:') {
      console.warn(`Blocked attempt to open non-http(s) URL: ${url}`);
      return { success: false, error: 'Only HTTP and HTTPS URLs are allowed' };
    }

    console.log(`Opening external URL in default browser: ${url}`);
    await shell.openExternal(url);
    return { success: true };
  } catch (error) {
    console.error(`Failed to open external URL: ${url}`, error);
    return { success: false, error: String(error) };
  }
});

// IPC handler for showing About dialog
ipcMain.handle('show-about-dialog', () => {
  showAboutDialog();
});

// IPC handler for selecting a folder
ipcMain.handle('dialog:select-folder', async () => {
  try {
    const result = await dialog.showOpenDialog(
      mainWindow && !mainWindow.isDestroyed() ? mainWindow : (undefined as any),
      {
        properties: ['openDirectory', 'createDirectory'],
        title: 'Select Calculations Directory',
        buttonLabel: 'Select',
      }
    );

    if (result.canceled || result.filePaths.length === 0) {
      return { canceled: true, filePath: null };
    }

    return { canceled: false, filePath: result.filePaths[0] };
  } catch (error) {
    console.error('Failed to open folder selection dialog:', error);
    return { canceled: true, filePath: null, error: String(error) };
  }
});
