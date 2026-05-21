import { spawn, ChildProcess } from 'child_process';
import http from 'node:http';
import path from 'node:path';
import fs from 'fs';
import { app, dialog } from 'electron';
import {
  detectPythonEnvironmentPath,
  createCleanEnvironment,
} from './python-env';
import { findAvailablePort } from './port-manager';
import { updateSplashStatus } from './splash-window-manager';
import type { ServerConfig } from './config';

const buildDiagnosticMessage = (
  port: number,
  url: string,
  retries: number
): string => {
  return app.isPackaged
    ? `Python backend failed to start after ${retries} attempts.\n\nDiagnostic information:\n- Port: ${port}\n- Health endpoint: ${url}\n\nThis may indicate:\n1. Bundled Python environment is corrupted\n2. Port ${port} is blocked by firewall\n3. Python dependencies are missing\n\nPlease report this issue with the console output.`
    : `Python backend failed to start after ${retries} attempts.\n\nDiagnostic information:\n- Port: ${port}\n- Health endpoint: ${url}\n- Environment: Development mode\n\nTroubleshooting steps:\n1. Check if conda environment 'pyscf-env' is activated\n2. Verify all dependencies are installed: conda env create -f .github/environment.yml\n3. Test the backend manually: cd src/python && python app.py\n4. Check if port ${port} is available\n\nFor more details, see CLAUDE.md`;
};
const resolveServerPort = async (defaultPort: number): Promise<number> => {
  const portRangeEnd = 5100;
  try {
    console.log(
      `Auto-detecting available port in range ${defaultPort}-${portRangeEnd}...`
    );
    const port = await findAvailablePort(defaultPort, portRangeEnd);
    console.log(`✓ Found available port: ${port}`);
    return port;
  } catch (error) {
    console.log(`⚠️  Auto-detection failed: ${error}`);
    console.log(`Attempting to use fallback port: ${defaultPort}`);
    try {
      await findAvailablePort(defaultPort, defaultPort);
      console.log(`✓ Fallback port ${defaultPort} is available`);
      return defaultPort;
    } catch (fallbackError) {
      console.log(
        `✗ CRITICAL: Fallback port ${defaultPort} is also unavailable`
      );
      try {
        console.log(
          `Searching in extended range ${portRangeEnd + 1}-${portRangeEnd + 100}...`
        );
        const port = await findAvailablePort(
          portRangeEnd + 1,
          portRangeEnd + 100
        );
        console.log(`✓ Found port in extended range: ${port}`);
        return port;
      } catch (extendedError) {
        throw new Error(
          `CRITICAL: No available ports found in any range. This may indicate:\n1. Too many services running on localhost\n2. Firewall blocking port access\n3. System resource limitations\n\nTried ranges: ${defaultPort}-${portRangeEnd}, ${portRangeEnd + 1}-${portRangeEnd + 100}`
        );
      }
    }
  }
};
const buildGunicornArgs = (
  serverConfig: ServerConfig,
  serverPort: number
): string[] => {
  const serverSettings = serverConfig.server;
  const gunicornSettings = serverConfig.gunicorn;
  return [
    '-m',
    'gunicorn',
    '--workers',
    String(gunicornSettings.workers),
    '--worker-class',
    gunicornSettings.worker_class,
    '--bind',
    `${serverSettings.host}:${serverPort}`,
    '--timeout',
    String(gunicornSettings.timeout),
    '--keep-alive',
    String(gunicornSettings.keep_alive),
    ...(gunicornSettings.access_logfile !== null &&
    gunicornSettings.access_logfile !== undefined
      ? ['--access-logfile', gunicornSettings.access_logfile]
      : []),
    '--log-level',
    gunicornSettings.log_level,
    ...(gunicornSettings.preload_app === true ? ['--preload'] : []),
    'app:app',
  ];
};

const attachOutputHandlers = (proc: ChildProcess, serverPort: number): void => {
  proc.stdout?.on('data', data => {
    const output = data.toString().trim();
    console.log(`[PYTHON STDOUT] ${output}`);
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
      console.log('✗ CRITICAL: Python/FastAPI backend object error detected');
    }
    if (
      output.includes('ModuleNotFoundError') ||
      output.includes('ImportError')
    ) {
      console.log(`✗ CRITICAL: Python import error detected - ${output}`);
    }
  });

  proc.stderr?.on('data', data => {
    const errorOutput = data.toString().trim();
    console.log(`[PYTHON STDERR] ${errorOutput}`);
    if (errorOutput.includes('ModuleNotFoundError')) {
      console.log(`✗ CRITICAL: Missing Python module - ${errorOutput}`);
    }
    if (errorOutput.includes('ImportError')) {
      console.log(`✗ CRITICAL: Python import error - ${errorOutput}`);
    }
    if (errorOutput.includes('gunicorn')) {
      console.log(`⚠️  Gunicorn-related error - ${errorOutput}`);
    }
    if (errorOutput.toLowerCase().includes('fastapi')) {
      console.log(`⚠️  FastAPI-related error - ${errorOutput}`);
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
};

type ProcessContext = {
  pythonExecutablePath: string;
  pythonPath: string;
  serverPort: number;
};

const attachLifecycleHandlers = (
  proc: ChildProcess,
  context: ProcessContext,
  reject: (err: Error) => void
): void => {
  const { pythonExecutablePath, pythonPath, serverPort } = context;

  proc.on('error', error => {
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

  proc.on('close', (code, signal) => {
    console.log(`=== Python Server Process Terminated ===`);
    console.log(`Exit code: ${code}`);
    console.log(`Signal: ${signal || 'None'}`);
    console.log(`Was quitting: ${isQuitting}`);
    console.log(`Server port: ${serverPort}`);
    console.log(`Python executable: ${pythonExecutablePath}`);
    console.log(`Working directory: ${pythonPath}`);

    pythonProcess = null;

    if (!isQuitting) {
      const errorMessage = `The Python backend process has unexpectedly stopped (exit code: ${code})${signal ? `, signal: ${signal}` : ''}.\n\nDebugging Information:\n• Python executable: ${pythonExecutablePath}\n• Working directory: ${pythonPath}\n• Server port: ${serverPort}\n• Packaged mode: ${app.isPackaged}\n\nPlease check the console output for detailed error messages and restart the application.`;
      console.log(`✗ CRITICAL: Showing error dialog to user`);
      dialog.showErrorBox('Backend Process Error', errorMessage);
      app.quit();
    }
  });
};

let pythonProcess: ChildProcess | null = null;
let flaskPort: number | null = null;
let isQuitting = false;

/**
 * Pythonサーバーのヘルスチェックを行い、起動完了を待つ関数
 * @param port チェック対象のポート番号
 * @param authToken 認証トークン
 * @param retries リトライ回数
 * @param delay リトライ間隔（ミリ秒）
 * @returns サーバーが正常に応答すれば解決されるPromise
 */
export const checkServerHealth = (
  port: number,
  authToken: string,
  retries = 20,
  delay = 500
): Promise<void> => {
  return new Promise((resolve, reject) => {
    let attempts = 0;
    const url = `http://127.0.0.1:${port}/health`;
    const options = {
      headers: {
        'X-Auth-Token': authToken,
      },
    };
    const interval = setInterval(() => {
      http
        .get(url, options, res => {
          if (res.statusCode === 200) {
            clearInterval(interval);
            console.log('Python server is healthy.');
            resolve();
          } else {
            res.resume(); // ストリームリーク防止
            attempts++;
            console.log(
              `Health check attempt ${attempts}/${retries} failed with status: ${res.statusCode}`
            );
            updateSplashStatus(
              'health-check',
              'Waiting for server...',
              attempts
            );
            if (attempts >= retries) {
              clearInterval(interval);
              const diagnosticMessage = buildDiagnosticMessage(
                port,
                url,
                retries
              );
              reject(new Error(diagnosticMessage));
            }
          }
        })
        .on('error', _err => {
          attempts++;
          console.log(
            `Health check attempt ${attempts}/${retries} failed for ${url}`
          );

          // スプラッシュにリトライカウントを表示
          updateSplashStatus('health-check', 'Waiting for server...', attempts);

          if (attempts >= retries) {
            clearInterval(interval);
            const diagnosticMessage = buildDiagnosticMessage(
              port,
              url,
              retries
            );
            reject(new Error(diagnosticMessage));
          }
        });
    }, delay);
  });
};

/**
 * Python/FastAPI backend serverを起動する統一関数
 * 設定ファイルに基づいて開発・本番環境で同一の起動方法を使用
 */
export const startPythonServer = async (
  serverConfig: ServerConfig,
  authToken: string
): Promise<number> => {
  return new Promise(async (resolve, reject) => {
    // 重複実行防止: 既にサーバーが起動している場合はスキップ
    if (pythonProcess && !pythonProcess.killed && flaskPort) {
      console.log('Python server already running, skipping startup');
      resolve(flaskPort);
      return;
    }

    console.log('Starting Python/FastAPI backend...');

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

    /**
     * Port Detection Strategy:
     * Electron is solely responsible for port detection.
     * Python passively receives the port via PYSCF_SERVER_PORT.
     */

    // ポート検出の進捗を通知
    updateSplashStatus('finding-port', 'Finding available port...');

    const defaultPort =
      typeof serverConfig.server.port === 'number'
        ? serverConfig.server.port
        : 5000;

    try {
      flaskPort = await resolveServerPort(defaultPort);
    } catch (e) {
      reject(e);
      return;
    }
    const serverPort = flaskPort;
    console.log(`Using server port: ${serverPort}`);

    // Python作業ディレクトリの決定
    const pythonPath = app.isPackaged
      ? path.join(process.resourcesPath, 'src', 'python')
      : path.join(__dirname, '..', 'src', 'python');

    // パッケージ環境では初回起動時の初期化（大規模ライブラリロード、バイトコードコンパイル等）に時間がかかるため、
    // 遅延とタイムアウトを大幅に延長。2回目以降はキャッシュにより高速化されるため、実際の待機時間は短くなる
    const initialDelay = app.isPackaged ? 10000 : 2000;
    const healthCheckRetries = app.isPackaged ? 120 : 20; // パッケージ環境: 初回最大130秒, 開発環境: 10秒
    const healthCheckInterval = app.isPackaged ? 1000 : 500; // パッケージ環境: 1秒間隔, 開発環境: 0.5秒間隔

    // 本番環境でもGunicornを使用する統一ロジック
    if (!serverConfig.production.use_gunicorn) {
      const errorMessage =
        'Gunicorn is disabled in config. Direct execution fallback has been removed.';
      console.error(errorMessage);
      reject(new Error(errorMessage));
      return;
    }

    console.log('Starting server with Gunicorn (unified mode)');

    // サーバー起動の進捗を通知
    updateSplashStatus('starting-server', 'Starting Python backend...');

    const gunicornArgs = buildGunicornArgs(serverConfig, serverPort);

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
    const envVars = createCleanEnvironment(condaBinDir, serverPort, authToken);

    pythonProcess = spawn(pythonExecutablePath, gunicornArgs, {
      cwd: pythonPath,
      stdio: ['pipe', 'pipe', 'pipe'],
      env: envVars,
    });

    attachOutputHandlers(pythonProcess, serverPort);
    attachLifecycleHandlers(
      pythonProcess,
      { pythonExecutablePath, pythonPath, serverPort },
      reject
    );

    // Gunicorn使用時は事前にポートが決まっているので、少し待ってからヘルスチェック開始
    setTimeout(() => {
      console.log(
        `Starting health check for Gunicorn server on port ${flaskPort}`
      );
      checkServerHealth(
        flaskPort!,
        authToken,
        healthCheckRetries,
        healthCheckInterval
      )
        .then(() => resolve(flaskPort!))
        .catch(reject);
    }, initialDelay);
  });
};

/**
 * Python/FastAPI backend serverを停止する関数
 */
export const stopPythonServer = (): void => {
  if (pythonProcess && !pythonProcess.killed) {
    console.log('Stopping Python/FastAPI backend...');
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
