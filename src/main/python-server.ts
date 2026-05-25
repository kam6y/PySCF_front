import { spawn } from 'child_process';
import path from 'node:path';
import { app } from 'electron';
import { checkBackendHealth } from './backend-health';
import {
  buildGunicornArgs,
  logGunicornStartupDiagnostics,
} from './backend-gunicorn';
import { resolveBackendPort } from './backend-port';
import { BackendProcessController } from './backend-process';
import {
  detectPythonEnvironmentPath,
  createCleanEnvironment,
} from './python-env';
import { updateSplashStatus } from './splash-window-manager';
import type { ServerConfig } from './config';

type StartupTiming = {
  initialDelay: number;
  healthCheckRetries: number;
  healthCheckInterval: number;
};

const backendProcess = new BackendProcessController();
let backendPort: number | null = null;

const getStartupTiming = (): StartupTiming => {
  if (app.isPackaged) {
    return {
      initialDelay: 10000,
      healthCheckRetries: 120,
      healthCheckInterval: 1000,
    };
  }

  return {
    initialDelay: 2000,
    healthCheckRetries: 20,
    healthCheckInterval: 500,
  };
};

const getPythonBackendPath = (): string => {
  return app.isPackaged
    ? path.join(process.resourcesPath, 'src', 'python')
    : path.join(__dirname, '..', 'src', 'python');
};

const getDefaultBackendPort = (serverConfig: ServerConfig): number => {
  return typeof serverConfig.server.port === 'number'
    ? serverConfig.server.port
    : 5000;
};

const buildMissingPythonEnvironmentMessage = (): string => {
  if (app.isPackaged) {
    return `Python environment not found.\n\nThe bundled conda environment is missing or incomplete.\nThis appears to be a packaging issue. Please report this as a bug.\n\nRequired location: ${path.join(
      process.resourcesPath,
      'conda_env'
    )}\nRequired components: python, gunicorn, and all dependencies`;
  }

  return `Python environment not found.\n\nSetup instructions:\n1. Run automated setup: npm run setup-env\n2. Or set environment variable: export CONDA_ENV_PATH=/path/to/your/pyscf-env\n3. Verify setup: npm run verify-env\n\nFor detailed setup instructions, see CLAUDE.md`;
};

/**
 * Python/FastAPI backend serverを起動する統一関数
 * 設定ファイルに基づいて開発・本番環境で同一の起動方法を使用
 */
export const startPythonServer = async (
  serverConfig: ServerConfig,
  authToken: string
): Promise<number> => {
  if (backendProcess.isRunning() && backendPort) {
    console.log('Python backend already running, skipping startup');
    return backendPort;
  }

  console.log('Starting Python/FastAPI backend...');

  const pythonExecutablePath = await detectPythonEnvironmentPath();
  if (!pythonExecutablePath) {
    const errorMessage = buildMissingPythonEnvironmentMessage();
    console.error(errorMessage);
    throw new Error(errorMessage);
  }

  console.log(`Using Python environment: ${pythonExecutablePath}`);

  updateSplashStatus('finding-port', 'Finding available port...');
  backendPort = await resolveBackendPort(getDefaultBackendPort(serverConfig));
  const serverPort = backendPort;
  console.log(`Using server port: ${serverPort}`);

  if (!serverConfig.production.use_gunicorn) {
    const errorMessage =
      'Gunicorn is disabled in config. Direct execution fallback has been removed.';
    console.error(errorMessage);
    throw new Error(errorMessage);
  }

  const pythonPath = getPythonBackendPath();
  const { initialDelay, healthCheckRetries, healthCheckInterval } =
    getStartupTiming();

  console.log('Starting server with Gunicorn (unified mode)');
  updateSplashStatus('starting-server', 'Starting Python backend...');

  const gunicornArgs = buildGunicornArgs(serverConfig, serverPort);
  logGunicornStartupDiagnostics({
    pythonExecutablePath,
    pythonPath,
    gunicornArgs,
  });

  const condaBinDir = path.dirname(pythonExecutablePath);
  const envVars = createCleanEnvironment(condaBinDir, serverPort, authToken);

  return new Promise((resolve, reject) => {
    const pythonProcess = spawn(pythonExecutablePath, gunicornArgs, {
      cwd: pythonPath,
      stdio: ['pipe', 'pipe', 'pipe'],
      env: envVars,
    });

    backendProcess.attach(
      pythonProcess,
      { pythonExecutablePath, pythonPath, serverPort },
      reject
    );

    setTimeout(() => {
      console.log(
        `Starting health check for Gunicorn server on port ${serverPort}`
      );
      checkBackendHealth({
        port: serverPort,
        authToken,
        retries: healthCheckRetries,
        delay: healthCheckInterval,
        isPackaged: app.isPackaged,
      })
        .then(() => resolve(serverPort))
        .catch(reject);
    }, initialDelay);
  });
};

/**
 * Python/FastAPI backend serverを停止する関数
 */
export const stopPythonServer = (): void => {
  backendProcess.stop();
};
