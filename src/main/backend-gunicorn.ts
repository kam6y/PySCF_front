import fs from 'fs';
import path from 'node:path';
import type { ServerConfig } from './config';

export type GunicornStartupDiagnostics = {
  pythonExecutablePath: string;
  pythonPath: string;
  gunicornArgs: string[];
};

export const buildGunicornArgs = (
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

export const logGunicornStartupDiagnostics = ({
  pythonExecutablePath,
  pythonPath,
  gunicornArgs,
}: GunicornStartupDiagnostics): void => {
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

  console.log(`File checks:`);
  console.log(
    `• Python executable exists: ${fs.existsSync(pythonExecutablePath)}`
  );
  console.log(`• Working directory exists: ${fs.existsSync(pythonPath)}`);
  const appPyPath = path.join(pythonPath, 'app.py');
  console.log(`• app.py exists: ${fs.existsSync(appPyPath)}`);
};
