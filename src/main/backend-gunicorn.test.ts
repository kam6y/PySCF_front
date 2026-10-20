const assert: typeof import('node:assert/strict') = require('node:assert/strict');

import type { ServerConfig } from './config';

type BackendGunicornModule = {
  buildGunicornArgs: (
    serverConfig: ServerConfig,
    serverPort: number
  ) => string[];
};

const loadBackendGunicorn = (): BackendGunicornModule => {
  return require('./backend-gunicorn') as BackendGunicornModule;
};

const createServerConfig = (
  overrides: Partial<ServerConfig['gunicorn']> = {}
): ServerConfig => {
  return {
    server: {
      host: '127.0.0.1',
      port: 5000,
    },
    gunicorn: {
      workers: 1,
      worker_class: 'uvicorn.workers.UvicornWorker',
      timeout: 0,
      keep_alive: 2,
      access_logfile: null,
      log_level: 'info',
      preload_app: false,
      ...overrides,
    },
    production: {
      use_gunicorn: true,
    },
  };
};

const testBuildGunicornArgsUsesBackendPort = (): void => {
  const { buildGunicornArgs } = loadBackendGunicorn();

  const args = buildGunicornArgs(createServerConfig(), 5091);

  assert.deepEqual(args, [
    '-m',
    'gunicorn',
    '--workers',
    '1',
    '--worker-class',
    'uvicorn.workers.UvicornWorker',
    '--bind',
    '127.0.0.1:5091',
    '--timeout',
    '0',
    '--keep-alive',
    '2',
    '--log-level',
    'info',
    'app:app',
  ]);
};

const testBuildGunicornArgsKeepsOptionalSettings = (): void => {
  const { buildGunicornArgs } = loadBackendGunicorn();

  const args = buildGunicornArgs(
    createServerConfig({
      access_logfile: '-',
      preload_app: true,
    }),
    5000
  );

  assert.ok(args.includes('--access-logfile'));
  assert.ok(args.includes('-'));
  assert.ok(args.includes('--preload'));
};

const run = (): void => {
  testBuildGunicornArgsUsesBackendPort();
  testBuildGunicornArgsKeepsOptionalSettings();
  console.log('backend gunicorn tests passed');
};

run();

export {};
