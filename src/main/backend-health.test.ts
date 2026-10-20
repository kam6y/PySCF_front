const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type BackendHealthModule = {
  buildBackendHealthDiagnosticMessage: (params: {
    isPackaged: boolean;
    port: number;
    url: string;
    retries: number;
  }) => string;
};

const loadBackendHealth = (): BackendHealthModule => {
  return require('./backend-health') as BackendHealthModule;
};

const testBuildBackendHealthDiagnosticMessageUsesFastApiNaming = (): void => {
  const { buildBackendHealthDiagnosticMessage } = loadBackendHealth();

  const message = buildBackendHealthDiagnosticMessage({
    isPackaged: false,
    port: 5055,
    url: 'http://127.0.0.1:5055/health',
    retries: 3,
  });

  assert.match(message, /Python backend failed to start after 3 attempts/);
  assert.match(message, /FastAPI backend/);
  assert.match(message, /Port: 5055/);
  assert.doesNotMatch(message, /Flask/i);
};

const testPackagedDiagnosticMentionsBundledEnvironment = (): void => {
  const { buildBackendHealthDiagnosticMessage } = loadBackendHealth();

  const message = buildBackendHealthDiagnosticMessage({
    isPackaged: true,
    port: 5000,
    url: 'http://127.0.0.1:5000/health',
    retries: 120,
  });

  assert.match(message, /Bundled Python environment/);
  assert.match(message, /Health endpoint: http:\/\/127\.0\.0\.1:5000\/health/);
};

const run = (): void => {
  testBuildBackendHealthDiagnosticMessageUsesFastApiNaming();
  testPackagedDiagnosticMentionsBundledEnvironment();
  console.log('backend health tests passed');
};

run();

export {};
