const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type RuntimeConfigModule = {
  getBackendRuntimeConfig: () => {
    backendPort: number | null;
    backendBaseUrl: string;
  };
};

const loadRuntimeConfig = (): RuntimeConfigModule => {
  return require('./runtime-config') as RuntimeConfigModule;
};

const setMockWindow = (backendPort: number | null): void => {
  (globalThis as any).window = {
    electronAPI: {
      backendPort,
    },
  };
};

const testGetBackendRuntimeConfigReadsSingleSource = (): void => {
  setMockWindow(5077);
  const { getBackendRuntimeConfig } = loadRuntimeConfig();

  assert.deepEqual(getBackendRuntimeConfig(), {
    backendPort: 5077,
    backendBaseUrl: 'http://127.0.0.1:5077',
  });
};

const testGetBackendRuntimeConfigFailsExplicitlyWithoutPort = (): void => {
  setMockWindow(null);
  const { getBackendRuntimeConfig } = loadRuntimeConfig();

  assert.deepEqual(getBackendRuntimeConfig(), {
    backendPort: null,
    backendBaseUrl: '',
  });
};

const run = (): void => {
  testGetBackendRuntimeConfigReadsSingleSource();
  testGetBackendRuntimeConfigFailsExplicitlyWithoutPort();
  console.log('runtime config tests passed');
};

run();

export {};
