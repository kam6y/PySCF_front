const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type RendererEntryModule = {
  getMainRendererEntry: (params: {
    backendPort: number;
    htmlPath: string;
    isPackaged: boolean;
    rendererUrl?: string;
  }) => {
    type: 'file' | 'url';
    path?: string;
    query?: Record<string, string>;
    url?: string;
  };
  getSplashRendererEntry: (params: {
    htmlPath: string;
    isPackaged: boolean;
    rendererUrl?: string;
  }) => {
    type: 'file' | 'url';
    path?: string;
    url?: string;
  };
};

const loadRendererEntry = (): RendererEntryModule => {
  return require('./renderer-entry') as RendererEntryModule;
};

const testMainRendererUsesDevServerWithBackendPort = (): void => {
  const { getMainRendererEntry } = loadRendererEntry();

  const entry = getMainRendererEntry({
    backendPort: 5060,
    htmlPath: '/app/dist/index.html',
    isPackaged: false,
    rendererUrl: 'http://localhost:5173/',
  });

  assert.deepEqual(entry, {
    type: 'url',
    url: 'http://localhost:5173/?backend_port=5060',
  });
};

const testMainRendererUsesFileWhenPackaged = (): void => {
  const { getMainRendererEntry } = loadRendererEntry();

  const entry = getMainRendererEntry({
    backendPort: 5061,
    htmlPath: '/app/dist/index.html',
    isPackaged: true,
    rendererUrl: 'http://localhost:5173/',
  });

  assert.deepEqual(entry, {
    type: 'file',
    path: '/app/dist/index.html',
    query: {
      backend_port: '5061',
    },
  });
};

const testSplashRendererUsesSplashHtmlOnDevServer = (): void => {
  const { getSplashRendererEntry } = loadRendererEntry();

  const entry = getSplashRendererEntry({
    htmlPath: '/app/dist/splash.html',
    isPackaged: false,
    rendererUrl: 'http://localhost:5173',
  });

  assert.deepEqual(entry, {
    type: 'url',
    url: 'http://localhost:5173/splash.html',
  });
};

const testSplashRendererUsesFileWithoutDevServer = (): void => {
  const { getSplashRendererEntry } = loadRendererEntry();

  const entry = getSplashRendererEntry({
    htmlPath: '/app/dist/splash.html',
    isPackaged: false,
  });

  assert.deepEqual(entry, {
    type: 'file',
    path: '/app/dist/splash.html',
  });
};

const run = (): void => {
  testMainRendererUsesDevServerWithBackendPort();
  testMainRendererUsesFileWhenPackaged();
  testSplashRendererUsesSplashHtmlOnDevServer();
  testSplashRendererUsesFileWithoutDevServer();
  console.log('renderer entry tests passed');
};

run();

export {};
