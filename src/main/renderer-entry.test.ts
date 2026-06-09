const assert: typeof import('node:assert/strict') = require('node:assert/strict');

import type { RendererEntry } from './renderer-entry';

type RendererEntryType = RendererEntry;

type NavigationGuardTargetType = {
  on(
    event: 'will-navigate' | 'will-redirect',
    listener: (
      event: { preventDefault(): void; url?: string },
      url: string
    ) => void
  ): void;
  setWindowOpenHandler(
    handler: (details: { url: string }) => { action: 'deny' | 'allow' }
  ): void;
};

type RendererEntryModule = {
  getMainRendererEntry: (params: {
    backendPort: number;
    htmlPath: string;
    isPackaged: boolean;
    rendererUrl?: string;
  }) => RendererEntryType;
  getSplashRendererEntry: (params: {
    htmlPath: string;
    isPackaged: boolean;
    rendererUrl?: string;
  }) => RendererEntryType;
  isAllowedDevRendererUrl: (rendererUrl: string) => boolean;
  isAllowedNavigation: (
    targetUrl: string,
    rendererEntry: RendererEntryType
  ) => boolean;
  installNavigationGuards: (
    webContents: NavigationGuardTargetType,
    rendererEntry: RendererEntryType,
    deps: { openExternal: (url: string) => Promise<void> }
  ) => void;
};

const loadRendererEntry = (): RendererEntryModule => {
  return require('./renderer-entry') as RendererEntryModule;
};

// ============================================================
// Existing tests (preserved)
// ============================================================

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

// ============================================================
// SEC-002: isAllowedDevRendererUrl tests
// ============================================================

const testAllowedDevUrl_localhostHttp = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(isAllowedDevRendererUrl('http://localhost:5173'), true);
  assert.equal(isAllowedDevRendererUrl('http://localhost:5173/'), true);
  assert.equal(isAllowedDevRendererUrl('http://localhost:3000'), true);
};

const testAllowedDevUrl_127001Http = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(isAllowedDevRendererUrl('http://127.0.0.1:5173'), true);
  assert.equal(isAllowedDevRendererUrl('http://127.0.0.1:8080/'), true);
};

const testDeniedDevUrl_https = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(isAllowedDevRendererUrl('https://localhost:5173'), false);
  assert.equal(isAllowedDevRendererUrl('https://127.0.0.1:5173'), false);
};

const testDeniedDevUrl_remoteHost = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(isAllowedDevRendererUrl('http://evil.com:5173'), false);
  assert.equal(isAllowedDevRendererUrl('http://192.168.1.1:5173'), false);
  assert.equal(isAllowedDevRendererUrl('http://0.0.0.0:5173'), false);
};

const testDeniedDevUrl_credentialsInUrl = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(
    isAllowedDevRendererUrl('http://user:pass@localhost:5173'),
    false
  );
};

const testDeniedDevUrl_nonHttpScheme = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(isAllowedDevRendererUrl('file:///app/dist/index.html'), false);
  assert.equal(isAllowedDevRendererUrl('ftp://localhost:21'), false);
  assert.equal(isAllowedDevRendererUrl('javascript:alert(1)'), false);
};

const testDeniedDevUrl_portless = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  // No port → would silently target port 80 — reject
  assert.equal(isAllowedDevRendererUrl('http://localhost'), false);
  assert.equal(isAllowedDevRendererUrl('http://localhost/'), false);
  assert.equal(isAllowedDevRendererUrl('http://127.0.0.1'), false);
  // Explicit :80 is normalized to '' by WHATWG URL parser — also rejected
  assert.equal(isAllowedDevRendererUrl('http://localhost:80'), false);
};

const testDeniedDevUrl_port0 = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  // Port 0 = ephemeral/invalid target — reject
  assert.equal(isAllowedDevRendererUrl('http://localhost:0'), false);
};

const testDeniedDevUrl_malformed = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  assert.equal(isAllowedDevRendererUrl('not-a-url'), false);
  assert.equal(isAllowedDevRendererUrl(''), false);
  assert.equal(isAllowedDevRendererUrl('//localhost:5173'), false);
};

// SEC-002: Integration — invalid dev URL falls back to file entry

const testMainRendererFallsBackOnInvalidDevUrl = (): void => {
  const { getMainRendererEntry } = loadRendererEntry();

  const entry = getMainRendererEntry({
    backendPort: 5060,
    htmlPath: '/app/dist/index.html',
    isPackaged: false,
    rendererUrl: 'http://evil.com:5173/',
  });

  // Should fall back to file entry, not produce a url entry
  assert.equal(entry.type, 'file');
  if (entry.type === 'file') {
    assert.equal(entry.path, '/app/dist/index.html');
  }
};

const testSplashRendererFallsBackOnInvalidDevUrl = (): void => {
  const { getSplashRendererEntry } = loadRendererEntry();

  const entry = getSplashRendererEntry({
    htmlPath: '/app/dist/splash.html',
    isPackaged: false,
    rendererUrl: 'https://evil.com:5173',
  });

  assert.equal(entry.type, 'file');
  if (entry.type === 'file') {
    assert.equal(entry.path, '/app/dist/splash.html');
  }
};

// ============================================================
// SEC-001: isAllowedNavigation tests
// ============================================================

const testNavAllowed_devModeSameOrigin = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'url' as const,
    url: 'http://localhost:5173/?backend_port=5060',
  };

  // Same origin, different path — allowed
  assert.equal(
    isAllowedNavigation('http://localhost:5173/other-page', entry),
    true
  );
  // Same origin, same root — allowed
  assert.equal(isAllowedNavigation('http://localhost:5173/', entry), true);
  // Same origin with query — allowed
  assert.equal(
    isAllowedNavigation('http://localhost:5173/?foo=bar', entry),
    true
  );
};

const testNavDenied_devModeDifferentOrigin = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'url' as const,
    url: 'http://localhost:5173/?backend_port=5060',
  };

  // Different host
  assert.equal(isAllowedNavigation('http://evil.com/steal', entry), false);
  // Different port
  assert.equal(isAllowedNavigation('http://localhost:9999/', entry), false);
  // Different protocol
  assert.equal(isAllowedNavigation('https://localhost:5173/', entry), false);
};

const testNavAllowed_packagedModeSameDir = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  // The loaded entry itself — must never lock out the packaged app
  assert.equal(isAllowedNavigation('file:///app/dist/index.html', entry), true);
  // File in same directory
  assert.equal(isAllowedNavigation('file:///app/dist/other.html', entry), true);
  // File in subdirectory
  assert.equal(
    isAllowedNavigation('file:///app/dist/sub/page.html', entry),
    true
  );
};

const testNavDenied_packagedModeHttpUrl = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  // HTTP URL in packaged mode — denied
  assert.equal(isAllowedNavigation('http://evil.com/steal', entry), false);
  assert.equal(isAllowedNavigation('https://phishing.com/', entry), false);
};

const testNavDenied_packagedModeTraversalAttempt = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  // Path traversal attempt — denied
  assert.equal(
    isAllowedNavigation('file:///app/dist/../secret/keys.txt', entry),
    false
  );
  // Completely different path
  assert.equal(isAllowedNavigation('file:///etc/passwd', entry), false);
};

const testNavDenied_malformedUrl = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = { type: 'url' as const, url: 'http://localhost:5173/' };

  assert.equal(isAllowedNavigation('not-a-url', entry), false);
  assert.equal(isAllowedNavigation('', entry), false);
};

// --- F1: about:blank / about:srcdoc ---

const testNavAllowed_aboutBlank = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const fileEntry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };
  const urlEntry = { type: 'url' as const, url: 'http://localhost:5173/' };

  // about:blank must be allowed in both modes
  assert.equal(isAllowedNavigation('about:blank', fileEntry), true);
  assert.equal(isAllowedNavigation('about:blank', urlEntry), true);
  // about:srcdoc must be allowed in both modes
  assert.equal(isAllowedNavigation('about:srcdoc', fileEntry), true);
  assert.equal(isAllowedNavigation('about:srcdoc', urlEntry), true);
};

const testNavDenied_aboutOther = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = { type: 'url' as const, url: 'http://localhost:5173/' };

  // about: URLs other than blank/srcdoc must be denied
  assert.equal(isAllowedNavigation('about:invalid', entry), false);
  assert.equal(isAllowedNavigation('about:config', entry), false);
};

// --- F5: sibling-prefix path must be denied ---

const testNavDenied_packagedModeSiblingPrefix = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  // /app/distractor/ shares the prefix "/app/dist" but is NOT inside /app/dist/
  assert.equal(
    isAllowedNavigation('file:///app/distractor/x.html', entry),
    false
  );
  // /app/dist-evil/ similar sibling-prefix attack
  assert.equal(
    isAllowedNavigation('file:///app/dist-evil/x.html', entry),
    false
  );
};

// --- F3b regression: file URL with query string must be allowed ---

const testNavAllowed_packagedModeFileWithQuery = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  // Electron appends query params to the loaded file URL;
  // navigation guards pass this URL through isAllowedNavigation.
  assert.equal(
    isAllowedNavigation('file:///app/dist/index.html?backend_port=5060', entry),
    true
  );
};

// ============================================================
// F4: installNavigationGuards wiring tests
// ============================================================

type FakeEvent = { defaultPrevented: boolean; preventDefault: () => void };
type EventHandler = (
  event: { preventDefault(): void; url?: string },
  url: string
) => void;
type WindowOpenHandler = (details: { url: string }) => {
  action: 'deny' | 'allow';
};

const createFakeEvent = (): FakeEvent => {
  const evt: FakeEvent = {
    defaultPrevented: false,
    preventDefault() {
      evt.defaultPrevented = true;
    },
  };
  return evt;
};

const createFakeWebContents = () => {
  const handlers: Partial<
    Record<'will-navigate' | 'will-redirect', EventHandler>
  > = {};
  let windowOpenHandler: WindowOpenHandler | null = null;

  return {
    on(event: 'will-navigate' | 'will-redirect', handler: EventHandler) {
      handlers[event] = handler;
    },
    setWindowOpenHandler(handler: WindowOpenHandler) {
      windowOpenHandler = handler;
    },
    getHandler(
      event: 'will-navigate' | 'will-redirect'
    ): EventHandler | undefined {
      return handlers[event];
    },
    getWindowOpenHandler(): WindowOpenHandler | null {
      return windowOpenHandler;
    },
  };
};

const testGuards_willNavigateBlocksDisallowedUrl = (): void => {
  const { installNavigationGuards } = loadRendererEntry();
  const fake = createFakeWebContents();
  const entry = { type: 'url' as const, url: 'http://localhost:5173/' };
  const openExternalCalls: string[] = [];

  installNavigationGuards(fake, entry, {
    openExternal: async url => {
      openExternalCalls.push(url);
    },
  });

  // Disallowed URL must call preventDefault
  const handler = fake.getHandler('will-navigate');
  assert.ok(handler, 'will-navigate handler must be registered');
  const evt = createFakeEvent();
  handler(evt, 'http://evil.com/steal');
  assert.equal(
    evt.defaultPrevented,
    true,
    'will-navigate should block disallowed URL'
  );

  // Allowed URL must NOT call preventDefault
  const evt2 = createFakeEvent();
  handler(evt2, 'http://localhost:5173/page');
  assert.equal(
    evt2.defaultPrevented,
    false,
    'will-navigate should allow same-origin URL'
  );
};

const testGuards_willRedirectBlocksDisallowedUrl = (): void => {
  const { installNavigationGuards } = loadRendererEntry();
  const fake = createFakeWebContents();
  const entry = { type: 'url' as const, url: 'http://localhost:5173/' };

  installNavigationGuards(fake, entry, {
    openExternal: async () => {},
  });

  const handler = fake.getHandler('will-redirect');
  assert.ok(handler, 'will-redirect handler must be registered');

  // Disallowed redirect must be blocked
  const evt = createFakeEvent();
  handler(evt, 'http://evil.com/redirect-target');
  assert.equal(
    evt.defaultPrevented,
    true,
    'will-redirect should block disallowed URL'
  );

  // Allowed redirect must pass through
  const evt2 = createFakeEvent();
  handler(evt2, 'http://localhost:5173/redirected');
  assert.equal(
    evt2.defaultPrevented,
    false,
    'will-redirect should allow same-origin URL'
  );
};

const testGuards_windowOpenDeniesAllAndOpensExternalForHttp = (): void => {
  const { installNavigationGuards } = loadRendererEntry();
  const fake = createFakeWebContents();
  const entry = { type: 'url' as const, url: 'http://localhost:5173/' };
  const openExternalCalls: string[] = [];

  installNavigationGuards(fake, entry, {
    openExternal: async url => {
      openExternalCalls.push(url);
    },
  });

  const handler = fake.getWindowOpenHandler();
  assert.ok(handler, 'window-open handler must be registered');

  // http URL: deny + call openExternal
  const result1 = handler({ url: 'http://example.com/page' });
  assert.equal(result1.action, 'deny');
  assert.equal(openExternalCalls.length, 1);
  assert.equal(openExternalCalls[0], 'http://example.com/page');

  // https URL: deny + call openExternal
  const result2 = handler({ url: 'https://example.com/secure' });
  assert.equal(result2.action, 'deny');
  assert.equal(openExternalCalls.length, 2);
  assert.equal(openExternalCalls[1], 'https://example.com/secure');

  // file URL: deny + do NOT call openExternal
  const result3 = handler({ url: 'file:///etc/passwd' });
  assert.equal(result3.action, 'deny');
  assert.equal(
    openExternalCalls.length,
    2,
    'openExternal must not be called for file: URLs'
  );

  // about:blank: deny + do NOT call openExternal
  const result4 = handler({ url: 'about:blank' });
  assert.equal(result4.action, 'deny');
  assert.equal(
    openExternalCalls.length,
    2,
    'openExternal must not be called for about: URLs'
  );

  // Malformed URL: deny + do NOT call openExternal
  const result5 = handler({ url: 'not-a-valid-url' });
  assert.equal(result5.action, 'deny');
  assert.equal(
    openExternalCalls.length,
    2,
    'openExternal must not be called for malformed URLs'
  );
};

// --- F2: file-entry variants for will-navigate / will-redirect ---

const testGuards_willNavigateBlocksDisallowedUrl_fileEntry = (): void => {
  const { installNavigationGuards } = loadRendererEntry();
  const fake = createFakeWebContents();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  installNavigationGuards(fake, entry, {
    openExternal: async () => {},
  });

  const handler = fake.getHandler('will-navigate');
  assert.ok(handler, 'will-navigate handler must be registered (file entry)');

  // Disallowed: HTTP URL in packaged/file mode
  const evt1 = createFakeEvent();
  handler(evt1, 'http://evil.com/steal');
  assert.equal(
    evt1.defaultPrevented,
    true,
    'will-navigate should block HTTP URL for file entry'
  );

  // Disallowed: file outside allowed directory
  const evt2 = createFakeEvent();
  handler(evt2, 'file:///etc/passwd');
  assert.equal(
    evt2.defaultPrevented,
    true,
    'will-navigate should block file outside allowed dir'
  );

  // Allowed: file in same directory
  const evt3 = createFakeEvent();
  handler(evt3, 'file:///app/dist/other.html');
  assert.equal(
    evt3.defaultPrevented,
    false,
    'will-navigate should allow same-dir file URL'
  );
};

const testGuards_willRedirectBlocksDisallowedUrl_fileEntry = (): void => {
  const { installNavigationGuards } = loadRendererEntry();
  const fake = createFakeWebContents();
  const entry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };

  installNavigationGuards(fake, entry, {
    openExternal: async () => {},
  });

  const handler = fake.getHandler('will-redirect');
  assert.ok(handler, 'will-redirect handler must be registered (file entry)');

  // Disallowed: HTTP URL in packaged/file mode
  const evt1 = createFakeEvent();
  handler(evt1, 'http://evil.com/redirect-target');
  assert.equal(
    evt1.defaultPrevented,
    true,
    'will-redirect should block HTTP URL for file entry'
  );

  // Disallowed: file outside allowed directory
  const evt2 = createFakeEvent();
  handler(evt2, 'file:///etc/passwd');
  assert.equal(
    evt2.defaultPrevented,
    true,
    'will-redirect should block file outside allowed dir'
  );

  // Allowed: file in same directory
  const evt3 = createFakeEvent();
  handler(evt3, 'file:///app/dist/other.html');
  assert.equal(
    evt3.defaultPrevented,
    false,
    'will-redirect should allow same-dir file URL'
  );
};

// ============================================================
// F3: Boundary tests — fail-closed edge cases
// ============================================================

const testDeniedDevUrl_ipv6Loopback = (): void => {
  const { isAllowedDevRendererUrl } = loadRendererEntry();
  // IPv6 loopback is NOT in ALLOWED_DEV_HOSTNAMES — must be rejected (fail-closed)
  assert.equal(isAllowedDevRendererUrl('http://[::1]:5173'), false);
};

const testNavDenied_urlEntryMalformedUrl = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  // url entry with empty/malformed URL — inner new URL('') throws, catch→deny
  const entry = { type: 'url' as const, url: '' };
  assert.equal(isAllowedNavigation('http://localhost:5173/', entry), false);

  const entry2 = { type: 'url' as const, url: 'not-a-valid-url' };
  assert.equal(isAllowedNavigation('http://localhost:5173/', entry2), false);
};

// --- F1: data: and javascript: URI denial (XSS/navigation-bypass vectors) ---

const testNavDenied_dataUri = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const fileEntry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };
  const urlEntry = { type: 'url' as const, url: 'http://localhost:5173/' };

  // data: URIs must be denied in both file and url modes
  assert.equal(
    isAllowedNavigation('data:text/html,<script>alert(1)</script>', fileEntry),
    false,
    'data: URI must be denied for file entry'
  );
  assert.equal(
    isAllowedNavigation('data:text/html,<script>alert(1)</script>', urlEntry),
    false,
    'data: URI must be denied for url entry'
  );
};

const testNavDenied_javascriptUri = (): void => {
  const { isAllowedNavigation } = loadRendererEntry();
  const fileEntry = {
    type: 'file' as const,
    path: '/app/dist/index.html',
    query: { backend_port: '5060' },
  };
  const urlEntry = { type: 'url' as const, url: 'http://localhost:5173/' };

  // javascript: URIs must be denied in both file and url modes
  assert.equal(
    isAllowedNavigation('javascript:alert(1)', fileEntry),
    false,
    'javascript: URI must be denied for file entry'
  );
  assert.equal(
    isAllowedNavigation('javascript:alert(1)', urlEntry),
    false,
    'javascript: URI must be denied for url entry'
  );
};

// ============================================================
// Runner
// ============================================================

const run = (): void => {
  // Existing tests
  testMainRendererUsesDevServerWithBackendPort();
  testMainRendererUsesFileWhenPackaged();
  testSplashRendererUsesSplashHtmlOnDevServer();
  testSplashRendererUsesFileWithoutDevServer();

  // SEC-002: isAllowedDevRendererUrl
  testAllowedDevUrl_localhostHttp();
  testAllowedDevUrl_127001Http();
  testDeniedDevUrl_https();
  testDeniedDevUrl_remoteHost();
  testDeniedDevUrl_credentialsInUrl();
  testDeniedDevUrl_nonHttpScheme();
  testDeniedDevUrl_portless();
  testDeniedDevUrl_port0();
  testDeniedDevUrl_malformed();
  testMainRendererFallsBackOnInvalidDevUrl();
  testSplashRendererFallsBackOnInvalidDevUrl();

  // SEC-001: isAllowedNavigation
  testNavAllowed_devModeSameOrigin();
  testNavDenied_devModeDifferentOrigin();
  testNavAllowed_packagedModeSameDir();
  testNavDenied_packagedModeHttpUrl();
  testNavDenied_packagedModeTraversalAttempt();
  testNavDenied_malformedUrl();

  // F1: about:blank / about:srcdoc
  testNavAllowed_aboutBlank();
  testNavDenied_aboutOther();

  // F3b regression: file URL with query string
  testNavAllowed_packagedModeFileWithQuery();

  // F5: sibling-prefix guard
  testNavDenied_packagedModeSiblingPrefix();

  // F4: installNavigationGuards wiring
  testGuards_willNavigateBlocksDisallowedUrl();
  testGuards_willRedirectBlocksDisallowedUrl();
  testGuards_windowOpenDeniesAllAndOpensExternalForHttp();

  // F2: file-entry guard variants
  testGuards_willNavigateBlocksDisallowedUrl_fileEntry();
  testGuards_willRedirectBlocksDisallowedUrl_fileEntry();

  // F3: Boundary tests
  testDeniedDevUrl_ipv6Loopback();
  testNavDenied_urlEntryMalformedUrl();

  // F1: data: and javascript: URI denial (XSS/navigation-bypass vectors)
  testNavDenied_dataUri();
  testNavDenied_javascriptUri();

  console.log('renderer entry tests passed (34 tests)');
};

run();

export {};
