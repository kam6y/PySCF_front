const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type PermissionCallback = (granted: boolean) => void;
type PermissionRequestHandler = (
  webContents: unknown,
  permission: string,
  callback: PermissionCallback
) => void;

type HeadersReceivedDetails = {
  responseHeaders?: Record<string, string[]>;
};
type HeadersReceivedCallback = (result: {
  responseHeaders?: Record<string, string[]>;
}) => void;
type HeadersReceivedHandler = (
  details: HeadersReceivedDetails,
  callback: HeadersReceivedCallback
) => void;

type FakeSession = {
  setPermissionRequestHandler: (handler: PermissionRequestHandler) => void;
  webRequest: {
    onHeadersReceived: (handler: HeadersReceivedHandler) => void;
  };
  // Accessors for test assertions
  getPermissionRequestHandler: () => PermissionRequestHandler | null;
  getHeadersReceivedHandler: () => HeadersReceivedHandler | null;
  getPermissionRequestHandlerCallCount: () => number;
  getHeadersReceivedCallCount: () => number;
};

const createFakeSession = (): FakeSession => {
  let permissionHandler: PermissionRequestHandler | null = null;
  let headersHandler: HeadersReceivedHandler | null = null;
  let permissionCallCount = 0;
  let headersCallCount = 0;

  return {
    setPermissionRequestHandler(handler: PermissionRequestHandler) {
      permissionHandler = handler;
      permissionCallCount++;
    },
    webRequest: {
      onHeadersReceived(handler: HeadersReceivedHandler) {
        headersHandler = handler;
        headersCallCount++;
      },
    },
    getPermissionRequestHandler: () => permissionHandler,
    getHeadersReceivedHandler: () => headersHandler,
    getPermissionRequestHandlerCallCount: () => permissionCallCount,
    getHeadersReceivedCallCount: () => headersCallCount,
  };
};

type SessionHardeningModule = {
  hardenSession: (session: FakeSession, isPackaged: boolean, backendPort?: number) => void;
  buildCsp: (isPackaged: boolean, backendPort?: number) => string;
  isValidPort: (port: number | undefined) => boolean;
};

/**
 * Each test needs a fresh module instance because the WeakSet guard
 * (CL-I1 fix) persists across calls within the same module load.
 * We delete the require cache before each load to get a fresh WeakSet.
 */
const loadSessionHardening = (): SessionHardeningModule => {
  const modulePath = require.resolve('./session-hardening');
  delete require.cache[modulePath];
  return require('./session-hardening') as SessionHardeningModule;
};

// ============================================================
// T-I7(a): Permission request handler denies unlisted permissions
// ============================================================

const testUnlistedPermission_isDenied = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false);

  const handler = session.getPermissionRequestHandler();
  assert.ok(handler, 'setPermissionRequestHandler must be called');

  // Request an unlisted permission — should be denied
  let granted: boolean | undefined;
  handler(null, 'camera', (result) => {
    granted = result;
  });
  assert.equal(granted, false, 'Unlisted permission "camera" must be denied');

  // Another unlisted permission
  let granted2: boolean | undefined;
  handler(null, 'geolocation', (result) => {
    granted2 = result;
  });
  assert.equal(
    granted2,
    false,
    'Unlisted permission "geolocation" must be denied'
  );
};

// ============================================================
// T-I7(b): onHeadersReceived injects Content-Security-Policy
// ============================================================

const testHeadersReceived_injectsCSP = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false);

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler, 'onHeadersReceived handler must be registered');

  // Invoke the handler with synthetic details
  let resultHeaders: Record<string, string[]> | undefined;
  handler({ responseHeaders: { 'X-Existing': ['value'] } }, (result) => {
    resultHeaders = result.responseHeaders;
  });

  assert.ok(resultHeaders, 'Callback must receive responseHeaders');
  assert.ok(
    resultHeaders!['Content-Security-Policy'],
    'CSP header must be present'
  );
  assert.equal(
    resultHeaders!['Content-Security-Policy'].length,
    1,
    'CSP header must have exactly one value'
  );
  // Existing headers must be preserved
  assert.deepEqual(
    resultHeaders!['X-Existing'],
    ['value'],
    'Existing response headers must be preserved'
  );
};

// ============================================================
// T-I7(c): Dev CSP includes 'unsafe-eval', production does NOT
// ============================================================

const testDevCSP_includesUnsafeEval = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false); // dev mode

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler);

  let cspValue = '';
  handler({ responseHeaders: {} }, (result) => {
    cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });

  // Dev CSP must contain the quoted token 'unsafe-eval' (not just the
  // substring inside 'wasm-unsafe-eval')
  assert.match(
    cspValue,
    /'unsafe-eval'/,
    "Dev CSP must include 'unsafe-eval'"
  );

  // Dev script-src must also allow 'unsafe-inline' for the @vitejs/plugin-react
  // inline preamble (React Fast Refresh). Scope the match to the script-src
  // directive so we don't accidentally match the style-src 'unsafe-inline'.
  assert.match(
    cspValue,
    /script-src[^;]*'unsafe-inline'/,
    "Dev CSP script-src must include 'unsafe-inline' (Vite React preamble)"
  );
};

const testProductionCSP_doesNotIncludeUnsafeEval = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, true); // packaged mode

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler);

  let cspValue = '';
  handler({ responseHeaders: {} }, (result) => {
    cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });

  // Production CSP must NOT contain the quoted token 'unsafe-eval'.
  // Note: 'wasm-unsafe-eval' contains the bare substring "unsafe-eval",
  // so we match on the apostrophe-delimited form to discriminate.
  assert.doesNotMatch(
    cspValue,
    /'unsafe-eval'/,
    "Production CSP must NOT include 'unsafe-eval'"
  );

  // Verify that wasm-unsafe-eval IS present in production
  assert.match(
    cspValue,
    /'wasm-unsafe-eval'/,
    "Production CSP must include 'wasm-unsafe-eval'"
  );

  // Production script-src must NOT allow inline scripts. (style-src may include
  // 'unsafe-inline', so we scope the negative match to the script-src directive.)
  assert.doesNotMatch(
    cspValue,
    /script-src[^;]*'unsafe-inline'/,
    "Production CSP script-src must NOT include 'unsafe-inline'"
  );
};

// ============================================================
// SEC-L4: CSP baseline directives (object-src, base-uri, form-action)
// ============================================================

const testCSP_containsBaselineDirectives = (): void => {
  const { hardenSession } = loadSessionHardening();

  const directives = ["object-src 'none'", "base-uri 'none'", "form-action 'none'"];

  // Verify dev CSP
  const devSession = createFakeSession();
  hardenSession(devSession, false);
  const devHandler = devSession.getHeadersReceivedHandler();
  assert.ok(devHandler);
  let devCsp = '';
  devHandler({ responseHeaders: {} }, (result) => {
    devCsp = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });
  for (const directive of directives) {
    assert.ok(
      devCsp.includes(directive),
      `Dev CSP must include '${directive}'`
    );
  }

  // Verify production CSP (fresh module to reset WeakSet)
  const { hardenSession: hardenProd } = loadSessionHardening();
  const prodSession = createFakeSession();
  hardenProd(prodSession, true);
  const prodHandler = prodSession.getHeadersReceivedHandler();
  assert.ok(prodHandler);
  let prodCsp = '';
  prodHandler({ responseHeaders: {} }, (result) => {
    prodCsp = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });
  for (const directive of directives) {
    assert.ok(
      prodCsp.includes(directive),
      `Production CSP must include '${directive}'`
    );
  }
};

// ============================================================
// IMP-6: CSP header deduplication (case-insensitive)
// ============================================================

const testHeadersReceived_replacesLowercaseCspHeader = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false);

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler, 'onHeadersReceived handler must be registered');

  // Simulate an upstream response that already has a lowercase CSP header
  let resultHeaders: Record<string, string[]> | undefined;
  handler(
    {
      responseHeaders: {
        'content-security-policy': ["default-src 'none'"],
        'X-Other': ['keep-me'],
      },
    },
    (result) => {
      resultHeaders = result.responseHeaders;
    }
  );

  assert.ok(resultHeaders, 'Callback must receive responseHeaders');

  // The lowercase variant must be removed — only the canonical Title-Case key
  // should exist, and it must carry the injected CSP value (not the upstream one).
  const cspKeys = Object.keys(resultHeaders!).filter(
    (k) => k.toLowerCase() === 'content-security-policy'
  );
  assert.equal(
    cspKeys.length,
    1,
    'Exactly one CSP header key must exist (no duplicates)'
  );
  assert.equal(
    cspKeys[0],
    'Content-Security-Policy',
    'CSP key must be the canonical Title-Case form'
  );
  assert.ok(
    !resultHeaders!['Content-Security-Policy'][0].includes("default-src 'none'"),
    'Injected CSP must replace the upstream value, not preserve it'
  );

  // Other headers must be preserved
  assert.deepEqual(
    resultHeaders!['X-Other'],
    ['keep-me'],
    'Non-CSP headers must be preserved'
  );
};

const testHeadersReceived_undefinedResponseHeaders_isSafe = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false);

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler);

  // Simulate responseHeaders being undefined (some Electron versions)
  let resultHeaders: Record<string, string[]> | undefined;
  handler({ responseHeaders: undefined }, (result) => {
    resultHeaders = result.responseHeaders;
  });

  assert.ok(resultHeaders, 'Callback must receive responseHeaders even when input is undefined');
  assert.ok(
    resultHeaders!['Content-Security-Policy'],
    'CSP header must be present even when upstream headers are undefined'
  );
  assert.equal(
    resultHeaders!['Content-Security-Policy'].length,
    1,
    'CSP must have exactly one value'
  );
};

// ============================================================
// CL-I1: hardenSession is idempotent per session (WeakSet guard)
// ============================================================

const testHardenSession_calledTwice_registersHandlersOnlyOnce = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false);
  hardenSession(session, false); // second call — should be a no-op

  assert.equal(
    session.getPermissionRequestHandlerCallCount(),
    1,
    'setPermissionRequestHandler must be called exactly once per session'
  );
  assert.equal(
    session.getHeadersReceivedCallCount(),
    1,
    'onHeadersReceived must be called exactly once per session'
  );
};

// A2 fix: splash hardens session without port; main window upgrades it
const testHardenSession_portUpgrade_reRegistersHeaders = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  // First call: splash without port
  hardenSession(session, true);
  assert.equal(session.getPermissionRequestHandlerCallCount(), 1);
  assert.equal(session.getHeadersReceivedCallCount(), 1);

  // Verify initial CSP has no port
  let cspValue = '';
  const handler1 = session.getHeadersReceivedHandler();
  assert.ok(handler1);
  handler1({ responseHeaders: {} }, (result) => {
    cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });
  assert.ok(
    !cspValue.includes('127.0.0.1:'),
    'Initial CSP without port must not contain loopback port'
  );

  // Second call: main window with port — should upgrade
  hardenSession(session, true, 5060);
  assert.equal(
    session.getPermissionRequestHandlerCallCount(),
    1,
    'Permission handler must NOT be re-registered on port upgrade'
  );
  assert.equal(
    session.getHeadersReceivedCallCount(),
    2,
    'onHeadersReceived must be re-registered on port upgrade'
  );

  // Verify upgraded CSP has the port
  let upgradedCsp = '';
  const handler2 = session.getHeadersReceivedHandler();
  assert.ok(handler2);
  handler2({ responseHeaders: {} }, (result) => {
    upgradedCsp = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });
  assert.ok(
    upgradedCsp.includes('http://127.0.0.1:5060'),
    'Upgraded CSP must pin to the backend port'
  );

  // Third call with same port: should be a no-op
  hardenSession(session, true, 5060);
  assert.equal(
    session.getHeadersReceivedCallCount(),
    2,
    'onHeadersReceived must NOT be re-registered when already port-pinned'
  );
};

const testHardenSession_differentSessions_bothHardened = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session1 = createFakeSession();
  const session2 = createFakeSession();

  hardenSession(session1, false);
  hardenSession(session2, true);

  // Both sessions must have handlers registered
  assert.ok(
    session1.getPermissionRequestHandler(),
    'First session must have permission handler'
  );
  assert.ok(
    session2.getPermissionRequestHandler(),
    'Second session must have permission handler'
  );
  assert.equal(session1.getPermissionRequestHandlerCallCount(), 1);
  assert.equal(session2.getPermissionRequestHandlerCallCount(), 1);
};

// ============================================================
// M-003: Port-pinned CSP tests
// ============================================================

const testProductionCSP_withBackendPort_pinsConnectSrc = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, true, 5060); // production with known port

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler);

  let cspValue = '';
  handler({ responseHeaders: {} }, (result) => {
    cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });

  // Must contain the specific port, NOT the wildcard
  assert.ok(
    cspValue.includes('http://127.0.0.1:5060'),
    'Production CSP must pin HTTP connect-src to the backend port'
  );
  assert.ok(
    cspValue.includes('ws://127.0.0.1:5060'),
    'Production CSP must pin WS connect-src to the backend port'
  );
  assert.ok(
    !cspValue.includes('127.0.0.1:*'),
    'Production CSP must NOT contain wildcard port when backendPort is provided'
  );
};

const testProductionCSP_withoutBackendPort_noLoopbackWildcard = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, true); // production without backend port (e.g. splash)

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler);

  let cspValue = '';
  handler({ responseHeaders: {} }, (result) => {
    cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });

  // Must NOT contain any loopback wildcard or specific port
  assert.ok(
    !cspValue.includes('127.0.0.1:*'),
    'Production CSP without backendPort must NOT contain wildcard loopback'
  );
  assert.ok(
    !cspValue.includes('127.0.0.1:'),
    'Production CSP without backendPort must NOT contain any loopback port'
  );
  // Must still have connect-src 'self'
  assert.ok(
    cspValue.includes("connect-src 'self'"),
    'Production CSP without backendPort must include connect-src self'
  );
};

// E1: Invalid backendPort values are treated as absent (fail-closed)
const testBuildCsp_invalidPort_failsClosed = (): void => {
  // J12(c): removed unused `hardenSession` destructure

  const invalidPorts = [0, -1, 65536, NaN, 1.5, Infinity, -Infinity];
  for (const port of invalidPorts) {
    const session = createFakeSession();
    // Each needs a fresh module to reset the WeakMap
    const { hardenSession: harden } = loadSessionHardening();
    harden(session, true, port);

    const handler = session.getHeadersReceivedHandler();
    assert.ok(handler, `Handler must be registered for port=${port}`);

    let cspValue = '';
    handler({ responseHeaders: {} }, (result) => {
      cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
    });

    // Invalid port should produce the "no loopback" CSP (fail-closed)
    assert.ok(
      !cspValue.includes('127.0.0.1:'),
      `Invalid port ${port} must not appear in CSP. Got: ${cspValue}`
    );
    assert.ok(
      cspValue.includes("connect-src 'self'"),
      `Invalid port ${port} must fall back to connect-src 'self'`
    );
  }
};

// H1 regression: hardenSession(session, true, NaN) then hardenSession(session, true, 5060)
// must end up with a port-pinned CSP (NaN must NOT lock portPinned:true)
const testHardenSession_nanThenValidPort_upgradesCorrectly = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  // First call with NaN — invalid port, must NOT record portPinned:true
  hardenSession(session, true, NaN);
  assert.equal(session.getHeadersReceivedCallCount(), 1);

  // Verify CSP has no port (NaN treated as absent)
  let cspAfterNaN = '';
  const handler1 = session.getHeadersReceivedHandler();
  assert.ok(handler1);
  handler1({ responseHeaders: {} }, (result) => {
    cspAfterNaN = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });
  assert.ok(
    !cspAfterNaN.includes('127.0.0.1:'),
    'CSP after NaN port must not contain loopback port'
  );

  // Second call with valid port — must upgrade
  hardenSession(session, true, 5060);
  assert.equal(
    session.getHeadersReceivedCallCount(),
    2,
    'onHeadersReceived must be re-registered when upgrading from NaN to valid port'
  );

  // Verify CSP now has the pinned port
  let cspAfterUpgrade = '';
  const handler2 = session.getHeadersReceivedHandler();
  assert.ok(handler2);
  handler2({ responseHeaders: {} }, (result) => {
    cspAfterUpgrade = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });
  assert.ok(
    cspAfterUpgrade.includes('http://127.0.0.1:5060'),
    'CSP must be port-pinned after upgrading from NaN to valid port'
  );
  assert.ok(
    cspAfterUpgrade.includes('ws://127.0.0.1:5060'),
    'CSP must include WS port-pinned entry after upgrade'
  );
};

const testDevCSP_keepsWildcardRegardlessOfPort = (): void => {
  const { hardenSession } = loadSessionHardening();
  const session = createFakeSession();

  hardenSession(session, false, 5060); // dev mode with port

  const handler = session.getHeadersReceivedHandler();
  assert.ok(handler);

  let cspValue = '';
  handler({ responseHeaders: {} }, (result) => {
    cspValue = result.responseHeaders?.['Content-Security-Policy']?.[0] ?? '';
  });

  // Dev mode must keep the wildcard for Vite HMR compatibility
  assert.ok(
    cspValue.includes('127.0.0.1:*'),
    'Dev CSP must keep wildcard port even when backendPort is provided'
  );
};

// ============================================================
// H11(b): Direct buildCsp unit tests
// ============================================================

const testBuildCsp_productionWithPort = (): void => {
  const { buildCsp } = loadSessionHardening();
  const csp = buildCsp(true, 5060);
  assert.ok(csp.includes('http://127.0.0.1:5060'), 'Must pin HTTP to port 5060');
  assert.ok(csp.includes('ws://127.0.0.1:5060'), 'Must pin WS to port 5060');
  assert.ok(!csp.includes('127.0.0.1:*'), 'Must NOT contain wildcard');
  assert.ok(!csp.includes("'unsafe-eval'"), 'Production must NOT include unsafe-eval');
};

const testBuildCsp_productionWithoutPort = (): void => {
  const { buildCsp } = loadSessionHardening();
  const csp = buildCsp(true, undefined);
  assert.ok(csp.includes("connect-src 'self'"), 'Must include connect-src self');
  assert.ok(!csp.includes('127.0.0.1'), 'Must NOT contain any loopback');
};

const testBuildCsp_devMode = (): void => {
  const { buildCsp } = loadSessionHardening();
  const csp = buildCsp(false, 5060);
  assert.ok(csp.includes('127.0.0.1:*'), 'Dev must keep wildcard');
  assert.ok(csp.includes("'unsafe-eval'"), 'Dev must include unsafe-eval');
  assert.ok(csp.includes("'unsafe-inline'"), 'Dev script-src must include unsafe-inline');
};

// ============================================================
// J12(a): Direct isValidPort boundary tests
// ============================================================

const testIsValidPort_boundaries = (): void => {
  const { isValidPort } = loadSessionHardening();

  // Invalid values
  assert.equal(isValidPort(undefined), false, 'undefined must be invalid');
  assert.equal(isValidPort(NaN), false, 'NaN must be invalid');
  assert.equal(isValidPort(0), false, '0 must be invalid');
  assert.equal(isValidPort(65536), false, '65536 must be invalid');
  assert.equal(isValidPort(1.5), false, '1.5 must be invalid');
  assert.equal(isValidPort(-1), false, '-1 must be invalid');
  assert.equal(isValidPort(Infinity), false, 'Infinity must be invalid');

  // Valid boundary values
  assert.equal(isValidPort(1), true, '1 must be valid');
  assert.equal(isValidPort(65535), true, '65535 must be valid');
  assert.equal(isValidPort(8080), true, '8080 must be valid');
};

// ============================================================
// Runner
// ============================================================

const run = (): void => {
  // T-I7(a): Permission denial
  testUnlistedPermission_isDenied();

  // T-I7(b): CSP header injection
  testHeadersReceived_injectsCSP();

  // T-I7(c): Dev vs production CSP
  testDevCSP_includesUnsafeEval();
  testProductionCSP_doesNotIncludeUnsafeEval();

  // SEC-L4: CSP baseline directives
  testCSP_containsBaselineDirectives();

  // IMP-6: CSP header deduplication
  testHeadersReceived_replacesLowercaseCspHeader();
  testHeadersReceived_undefinedResponseHeaders_isSafe();

  // CL-I1: Idempotency guard
  testHardenSession_calledTwice_registersHandlersOnlyOnce();
  testHardenSession_differentSessions_bothHardened();

  // M-003: Port-pinned CSP
  testProductionCSP_withBackendPort_pinsConnectSrc();
  testProductionCSP_withoutBackendPort_noLoopbackWildcard();
  testDevCSP_keepsWildcardRegardlessOfPort();

  // E1: Invalid port validation
  testBuildCsp_invalidPort_failsClosed();

  // A2: Port-upgrade re-hardening
  testHardenSession_portUpgrade_reRegistersHeaders();

  // H1: NaN-then-valid-port regression
  testHardenSession_nanThenValidPort_upgradesCorrectly();

  // H11(b): Direct buildCsp unit tests
  testBuildCsp_productionWithPort();
  testBuildCsp_productionWithoutPort();
  testBuildCsp_devMode();

  // J12(a): Direct isValidPort boundary tests
  testIsValidPort_boundaries();

  console.log('session-hardening tests passed (19 tests)');
};

run();

export {};
