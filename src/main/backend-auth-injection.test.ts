const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type BeforeSendHandler = (
  details: { url: string; requestHeaders: Record<string, string> },
  callback: (response: { requestHeaders?: Record<string, string> }) => void
) => void;

type FakeSession = {
  webRequest: {
    onBeforeSendHeaders: (
      filter: { urls: string[] },
      handler: BeforeSendHandler
    ) => void;
  };
  getCapturedFilter: () => { urls: string[] } | null;
  getCapturedHandler: () => BeforeSendHandler | null;
  getRegistrationCount: () => number;
};

const createFakeSession = (): FakeSession => {
  let capturedFilter: { urls: string[] } | null = null;
  let capturedHandler: BeforeSendHandler | null = null;
  let registrationCount = 0;

  return {
    webRequest: {
      onBeforeSendHeaders(
        filter: { urls: string[] },
        handler: BeforeSendHandler
      ) {
        capturedFilter = filter;
        capturedHandler = handler;
        registrationCount++;
      },
    },
    getCapturedFilter: () => capturedFilter,
    getCapturedHandler: () => capturedHandler,
    getRegistrationCount: () => registrationCount,
  };
};

type BackendAuthInjectionModule = {
  registerBackendAuthInjection: (
    session: FakeSession,
    backendPort: number,
    authToken: string
  ) => void;
};

const loadModule = (): BackendAuthInjectionModule => {
  const modulePath = require.resolve('./backend-auth-injection');
  delete require.cache[modulePath];
  return require('./backend-auth-injection') as BackendAuthInjectionModule;
};

// ============================================================
// T-SEC002(a): URL filter is port-agnostic (covers all loopback)
// ============================================================

const testUrlFilterIsPortAgnostic = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'test-token');

  const filter = session.getCapturedFilter();
  assert.ok(filter, 'onBeforeSendHeaders must be called');
  assert.deepEqual(
    filter.urls,
    ['http://127.0.0.1/*', 'http://127.0.0.1:*/*'],
    'URL filter must be port-agnostic to avoid Chromium match-pattern port issues'
  );
};

// ============================================================
// T-SEC002(b): Handler preserves existing headers and injects token
// ============================================================

const testHandlerInjectsTokenAndPreservesHeaders = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'my-secret-token');

  const handler = session.getCapturedHandler();
  assert.ok(handler, 'Handler must be registered');

  let result: { requestHeaders?: Record<string, string> } | undefined;
  handler(
    { url: 'http://127.0.0.1:5057/api/health', requestHeaders: { Accept: 'text/event-stream' } },
    (response) => {
      result = response;
    }
  );

  assert.ok(result, 'Callback must be invoked');
  assert.ok(result.requestHeaders, 'requestHeaders must be present');
  assert.equal(
    result.requestHeaders!['Accept'],
    'text/event-stream',
    'Existing Accept header must be preserved'
  );
  assert.equal(
    result.requestHeaders!['X-Auth-Token'],
    'my-secret-token',
    'X-Auth-Token must be injected'
  );
};

// ============================================================
// T-SEC002(c): Injected token value matches the authToken argument
// ============================================================

const testInjectedTokenMatchesArgument = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  const token = 'unique-token-value-abc123';
  registerBackendAuthInjection(session, 8080, token);

  const handler = session.getCapturedHandler();
  assert.ok(handler);

  let result: { requestHeaders?: Record<string, string> } | undefined;
  handler({ url: 'http://127.0.0.1:8080/api/calc', requestHeaders: {} }, (response) => {
    result = response;
  });

  assert.equal(
    result!.requestHeaders!['X-Auth-Token'],
    token,
    'Token value must exactly match the authToken argument'
  );
};

// ============================================================
// T-SEC002(d): Idempotency — second call on same session is a no-op
// ============================================================

const testIdempotency_calledTwice_registersOnlyOnce = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'tok');
  registerBackendAuthInjection(session, 5057, 'tok');

  assert.equal(
    session.getRegistrationCount(),
    1,
    'onBeforeSendHeaders must be called exactly once per session'
  );
};

// ============================================================
// T-SEC002(e): Idempotency — different sessions are both registered
// ============================================================

const testIdempotency_differentSessions_bothRegistered = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session1 = createFakeSession();
  const session2 = createFakeSession();

  registerBackendAuthInjection(session1, 5057, 'tok1');
  registerBackendAuthInjection(session2, 5058, 'tok2');

  assert.equal(session1.getRegistrationCount(), 1);
  assert.equal(session2.getRegistrationCount(), 1);
  assert.ok(session1.getCapturedHandler(), 'First session must have handler');
  assert.ok(session2.getCapturedHandler(), 'Second session must have handler');
};

// ============================================================
// T-SEC002(f): Handler with empty requestHeaders still injects token
// ============================================================

const testHandlerWithEmptyHeaders_injectsToken = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'my-token');

  const handler = session.getCapturedHandler();
  assert.ok(handler);

  // Simulate undefined requestHeaders (some Electron versions may pass this)
  let result: { requestHeaders?: Record<string, string> } | undefined;
  handler(
    { url: 'http://127.0.0.1:5057/api/health', requestHeaders: undefined as unknown as Record<string, string> },
    (response) => {
      result = response;
    }
  );

  assert.ok(result, 'Callback must be invoked');
  assert.ok(result.requestHeaders, 'requestHeaders must be present');
  assert.equal(
    result.requestHeaders!['X-Auth-Token'],
    'my-token',
    'X-Auth-Token must be injected even when original headers are undefined'
  );
};

// ============================================================
// T-SEC002(g): Invalid arguments throw
// ============================================================

const testInvalidPort_throws = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  assert.throws(
    () => registerBackendAuthInjection(session, NaN, 'tok'),
    { message: /backendPort must be a positive integer/ },
    'NaN port must throw'
  );

  assert.throws(
    () => registerBackendAuthInjection(session, 0, 'tok'),
    { message: /backendPort must be a positive integer/ },
    'Zero port must throw'
  );

  assert.throws(
    () => registerBackendAuthInjection(session, -1, 'tok'),
    { message: /backendPort must be a positive integer/ },
    'Negative port must throw'
  );

  assert.throws(
    () => registerBackendAuthInjection(session, 3.14, 'tok'),
    { message: /backendPort must be a positive integer/ },
    'Non-integer port must throw'
  );
};

const testEmptyToken_throws = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  assert.throws(
    () => registerBackendAuthInjection(session, 5057, ''),
    { message: /authToken must be a non-empty string/ },
    'Empty token must throw'
  );
};

// ============================================================
// T-SEC002(h): Port gate — matching port injects token
// ============================================================

const testPortGate_matchingPort_injectsToken = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'secret-tok');

  const handler = session.getCapturedHandler();
  assert.ok(handler);

  let result: { requestHeaders?: Record<string, string> } | undefined;
  handler(
    { url: 'http://127.0.0.1:5057/health', requestHeaders: { Accept: 'application/json' } },
    (response) => {
      result = response;
    }
  );

  assert.ok(result);
  assert.equal(
    result.requestHeaders!['X-Auth-Token'],
    'secret-tok',
    'Token must be injected when port matches backendPort'
  );
  assert.equal(
    result.requestHeaders!['Accept'],
    'application/json',
    'Existing headers must be preserved'
  );
};

// ============================================================
// T-SEC002(i): Port gate — different port does NOT inject token
// ============================================================

const testPortGate_differentPort_passesThrough = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'secret-tok');

  const handler = session.getCapturedHandler();
  assert.ok(handler);

  const originalHeaders = { Accept: 'text/html', 'Content-Type': 'application/json' };
  let result: { requestHeaders?: Record<string, string> } | undefined;
  handler(
    { url: 'http://127.0.0.1:9999/x', requestHeaders: originalHeaders },
    (response) => {
      result = response;
    }
  );

  assert.ok(result);
  assert.equal(
    result.requestHeaders!['X-Auth-Token'],
    undefined,
    'Token must NOT be injected when port does not match'
  );
  assert.deepEqual(
    result.requestHeaders,
    originalHeaders,
    'Original headers must be passed through unchanged'
  );
};

// ============================================================
// T-SEC002(j): Malformed URL does not throw and does not inject
// ============================================================

const testPortGate_malformedUrl_passesThrough = (): void => {
  const { registerBackendAuthInjection } = loadModule();
  const session = createFakeSession();

  registerBackendAuthInjection(session, 5057, 'secret-tok');

  const handler = session.getCapturedHandler();
  assert.ok(handler);

  const originalHeaders = { Accept: 'text/html' };
  let result: { requestHeaders?: Record<string, string> } | undefined;
  handler(
    { url: 'not-a-valid-url', requestHeaders: originalHeaders },
    (response) => {
      result = response;
    }
  );

  assert.ok(result, 'Callback must be invoked even for malformed URL');
  assert.equal(
    result.requestHeaders!['X-Auth-Token'],
    undefined,
    'Token must NOT be injected for malformed URL'
  );
  assert.deepEqual(
    result.requestHeaders,
    originalHeaders,
    'Original headers must be passed through unchanged for malformed URL'
  );
};

// ============================================================
// Runner
// ============================================================

const run = (): void => {
  testUrlFilterIsPortAgnostic();
  testHandlerInjectsTokenAndPreservesHeaders();
  testInjectedTokenMatchesArgument();
  testIdempotency_calledTwice_registersOnlyOnce();
  testIdempotency_differentSessions_bothRegistered();
  testHandlerWithEmptyHeaders_injectsToken();
  testInvalidPort_throws();
  testEmptyToken_throws();
  testPortGate_matchingPort_injectsToken();
  testPortGate_differentPort_passesThrough();
  testPortGate_malformedUrl_passesThrough();
  console.log('backend-auth-injection tests passed (11 tests)');
};

run();

export {};
