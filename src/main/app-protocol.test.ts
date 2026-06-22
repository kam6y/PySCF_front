const assert: typeof import('node:assert/strict') = require('node:assert/strict');

/**
 * Tests for the app:// custom protocol handler (M-002).
 */

type AppProtocolHandlerDeps = {
  baseDir: string;
  isPackaged: boolean;
  getBackendPort: () => number | undefined;
  fetchImpl: (url: string) => Promise<Response>;
};

type AppProtocolModule = {
  APP_SCHEME: string;
  APP_HOST: string;
  buildAppUrl: (assetPath: string, query?: Record<string, string>) => string;
  getMimeType: (filePath: string) => string;
  validateRequestPath: (requestPath: string, baseDir: string) => string | null;
  sanitizeForLog: (input: string) => string;
  createAppProtocolRequestHandler: (
    deps: AppProtocolHandlerDeps
  ) => (request: { method: string; url: string }) => Response | Promise<Response>;
};

const loadModule = (): AppProtocolModule => {
  const modulePath = require.resolve('./app-protocol');
  delete require.cache[modulePath];
  return require('./app-protocol') as AppProtocolModule;
};

// ============================================================
// Constants
// ============================================================

const testSchemeAndHost = (): void => {
  const { APP_SCHEME, APP_HOST } = loadModule();
  assert.equal(APP_SCHEME, 'app');
  assert.equal(APP_HOST, 'renderer');
};

// ============================================================
// buildAppUrl
// ============================================================

const testBuildAppUrl_basic = (): void => {
  const { buildAppUrl } = loadModule();
  assert.equal(buildAppUrl('index.html'), 'app://renderer/index.html');
};

const testBuildAppUrl_withQuery = (): void => {
  const { buildAppUrl } = loadModule();
  const url = buildAppUrl('index.html', { backend_port: '5060' });
  assert.equal(url, 'app://renderer/index.html?backend_port=5060');
};

const testBuildAppUrl_splashHtml = (): void => {
  const { buildAppUrl } = loadModule();
  assert.equal(buildAppUrl('splash.html'), 'app://renderer/splash.html');
};

// ============================================================
// getMimeType
// ============================================================

const testMimeType_html = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/index.html'), 'text/html');
  assert.equal(getMimeType('/dist/splash.html'), 'text/html');
};

const testMimeType_javascript = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/assets/main.js'), 'application/javascript');
  assert.equal(getMimeType('/dist/assets/chunk.mjs'), 'application/javascript');
};

const testMimeType_css = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/assets/style.css'), 'text/css');
};

const testMimeType_fonts = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/assets/font.woff2'), 'font/woff2');
  assert.equal(getMimeType('/dist/assets/font.woff'), 'font/woff');
  assert.equal(getMimeType('/dist/assets/font.ttf'), 'font/ttf');
};

const testMimeType_images = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/assets/icon.png'), 'image/png');
  assert.equal(getMimeType('/dist/assets/icon.svg'), 'image/svg+xml');
  assert.equal(getMimeType('/dist/assets/photo.jpg'), 'image/jpeg');
};

const testMimeType_wasm = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/assets/module.wasm'), 'application/wasm');
};

const testMimeType_unknown = (): void => {
  const { getMimeType } = loadModule();
  assert.equal(getMimeType('/dist/assets/file.xyz'), 'application/octet-stream');
  assert.equal(getMimeType('/dist/noext'), 'application/octet-stream');
};

// ============================================================
// validateRequestPath — Path traversal rejection
// ============================================================

const BASE_DIR = '/app/dist';

const testValidPath_indexHtml = (): void => {
  const { validateRequestPath } = loadModule();
  const result = validateRequestPath('/index.html', BASE_DIR);
  assert.ok(result, 'index.html must be allowed');
  assert.equal(result, '/app/dist/index.html');
};

const testValidPath_splashHtml = (): void => {
  const { validateRequestPath } = loadModule();
  const result = validateRequestPath('/splash.html', BASE_DIR);
  assert.ok(result, 'splash.html must be allowed');
  assert.equal(result, '/app/dist/splash.html');
};

const testValidPath_assetsSubdir = (): void => {
  const { validateRequestPath } = loadModule();
  const result = validateRequestPath('/assets/main.js', BASE_DIR);
  assert.ok(result, 'assets/main.js must be allowed');
  assert.equal(result, '/app/dist/assets/main.js');
};

const testValidPath_emptyPathServesIndex = (): void => {
  const { validateRequestPath } = loadModule();
  const result = validateRequestPath('/', BASE_DIR);
  assert.ok(result, 'Empty path should serve index.html');
  assert.equal(result, '/app/dist/index.html');
};

const testTraversal_dotDot = (): void => {
  const { validateRequestPath } = loadModule();
  assert.equal(
    validateRequestPath('/../secret/keys.txt', BASE_DIR),
    null,
    'Path with .. must be rejected'
  );
  assert.equal(
    validateRequestPath('/assets/../../etc/passwd', BASE_DIR),
    null,
    'Path with embedded .. must be rejected'
  );
};

const testTraversal_encodedDotDot = (): void => {
  const { validateRequestPath } = loadModule();
  // %2e%2e is URL-encoded ".."
  assert.equal(
    validateRequestPath('/%2e%2e/secret', BASE_DIR),
    null,
    'URL-encoded .. (%2e%2e) must be rejected'
  );
  assert.equal(
    validateRequestPath('/%2E%2E/secret', BASE_DIR),
    null,
    'URL-encoded .. (%2E%2E) must be rejected'
  );
};

const testTraversal_mixedEncodedDot = (): void => {
  const { validateRequestPath } = loadModule();
  // Mixed encoding: one dot encoded, one not
  assert.equal(
    validateRequestPath('/%2e./secret', BASE_DIR),
    null,
    'Mixed-encoded dot (%2e.) must be rejected'
  );
  assert.equal(
    validateRequestPath('/.%2e/secret', BASE_DIR),
    null,
    'Mixed-encoded dot (.%2e) must be rejected'
  );
};

const testTraversal_absolutePathEscape = (): void => {
  const { validateRequestPath } = loadModule();
  // H10: Definitive assertions — these must NOT silently pass when result is null.

  // On POSIX, '/C:/Windows/System32' after leading-slash removal =
  // 'C:/Windows/System32' which is not absolute on POSIX, so it resolves
  // under baseDir as a literal directory name — this is safe.
  const result1 = validateRequestPath('/C:/Windows/System32/cmd.exe', BASE_DIR);
  assert.ok(
    result1 !== null,
    'Windows-style path on POSIX must resolve within baseDir (not be rejected)'
  );
  assert.equal(
    result1,
    '/app/dist/C:/Windows/System32/cmd.exe',
    'Windows-style path must resolve as literal directory name under baseDir'
  );

  // Double-slash prefix: '//evil-server/share/file' after leading-slash removal
  // becomes 'evil-server/share/file' (all leading slashes are stripped by
  // the /^\/+/ regex). On POSIX, 'evil-server/share/file' is relative, so
  // path.resolve yields '/app/dist/evil-server/share/file' — contained.
  // J11: Assert the exact resolved path deterministically (no if/else accept-both).
  const result2 = validateRequestPath('//evil-server/share/file', BASE_DIR);
  assert.ok(
    result2 !== null,
    'Double-slash path on POSIX must resolve within baseDir (not be rejected)'
  );
  assert.equal(
    result2,
    '/app/dist/evil-server/share/file',
    'Double-slash path must resolve as literal child under baseDir'
  );
};

const testTraversal_nullByteRejection = (): void => {
  const { validateRequestPath } = loadModule();
  // Null bytes in paths can cause truncation in native file APIs (D2).
  // decodeURIComponent('%00') produces a null byte character.
  // These must be definitively rejected.
  assert.equal(
    validateRequestPath('/index.html%00.txt', BASE_DIR),
    null,
    'Path with null byte (%00) must be rejected'
  );
  assert.equal(
    validateRequestPath('/%00', BASE_DIR),
    null,
    'Path starting with null byte must be rejected'
  );
  assert.equal(
    validateRequestPath('/assets/file%00.js', BASE_DIR),
    null,
    'Path with embedded null byte must be rejected'
  );
};

const testTraversal_doubleEncodedDotDot = (): void => {
  const { validateRequestPath } = loadModule();
  // Double-encoded: %252e%252e decodes to literal '%2e%2e' (not '..').
  // The raw path '/%252e%252e/secret' does NOT contain the pattern '%2e%2e'
  // directly (it contains '%252e%252e'). After decodeURIComponent it becomes
  // '/%2e%2e/secret'. Since path.resolve treats '%2e%2e' as a literal
  // directory name (not traversal), the resolved path stays within baseDir.
  // This is safe — no actual traversal occurs. Assert it resolves within baseDir.
  const result = validateRequestPath('/%252e%252e/secret', BASE_DIR);
  assert.ok(
    result !== null,
    'Double-encoded %252e%252e is not actual traversal and should be allowed'
  );
  assert.ok(
    result!.startsWith('/app/dist/'),
    'Double-encoded path must resolve within baseDir'
  );
  assert.equal(
    result,
    '/app/dist/%2e%2e/secret',
    'Double-encoded path must resolve to literal %2e%2e directory name'
  );
};

// J12(b): Backslash traversal pattern tests
const testTraversal_backslashVariants = (): void => {
  const { validateRequestPath } = loadModule();

  // Raw backslash: on POSIX this is a literal char, but our TRAVERSAL_PATTERNS
  // fast-reject catches encoded backslash variants for defense-in-depth.

  // ..%5c (encoded backslash after ..)
  assert.equal(
    validateRequestPath('/..%5cetc/passwd', BASE_DIR),
    null,
    '..%5c must be rejected'
  );
  // ..%5C (uppercase)
  assert.equal(
    validateRequestPath('/..%5Cetc/passwd', BASE_DIR),
    null,
    '..%5C must be rejected'
  );
  // %5c.. (backslash before ..)
  assert.equal(
    validateRequestPath('/%5c../etc/passwd', BASE_DIR),
    null,
    '%5c.. must be rejected'
  );
  // %5C.. (uppercase)
  assert.equal(
    validateRequestPath('/%5C../etc/passwd', BASE_DIR),
    null,
    '%5C.. must be rejected'
  );
};

// ============================================================
// sanitizeForLog (D3)
// ============================================================

const testSanitizeForLog_shortString = (): void => {
  const { sanitizeForLog } = loadModule();
  assert.equal(sanitizeForLog('normal/path'), 'normal/path');
};

const testSanitizeForLog_truncatesLongString = (): void => {
  const { sanitizeForLog } = loadModule();
  const long = 'a'.repeat(300);
  const result = sanitizeForLog(long);
  assert.ok(result.length < 300, 'Must be truncated');
  assert.ok(result.endsWith('...(truncated)'), 'Must end with truncation marker');
};

const testSanitizeForLog_stripsControlChars = (): void => {
  const { sanitizeForLog } = loadModule();
  const result = sanitizeForLog('path\x00with\nnewline\rand\ttab');
  assert.ok(!result.includes('\x00'), 'Must strip null byte');
  assert.ok(!result.includes('\n'), 'Must strip newline');
  assert.ok(!result.includes('\r'), 'Must strip carriage return');
  assert.ok(!result.includes('\t'), 'Must strip tab');
};

const testValidPath_deeplyNested = (): void => {
  const { validateRequestPath } = loadModule();
  const result = validateRequestPath('/assets/fonts/sub/deep/file.woff2', BASE_DIR);
  assert.ok(result, 'Deeply nested valid path must be allowed');
  assert.equal(result, '/app/dist/assets/fonts/sub/deep/file.woff2');
};

// ============================================================
// J5: createAppProtocolRequestHandler tests
// ============================================================

const createTestHandler = (
  overrides: Partial<AppProtocolHandlerDeps> = {}
): ((request: { method: string; url: string }) => Response | Promise<Response>) => {
  const { createAppProtocolRequestHandler } = loadModule();
  return createAppProtocolRequestHandler({
    baseDir: BASE_DIR,
    isPackaged: true,
    getBackendPort: () => 5060,
    fetchImpl: () => Promise.resolve(new Response('ok', { status: 200 })),
    ...overrides,
  });
};

const testHandler_nonGetMethod_returns405 = async (): Promise<void> => {
  const handler = createTestHandler();
  const response = await handler({ method: 'POST', url: 'app://renderer/index.html' });
  assert.equal(response.status, 405, 'Non-GET must return 405');
  assert.ok(
    response.headers.get('Content-Security-Policy'),
    '405 response must include CSP header'
  );
};

const testHandler_badUrl_returns400 = async (): Promise<void> => {
  const handler = createTestHandler();
  const response = await handler({ method: 'GET', url: ':::invalid' });
  assert.equal(response.status, 400, 'Bad URL must return 400');
};

const testHandler_traversalPath_returns403 = async (): Promise<void> => {
  const handler = createTestHandler();
  // Note: new URL() normalizes .. and %2e%2e, so standard traversal patterns
  // are resolved by the URL parser before reaching validateRequestPath.
  // Test with null byte injection which validateRequestPath rejects (D2).
  const response = await handler({ method: 'GET', url: 'app://renderer/index.html%00.txt' });
  assert.equal(response.status, 403, 'Null byte path must return 403');
};

const testHandler_validPath_fetchOk_returns200 = async (): Promise<void> => {
  const handler = createTestHandler({
    fetchImpl: () => Promise.resolve(new Response('file content', { status: 200 })),
  });
  const response = await handler({ method: 'GET', url: 'app://renderer/index.html' });
  assert.equal(response.status, 200, 'Valid path with ok fetch must return 200');
  assert.equal(response.headers.get('Content-Type'), 'text/html');
  assert.ok(response.headers.get('Content-Security-Policy'), 'Must include CSP');
  assert.equal(response.headers.get('X-Content-Type-Options'), 'nosniff');
};

const testHandler_validPath_fetchNotOk_returns404 = async (): Promise<void> => {
  const handler = createTestHandler({
    fetchImpl: () => Promise.resolve(new Response('not found', { status: 404 })),
  });
  const response = await handler({ method: 'GET', url: 'app://renderer/missing.js' });
  assert.equal(response.status, 404, 'Non-OK fetch must return 404');
};

const testHandler_validPath_fetchRejects_returns404 = async (): Promise<void> => {
  const handler = createTestHandler({
    fetchImpl: () => Promise.reject(new Error('EACCES: permission denied')),
  });
  const response = await handler({ method: 'GET', url: 'app://renderer/index.html' });
  assert.equal(response.status, 404, 'Fetch rejection must return 404');
  assert.ok(
    response.headers.get('Content-Security-Policy'),
    '404 from rejection must include CSP'
  );
};

const testHandler_cspIncludesBackendPort = async (): Promise<void> => {
  const handler = createTestHandler({
    getBackendPort: () => 9999,
    fetchImpl: () => Promise.resolve(new Response('ok', { status: 200 })),
  });
  const response = await handler({ method: 'GET', url: 'app://renderer/index.html' });
  const csp = response.headers.get('Content-Security-Policy') ?? '';
  assert.ok(csp.includes('127.0.0.1:9999'), 'CSP must pin to provided backend port');
};

const testHandler_cspWithoutPort_noLoopback = async (): Promise<void> => {
  const handler = createTestHandler({
    getBackendPort: () => undefined,
    fetchImpl: () => Promise.resolve(new Response('ok', { status: 200 })),
  });
  const response = await handler({ method: 'GET', url: 'app://renderer/index.html' });
  const csp = response.headers.get('Content-Security-Policy') ?? '';
  assert.ok(!csp.includes('127.0.0.1:'), 'CSP without port must not contain loopback port');
  assert.ok(csp.includes("connect-src 'self'"), 'Must have connect-src self');
};

const testHandler_mimeTypeForJs = async (): Promise<void> => {
  const handler = createTestHandler({
    fetchImpl: () => Promise.resolve(new Response('code', { status: 200 })),
  });
  const response = await handler({ method: 'GET', url: 'app://renderer/assets/main.js' });
  assert.equal(response.headers.get('Content-Type'), 'application/javascript');
};

// ============================================================
// Runner
// ============================================================

const run = async (): Promise<void> => {
  // Constants
  testSchemeAndHost();

  // buildAppUrl
  testBuildAppUrl_basic();
  testBuildAppUrl_withQuery();
  testBuildAppUrl_splashHtml();

  // getMimeType
  testMimeType_html();
  testMimeType_javascript();
  testMimeType_css();
  testMimeType_fonts();
  testMimeType_images();
  testMimeType_wasm();
  testMimeType_unknown();

  // validateRequestPath — valid paths
  testValidPath_indexHtml();
  testValidPath_splashHtml();
  testValidPath_assetsSubdir();
  testValidPath_emptyPathServesIndex();
  testValidPath_deeplyNested();

  // validateRequestPath — path traversal rejection
  testTraversal_dotDot();
  testTraversal_encodedDotDot();
  testTraversal_mixedEncodedDot();
  testTraversal_absolutePathEscape();
  testTraversal_nullByteRejection();
  testTraversal_doubleEncodedDotDot();

  // J12(b): Backslash traversal patterns
  testTraversal_backslashVariants();

  // sanitizeForLog (D3)
  testSanitizeForLog_shortString();
  testSanitizeForLog_truncatesLongString();
  testSanitizeForLog_stripsControlChars();

  // J5: createAppProtocolRequestHandler tests
  await testHandler_nonGetMethod_returns405();
  await testHandler_badUrl_returns400();
  await testHandler_traversalPath_returns403();
  await testHandler_validPath_fetchOk_returns200();
  await testHandler_validPath_fetchNotOk_returns404();
  await testHandler_validPath_fetchRejects_returns404();
  await testHandler_cspIncludesBackendPort();
  await testHandler_cspWithoutPort_noLoopback();
  await testHandler_mimeTypeForJs();

  console.log('app-protocol tests passed (39 tests)');
};

run().catch((err) => {
  console.error(err);
  process.exit(1);
});

export {};
