const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type ExternalUrlValidationResult =
  | { valid: true; url: string }
  | { valid: false; error: string };

type IpcSecurityModule = {
  assertAllowedIpcSender: (
    event: { sender: unknown },
    getMainWindow: () => { webContents: unknown; isDestroyed: () => boolean } | null
  ) => void;
  validateExternalUrl: (url: unknown) => ExternalUrlValidationResult;
};

const loadIpcSecurity = (): IpcSecurityModule => {
  return require('./ipc-security') as IpcSecurityModule;
};

// ============================================================
// T-I6: assertAllowedIpcSender tests
// ============================================================

const testSenderIsMainWindow_doesNotThrow = (): void => {
  const { assertAllowedIpcSender } = loadIpcSecurity();

  const sharedWebContents = { id: 1 };
  const event = { sender: sharedWebContents };
  const getMainWindow = () => ({
    webContents: sharedWebContents,
    isDestroyed: () => false,
  });

  // Same object identity — must NOT throw
  assert.doesNotThrow(() => {
    assertAllowedIpcSender(event, getMainWindow);
  });
};

const testSenderIsDifferentObject_throws = (): void => {
  const { assertAllowedIpcSender } = loadIpcSecurity();

  const mainWebContents = { id: 1 };
  const rogueWebContents = { id: 2 };
  const event = { sender: rogueWebContents };
  const getMainWindow = () => ({
    webContents: mainWebContents,
    isDestroyed: () => false,
  });

  // Different object — must throw
  assert.throws(
    () => {
      assertAllowedIpcSender(event, getMainWindow);
    },
    { message: 'IPC rejected: unauthorized sender' }
  );
};

const testMainWindowIsNull_throws = (): void => {
  const { assertAllowedIpcSender } = loadIpcSecurity();

  const event = { sender: { id: 1 } };
  const getMainWindow = () => null;

  // No main window — must throw
  assert.throws(
    () => {
      assertAllowedIpcSender(event, getMainWindow);
    },
    { message: 'IPC rejected: main window is not available' }
  );
};

const testMainWindowIsDestroyed_throws = (): void => {
  const { assertAllowedIpcSender } = loadIpcSecurity();

  const sharedWebContents = { id: 1 };
  const event = { sender: sharedWebContents };
  const getMainWindow = () => ({
    webContents: sharedWebContents,
    isDestroyed: () => true,
  });

  // Destroyed window — must throw
  assert.throws(
    () => {
      assertAllowedIpcSender(event, getMainWindow);
    },
    { message: 'IPC rejected: main window is not available' }
  );
};

// ============================================================
// T-IMP5: validateExternalUrl tests
// ============================================================

const testValidateUrl_nonString_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl(42);
  assert.equal(result.valid, false);
  if (!result.valid) {
    assert.match(result.error, /non-empty string/);
  }

  const result2 = validateExternalUrl(null);
  assert.equal(result2.valid, false);

  const result3 = validateExternalUrl(undefined);
  assert.equal(result3.valid, false);
};

const testValidateUrl_emptyString_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('');
  assert.equal(result.valid, false);
  if (!result.valid) {
    assert.match(result.error, /non-empty string/);
  }
};

const testValidateUrl_overlong_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // 2049 characters — exceeds the 2048 limit
  const longUrl = 'https://example.com/' + 'a'.repeat(2029);
  assert.ok(longUrl.length > 2048);

  const result = validateExternalUrl(longUrl);
  assert.equal(result.valid, false);
  if (!result.valid) {
    assert.match(result.error, /maximum length/);
  }
};

const testValidateUrl_exactlyAtLimit_returnsValid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // Exactly 2048 characters — should pass
  const prefix = 'https://example.com/';
  const url = prefix + 'a'.repeat(2048 - prefix.length);
  assert.equal(url.length, 2048);

  const result = validateExternalUrl(url);
  assert.equal(result.valid, true);
};

const testValidateUrl_malformedUrl_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('not-a-url');
  assert.equal(result.valid, false);
  if (!result.valid) {
    assert.match(result.error, /malformed/);
  }
};

const testValidateUrl_nonHttpProtocol_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('file:///etc/passwd');
  assert.equal(result.valid, false);
  if (!result.valid) {
    assert.match(result.error, /HTTP and HTTPS/);
  }

  const result2 = validateExternalUrl('javascript:alert(1)');
  assert.equal(result2.valid, false);
};

const testValidateUrl_validHttpUrl_returnsValid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('http://example.com');
  assert.equal(result.valid, true);
  if (result.valid) {
    assert.equal(result.url, 'http://example.com');
  }
};

const testValidateUrl_validHttpsUrl_returnsValid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('https://example.com/path?q=1#frag');
  assert.equal(result.valid, true);
  if (result.valid) {
    assert.equal(result.url, 'https://example.com/path?q=1#frag');
  }
};

// ============================================================
// Runner
// ============================================================

const run = (): void => {
  // assertAllowedIpcSender tests
  testSenderIsMainWindow_doesNotThrow();
  testSenderIsDifferentObject_throws();
  testMainWindowIsNull_throws();
  testMainWindowIsDestroyed_throws();

  // validateExternalUrl tests (IMP-5)
  testValidateUrl_nonString_returnsInvalid();
  testValidateUrl_emptyString_returnsInvalid();
  testValidateUrl_overlong_returnsInvalid();
  testValidateUrl_exactlyAtLimit_returnsValid();
  testValidateUrl_malformedUrl_returnsInvalid();
  testValidateUrl_nonHttpProtocol_returnsInvalid();
  testValidateUrl_validHttpUrl_returnsValid();
  testValidateUrl_validHttpsUrl_returnsValid();

  console.log('ipc-security tests passed (12 tests)');
};

run();

export {};
