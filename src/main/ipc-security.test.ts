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
// M-001: Private/loopback URL blocking tests
// ============================================================

const testValidateUrl_localhost_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('http://localhost:8080/api');
  assert.equal(result.valid, false);
  if (!result.valid) {
    assert.match(result.error, /localhost|private/i);
  }

  const result2 = validateExternalUrl('https://localhost/path');
  assert.equal(result2.valid, false);
};

const testValidateUrl_loopback127_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('http://127.0.0.1:5000/api');
  assert.equal(result.valid, false);

  const result2 = validateExternalUrl('https://127.0.0.1/');
  assert.equal(result2.valid, false);
};

const testValidateUrl_ipv6Loopback_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  const result = validateExternalUrl('http://[::1]:8080/api');
  assert.equal(result.valid, false);
};

const testValidateUrl_privateRfc1918_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // 10.0.0.0/8
  assert.equal(validateExternalUrl('http://10.0.0.1/').valid, false);
  // 172.16.0.0/12
  assert.equal(validateExternalUrl('http://172.16.0.1/').valid, false);
  assert.equal(validateExternalUrl('http://172.31.255.255/').valid, false);
  // 192.168.0.0/16
  assert.equal(validateExternalUrl('http://192.168.1.1/').valid, false);
  // 169.254.0.0/16 (link-local)
  assert.equal(validateExternalUrl('http://169.254.1.1/').valid, false);
};

const testValidateUrl_publicIp_returnsValid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // Public IP should be allowed
  assert.equal(validateExternalUrl('https://8.8.8.8/').valid, true);
  // 172.32.x.x is outside RFC 1918 range
  assert.equal(validateExternalUrl('https://172.32.0.1/').valid, true);
};

const testValidateUrl_zeroAddress_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  assert.equal(validateExternalUrl('http://0.0.0.0/').valid, false);
};

// C1: Trailing-dot hostname bypass
const testValidateUrl_trailingDotHostname_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // 'localhost.' with trailing dot resolves to loopback — must be blocked
  assert.equal(
    validateExternalUrl('http://localhost.:8080/').valid,
    false,
    'localhost. (trailing dot) must be blocked'
  );
};

// C3: Entire 0.0.0.0/8 range is blocked
const testValidateUrl_zeroSlash8Range_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  assert.equal(
    validateExternalUrl('http://0.0.0.1/').valid,
    false,
    '0.0.0.1 (in 0/8 range) must be blocked'
  );
  assert.equal(
    validateExternalUrl('http://0.1.2.3/').valid,
    false,
    '0.1.2.3 (in 0/8 range) must be blocked'
  );
  assert.equal(
    validateExternalUrl('http://0.255.255.255/').valid,
    false,
    '0.255.255.255 (in 0/8 range) must be blocked'
  );
};

// F7: Non-loopback IPv6 is also blocked (blanket policy)
const testValidateUrl_nonLoopbackIpv6_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  assert.equal(
    validateExternalUrl('http://[2001:db8::1]/').valid,
    false,
    'Non-loopback IPv6 [2001:db8::1] must be blocked by blanket IPv6 policy'
  );
  assert.equal(
    validateExternalUrl('http://[fe80::1]/').valid,
    false,
    'Link-local IPv6 [fe80::1] must be blocked'
  );
};

// MINOR: URLs with embedded credentials are rejected
const testValidateUrl_embeddedCredentials_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  assert.equal(
    validateExternalUrl('https://user:pass@example.com/').valid,
    false,
    'URL with user:pass@ must be rejected'
  );
  assert.equal(
    validateExternalUrl('https://user@example.com/').valid,
    false,
    'URL with user@ must be rejected'
  );
};

// H9: Bare colon hostname (unbracket IPv6) is blocked
const testValidateUrl_bareColonHostname_returnsInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // Hostname with ':' but no brackets — isPrivateIpHostname must return true
  // new URL normalizes these to bracket-form, but testing the exported function
  // behavior via the public API: http://::1:8080/ is parsed by URL as
  // having hostname '[::1]' which is already blocked. Test a synthetic case
  // where a colon appears — URL('http://foo:bar@evil.com/') puts 'foo' in
  // username, so instead test a known IPv6 that URL brackets:
  assert.equal(
    validateExternalUrl('http://[::ffff:127.0.0.1]/').valid,
    false,
    'IPv6-mapped IPv4 loopback must be blocked'
  );
};

// ============================================================
// J15: IP format bypass tests — WHATWG URL parser canonicalizes
// ends-in-number hosts to dotted-decimal 127.0.0.1, so these
// are all blocked. Pin this behavior against parser changes.
// ============================================================

const testValidateUrl_ipFormatBypasses_allInvalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // Decimal (2130706433 = 127.0.0.1)
  assert.equal(
    validateExternalUrl('http://2130706433/').valid,
    false,
    'Decimal IP 2130706433 must be blocked (resolves to 127.0.0.1)'
  );
  // Hex
  assert.equal(
    validateExternalUrl('http://0x7f000001/').valid,
    false,
    'Hex IP 0x7f000001 must be blocked'
  );
  // Octal (full)
  assert.equal(
    validateExternalUrl('http://017700000001/').valid,
    false,
    'Octal IP 017700000001 must be blocked'
  );
  // Octal first octet
  assert.equal(
    validateExternalUrl('http://0177.0.0.1/').valid,
    false,
    'Octal-dotted 0177.0.0.1 must be blocked'
  );
  // Hex first octet
  assert.equal(
    validateExternalUrl('http://0x7f.0.0.1/').valid,
    false,
    'Hex-dotted 0x7f.0.0.1 must be blocked'
  );
};

// ============================================================
// J16: Multicast/reserved/broadcast IP ranges
// ============================================================

const testValidateUrl_multicastReservedBroadcast_invalid = (): void => {
  const { validateExternalUrl } = loadIpcSecurity();

  // Multicast 224.0.0.0/4
  assert.equal(
    validateExternalUrl('http://224.0.0.1/').valid,
    false,
    'Multicast 224.0.0.1 must be blocked'
  );
  assert.equal(
    validateExternalUrl('http://239.255.255.255/').valid,
    false,
    'Multicast 239.255.255.255 must be blocked'
  );
  // Reserved 240.0.0.0/4
  assert.equal(
    validateExternalUrl('http://240.0.0.1/').valid,
    false,
    'Reserved 240.0.0.1 must be blocked'
  );
  // Broadcast
  assert.equal(
    validateExternalUrl('http://255.255.255.255/').valid,
    false,
    'Broadcast 255.255.255.255 must be blocked'
  );
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

  // M-001: Private/loopback URL blocking
  testValidateUrl_localhost_returnsInvalid();
  testValidateUrl_loopback127_returnsInvalid();
  testValidateUrl_ipv6Loopback_returnsInvalid();
  testValidateUrl_privateRfc1918_returnsInvalid();
  testValidateUrl_publicIp_returnsValid();
  testValidateUrl_zeroAddress_returnsInvalid();

  // C1: Trailing-dot hostname bypass
  testValidateUrl_trailingDotHostname_returnsInvalid();

  // C3: Entire 0/8 range
  testValidateUrl_zeroSlash8Range_returnsInvalid();

  // F7: Non-loopback IPv6 blanket block
  testValidateUrl_nonLoopbackIpv6_returnsInvalid();

  // MINOR: Embedded credentials
  testValidateUrl_embeddedCredentials_returnsInvalid();

  // H9: Bare colon hostname
  testValidateUrl_bareColonHostname_returnsInvalid();

  // J15: IP format bypass tests
  testValidateUrl_ipFormatBypasses_allInvalid();

  // J16: Multicast/reserved/broadcast
  testValidateUrl_multicastReservedBroadcast_invalid();

  console.log('ipc-security tests passed (25 tests)');
};

run();

export {};
