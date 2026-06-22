import type { IpcMainInvokeEvent, BrowserWindow } from 'electron';

/**
 * Assert that the IPC event sender is the main window's webContents.
 *
 * This prevents rogue renderers (e.g. a compromised webview or devtools page)
 * from invoking privileged IPC handlers. The check uses object identity
 * (event.sender === mainWindow.webContents), which is unforgeable from the
 * renderer side.
 *
 * @throws Error when the sender is not the expected main window.
 */
export const assertAllowedIpcSender = (
  event: IpcMainInvokeEvent,
  getMainWindow: () => BrowserWindow | null
): void => {
  const mainWindow = getMainWindow();

  if (!mainWindow || mainWindow.isDestroyed()) {
    throw new Error('IPC rejected: main window is not available');
  }

  if (event.sender !== mainWindow.webContents) {
    console.warn(
      '[Security] IPC rejected: sender is not the main window webContents'
    );
    throw new Error('IPC rejected: unauthorized sender');
  }
};

/** Maximum length for URLs passed to open-external-url. */
const MAX_EXTERNAL_URL_LENGTH = 2048;

/** Protocols allowed for external URL opening. */
const ALLOWED_PROTOCOLS: ReadonlySet<string> = new Set(['http:', 'https:']);

/**
 * Hostnames that must never be opened via shell.openExternal.
 * These are loopback / link-local addresses that could reach local services
 * or leak data to an attacker-controlled local listener.
 */
const BLOCKED_HOSTNAMES: ReadonlySet<string> = new Set([
  'localhost',
  '127.0.0.1',
  '[::1]',
  '::1',
  '0.0.0.0',
  '[::0]',
  '::0',
]);

/**
 * Check whether a hostname refers to a private/reserved IP range.
 * Covers RFC 1918 (10.x, 172.16-31.x, 192.168.x), link-local (169.254.x),
 * and loopback (127.x). IPv6 private ranges are handled by BLOCKED_HOSTNAMES.
 */
const isPrivateIpHostname = (hostname: string): boolean => {
  // IPv6 or bracket-notation check.
  // H9 fix: any hostname containing ':' is treated as IPv6 (private/blocked)
  // even without brackets. new URL() always bracket-wraps IPv6, so bare ':'
  // is unexploitable via the current call path, but a bare ':' hostname is
  // never a valid public IPv4 address — block defensively.
  if (hostname.includes(':') || hostname.startsWith('[')) {
    return true;
  }

  const parts = hostname.split('.');
  if (parts.length !== 4) return false;

  const octets = parts.map(Number);
  if (octets.some((o) => !Number.isInteger(o) || o < 0 || o > 255)) {
    return false;
  }

  const [a, b] = octets;

  // 127.0.0.0/8 (loopback)
  if (a === 127) return true;
  // 10.0.0.0/8
  if (a === 10) return true;
  // 172.16.0.0/12
  if (a === 172 && b >= 16 && b <= 31) return true;
  // 192.168.0.0/16
  if (a === 192 && b === 168) return true;
  // 169.254.0.0/16 (link-local)
  if (a === 169 && b === 254) return true;
  // 0.0.0.0/8 — the entire 0/8 block is reserved (IANA) and can route to
  // localhost on Linux. Block the whole range, not just exact 0.0.0.0 (C3).
  if (a === 0) return true;
  // J16: Multicast (224.0.0.0/4), reserved (240.0.0.0/4), and broadcast
  // (255.255.255.255). a >= 224 covers all three ranges. These have no
  // legitimate use as shell.openExternal targets.
  if (a >= 224) return true;

  return false;
};

/**
 * Result of validating an external URL.
 * `valid: true` means the URL is safe to open; `valid: false` carries an
 * error message suitable for returning to the renderer.
 */
export type ExternalUrlValidationResult =
  | { valid: true; url: string }
  | { valid: false; error: string };

/**
 * Validate that an unknown value is a well-formed, safe-to-open external URL.
 *
 * Checks performed (in order):
 * 1. Must be a non-empty string.
 * 2. Must not exceed {@link MAX_EXTERNAL_URL_LENGTH} characters.
 * 3. Must parse as a valid URL.
 * 4. Must use an allowed protocol (http: or https:).
 *
 * Known limitation (J3): validation is purely syntactic — a DNS name resolving
 * to 127.0.0.1 (e.g. localtest.me) would pass. This is intentional because:
 * - The URL opens in the user's default system browser via shell.openExternal,
 *   NOT in an in-app fetch. Any web page the user visits can already navigate
 *   to such hosts.
 * - An interstitial confirmation dialog shows the full URL before opening.
 * - Async DNS pre-resolution would add latency/complexity and is subject to
 *   TOCTOU (DNS can change between check and browser open).
 *
 * This is a pure function with no Electron runtime dependency, so it can be
 * unit-tested under the ts-node harness.
 */
export const validateExternalUrl = (
  url: unknown
): ExternalUrlValidationResult => {
  if (typeof url !== 'string' || url.length === 0) {
    return { valid: false, error: 'URL must be a non-empty string' };
  }

  if (url.length > MAX_EXTERNAL_URL_LENGTH) {
    return {
      valid: false,
      error: `URL exceeds maximum length of ${MAX_EXTERNAL_URL_LENGTH} characters`,
    };
  }

  let parsedUrl: URL;
  try {
    parsedUrl = new URL(url);
  } catch {
    return { valid: false, error: 'URL is malformed' };
  }

  if (!ALLOWED_PROTOCOLS.has(parsedUrl.protocol)) {
    return { valid: false, error: 'Only HTTP and HTTPS URLs are allowed' };
  }

  // Reject URLs with embedded credentials (user:pass@host) — these can be
  // used for phishing (displaying a misleading domain in the URL bar).
  if (parsedUrl.username || parsedUrl.password) {
    return { valid: false, error: 'URLs with embedded credentials are not allowed' };
  }

  // Normalize hostname: strip trailing dot (DNS root label) which can
  // bypass hostname checks — 'localhost.' resolves to loopback but would
  // not match BLOCKED_HOSTNAMES without normalization (C1 fix).
  const rawHostname = parsedUrl.hostname.toLowerCase();
  const hostname = rawHostname.endsWith('.') ? rawHostname.slice(0, -1) : rawHostname;
  if (BLOCKED_HOSTNAMES.has(hostname) || isPrivateIpHostname(hostname)) {
    return {
      valid: false,
      error: 'URLs targeting localhost or private networks are not allowed',
    };
  }

  return { valid: true, url };
};
