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

// ============================================================
// External URL validation
// ============================================================

/** Maximum length for URLs passed to open-external-url. */
const MAX_EXTERNAL_URL_LENGTH = 2048;

/** Protocols allowed for external URL opening. */
const ALLOWED_PROTOCOLS: ReadonlySet<string> = new Set(['http:', 'https:']);

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

  return { valid: true, url };
};
