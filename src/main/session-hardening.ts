import type { Session } from 'electron';

/**
 * Track sessions that have already been hardened.
 * Both splash and main windows share session.defaultSession (no partition),
 * and Electron's onHeadersReceived / setPermissionRequestHandler are
 * single-handler APIs (last-write-wins). The WeakSet ensures we only
 * register handlers once per session, avoiding silent overwrites.
 */
const hardenedSessions = new WeakSet<Session>();

/**
 * Permissions that the renderer is allowed to request.
 * Everything not in this set is denied by default.
 * Currently empty — the app does not need any special permissions
 * (no camera, microphone, geolocation, notifications, etc.).
 */
const ALLOWED_PERMISSIONS: ReadonlySet<string> = new Set([
  // Add permissions here only when genuinely needed, e.g. 'clipboard-read'
]);

/**
 * Production CSP delivered as a response header.
 *
 * This mirrors the meta CSP in index.html but adds `frame-ancestors 'none'`
 * which cannot be enforced via a meta tag. The header-based CSP provides
 * defense-in-depth alongside the meta CSP and navigation guards.
 *
 * Note: For `file://` content in packaged mode, the browser may not apply
 * response headers (there is no HTTP response). The meta CSP in index.html
 * and splash.html remains the primary enforcement for packaged builds.
 * The header CSP covers dev mode (served over HTTP) and acts as a safety
 * net for any content loaded via HTTP in production.
 */
const PRODUCTION_CSP = [
  "default-src 'self'",
  "script-src 'self' 'wasm-unsafe-eval'",
  "style-src 'self' 'unsafe-inline'",
  "connect-src 'self' http://127.0.0.1:* ws://127.0.0.1:*",
  "img-src 'self' data:",
  "worker-src 'self' blob:",
  "child-src 'self' blob:",
  "frame-ancestors 'none'",
].join('; ');

/**
 * Dev CSP relaxes script-src to allow eval AND inline scripts for Vite HMR /
 * React Fast Refresh. In dev, @vitejs/plugin-react injects an inline <script>
 * preamble (head-prepend) which requires 'unsafe-inline' to execute. This is an
 * HTTP response-header CSP, so unlike the meta CSP it governs the whole document
 * regardless of script position. Dev-only; PRODUCTION_CSP stays strict.
 */
const DEV_CSP = [
  "default-src 'self'",
  "script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval' 'unsafe-inline'",
  "style-src 'self' 'unsafe-inline'",
  "connect-src 'self' http://127.0.0.1:* ws://127.0.0.1:*",
  "img-src 'self' data:",
  "worker-src 'self' blob:",
  "child-src 'self' blob:",
  "frame-ancestors 'none'",
].join('; ');

/**
 * Harden a session with permission restrictions and CSP response headers.
 *
 * Should be called on each BrowserWindow's session before content is loaded.
 * This is defense-in-depth — the existing meta CSP and navigation guards
 * remain as primary controls.
 */
export const hardenSession = (
  windowSession: Session,
  isPackaged: boolean
): void => {
  // Skip if this session has already been hardened (e.g. splash and main
  // windows sharing session.defaultSession).
  if (hardenedSessions.has(windowSession)) {
    return;
  }
  hardenedSessions.add(windowSession);

  // Deny all permission requests except explicitly allowed ones
  windowSession.setPermissionRequestHandler(
    (_webContents, permission, callback) => {
      if (ALLOWED_PERMISSIONS.has(permission)) {
        callback(true);
        return;
      }
      console.warn(`[Security] Denied permission request: ${permission}`);
      callback(false);
    }
  );

  // Inject CSP as a response header for HTTP-served content.
  // This enables frame-ancestors enforcement which meta CSP cannot provide.
  const cspValue = isPackaged ? PRODUCTION_CSP : DEV_CSP;

  windowSession.webRequest.onHeadersReceived((details, callback) => {
    // Build response headers, removing any pre-existing CSP header regardless
    // of case (e.g. 'content-security-policy', 'Content-Security-Policy') to
    // prevent duplicate/conflicting CSP directives.
    const incoming = details.responseHeaders ?? {};
    const filtered: Record<string, string[]> = {};
    for (const [key, value] of Object.entries(incoming)) {
      if (key.toLowerCase() !== 'content-security-policy') {
        filtered[key] = value;
      }
    }
    filtered['Content-Security-Policy'] = [cspValue];

    callback({ responseHeaders: filtered });
  });
};
