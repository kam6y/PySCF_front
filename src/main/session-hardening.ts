import type { Session } from 'electron';

/**
 * Track sessions that have already been hardened and whether a port-pinned
 * CSP has been applied. Both splash and main windows share
 * session.defaultSession (no partition), and Electron's onHeadersReceived /
 * setPermissionRequestHandler are single-handler APIs (last-write-wins).
 *
 * We allow re-hardening a session when a backendPort is provided for the
 * first time so the main window can upgrade the splash's initial
 * port-less CSP to a port-pinned one (A2 fix). Only the
 * onHeadersReceived handler is re-registered in that case (permission
 * handler stays).
 */
const hardenedSessions = new WeakMap<Session, { portPinned: boolean }>();

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
 * Validate that a port number is a usable integer in the 1..65535 range.
 * Shared by both buildCsp and hardenSession so the portPinned flag and the
 * CSP directive always agree on whether a port is valid (H1 fix).
 */
export const isValidPort = (port: number | undefined): port is number =>
  port !== undefined &&
  Number.isInteger(port) &&
  port >= 1 &&
  port <= 65535;

/**
 * Build a CSP string with optional backend port scoping.
 *
 * When a backendPort is provided, connect-src is locked to that specific port
 * on 127.0.0.1 (both HTTP and WS). When omitted (e.g. splash window which
 * does not contact the backend), connect-src allows only 'self' — no loopback
 * wildcard at all.
 *
 * Dev mode additionally requires 'unsafe-eval' and 'unsafe-inline' in
 * script-src for Vite HMR / React Fast Refresh, and keeps the loopback
 * wildcard in connect-src because Vite's dev server port may differ from the
 * backend port.
 */
export const buildCsp = (isPackaged: boolean, backendPort?: number): string => {
  // E1/J7: Use the shared isValidPort helper so buildCsp and hardenSession
  // always agree on what constitutes a valid port. Malformed values (NaN, 0,
  // negative, non-integer, out-of-range) produce the "no loopback" CSP (fail-closed).
  const validPort = isValidPort(backendPort) ? backendPort : undefined;

  const scriptSrc = isPackaged
    ? "script-src 'self' 'wasm-unsafe-eval'"
    : "script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval' 'unsafe-inline'";

  let connectSrc: string;
  if (!isPackaged) {
    // Dev mode: keep wildcard so Vite HMR and arbitrary dev-server ports work
    connectSrc = "connect-src 'self' http://127.0.0.1:* ws://127.0.0.1:*";
  } else if (validPort !== undefined) {
    // Production with known backend port: pin to that port only (M-003)
    connectSrc = `connect-src 'self' http://127.0.0.1:${validPort} ws://127.0.0.1:${validPort}`;
  } else {
    // Production without backend (e.g. splash): no loopback at all
    connectSrc = "connect-src 'self'";
  }

  return [
    "default-src 'self'",
    scriptSrc,
    "style-src 'self' 'unsafe-inline'",
    connectSrc,
    "img-src 'self' data:",
    "worker-src 'self' blob:",
    "child-src 'self' blob:",
    "frame-ancestors 'none'",
    "object-src 'none'",
    "base-uri 'none'",
    "form-action 'none'",
  ].join('; ');
};

/**
 * Harden a session with permission restrictions and CSP response headers.
 *
 * Should be called on each BrowserWindow's session before content is loaded.
 * This is defense-in-depth — the existing meta CSP and navigation guards
 * remain as primary controls.
 *
 * @param backendPort - When provided, the production CSP pins connect-src to
 *   this specific port on 127.0.0.1 instead of using a wildcard (M-003).
 *   Omit for windows that do not contact the backend (e.g. splash).
 */
export const hardenSession = (
  windowSession: Session,
  isPackaged: boolean,
  backendPort?: number
): void => {
  // H1 fix: use the same validation for the portPinned flag and the CSP
  // so that calling hardenSession(session, true, NaN) does NOT record
  // portPinned:true with a broken CSP.
  const portIsValid = isValidPort(backendPort);

  const existing = hardenedSessions.get(windowSession);
  const isPortUpgrade = portIsValid && existing && !existing.portPinned;

  // Skip if this session is already fully hardened (port-pinned or no port
  // needed). Allow re-hardening when a valid backendPort is now provided for
  // the first time (A2 fix: splash hardens the shared session without a port;
  // the main window upgrades it with the port-pinned CSP).
  if (existing && !isPortUpgrade) {
    return;
  }

  // Register permission handler only on the first hardening pass — it does
  // not depend on the backend port.
  if (!existing) {
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
  }

  // Inject CSP as a response header for HTTP-served content.
  // This enables frame-ancestors enforcement which meta CSP cannot provide.
  // On a port-upgrade pass, onHeadersReceived is re-registered (last-write-
  // wins) so the new port-pinned CSP takes effect for all subsequent
  // responses on this session.
  const cspValue = buildCsp(isPackaged, backendPort);

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

  hardenedSessions.set(windowSession, {
    portPinned: portIsValid,
  });
};
