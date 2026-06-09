import type {
  Session,
  OnBeforeSendHeadersListenerDetails,
  BeforeSendResponse,
} from 'electron';

/**
 * Track sessions that have already had auth injection registered.
 * Both splash and main windows may share session.defaultSession (no partition),
 * and createWindow can be called more than once (e.g. activate / second-instance
 * paths in main.ts). onBeforeSendHeaders is a single-handler (last-wins) API,
 * so re-registering is harmless today but the WeakSet makes it explicit and
 * consistent with hardenSession's pattern.
 */
const registeredSessions = new WeakSet<Session>();

/**
 * Inject the backend auth token into renderer->backend requests at the network
 * layer. The token is NEVER exposed to the renderer/DOM world, so a compromised
 * renderer cannot read it (e.g. via a removed getAuthToken bridge) and exfiltrate
 * it through openExternalUrl.
 *
 * The URL filter is intentionally port-agnostic (`http://127.0.0.1/*` and
 * `http://127.0.0.1:*​/*`) so the handler fires for ALL loopback requests.
 * The authoritative port gate lives inside the handler itself: it parses the
 * port from `details.url` and injects the token ONLY when it matches
 * `backendPort`. This design is immune to Electron/Chromium match-pattern
 * port semantics, which have differed across versions; the in-handler check
 * is fully unit-testable and deterministic.
 *
 * Non-matching requests (different port, malformed URL) are passed through
 * with their original headers untouched.
 *
 * NOTE: webRequest.onBeforeSendHeaders is a SINGLE-handler (last-wins) API per
 * session. This is the only onBeforeSendHeaders registration in the app; do not
 * add another on the same session or it will silently replace this one.
 */
export function registerBackendAuthInjection(
  windowSession: Session,
  backendPort: number,
  authToken: string
): void {
  if (!Number.isInteger(backendPort) || backendPort <= 0) {
    throw new Error(
      `registerBackendAuthInjection: backendPort must be a positive integer, got ${backendPort}`
    );
  }
  if (!authToken) {
    throw new Error(
      'registerBackendAuthInjection: authToken must be a non-empty string'
    );
  }

  if (registeredSessions.has(windowSession)) {
    return;
  }
  registeredSessions.add(windowSession);

  const backendPortString = String(backendPort);

  windowSession.webRequest.onBeforeSendHeaders(
    { urls: ['http://127.0.0.1/*', 'http://127.0.0.1:*/*'] },
    (
      details: OnBeforeSendHeadersListenerDetails,
      callback: (response: BeforeSendResponse) => void
    ) => {
      let port: string | undefined;
      try {
        port = new URL(details.url).port;
      } catch {
        port = undefined;
      }

      if (port === backendPortString) {
        callback({
          requestHeaders: {
            ...details.requestHeaders,
            'X-Auth-Token': authToken,
          },
        });
      } else {
        callback({ requestHeaders: details.requestHeaders });
      }
    }
  );
}
