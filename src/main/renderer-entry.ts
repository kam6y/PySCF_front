import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { buildAppUrl, sanitizeForLog } from './app-protocol';

type MainRendererEntryParams = {
  backendPort: number;
  htmlPath: string;
  isPackaged: boolean;
  rendererUrl?: string;
};

type SplashRendererEntryParams = {
  htmlPath: string;
  isPackaged: boolean;
  rendererUrl?: string;
};

type FileRendererEntry = {
  type: 'file';
  path: string;
  query?: Record<string, string>;
};

type UrlRendererEntry = {
  type: 'url';
  url: string;
};

/**
 * Custom app:// protocol entry for packaged mode.
 * Uses the app:// scheme registered by app-protocol.ts.
 */
type AppRendererEntry = {
  type: 'app';
  url: string;
};

export type RendererEntry = FileRendererEntry | UrlRendererEntry | AppRendererEntry;

const ALLOWED_DEV_HOSTNAMES: ReadonlySet<string> = new Set([
  'localhost',
  '127.0.0.1',
]);

/**
 * Validate that a renderer URL is a safe loopback dev-server origin.
 * Accepts only http://localhost:<port> and http://127.0.0.1:<port>.
 * Rejects https, non-loopback hosts, credentials in URL, non-HTTP schemes.
 */
export const isAllowedDevRendererUrl = (rendererUrl: string): boolean => {
  try {
    const parsed = new URL(rendererUrl);

    if (parsed.protocol !== 'http:') {
      return false;
    }

    if (parsed.username || parsed.password) {
      return false;
    }

    if (!ALLOWED_DEV_HOSTNAMES.has(parsed.hostname)) {
      return false;
    }

    // Require an explicit port — a portless URL (parsed.port === '')
    // would silently target port 80, which is virtually never intentional
    // for a dev server. Note: WHATWG URL normalizes :80 to '' for http,
    // so explicit :80 is also rejected, which is the desired behavior.
    if (parsed.port === '') {
      return false;
    }

    const portNum = Number(parsed.port);
    if (!Number.isInteger(portNum) || portNum < 1 || portNum > 65535) {
      return false;
    }

    return true;
  } catch {
    return false;
  }
};

/**
 * Inert about: URLs that Electron may navigate to internally.
 * Only exact pathnames are allowed — never the entire about: scheme.
 */
const ALLOWED_ABOUT_PATHNAMES = new Set(['blank', 'srcdoc']);

/**
 * Check whether a navigation target URL is allowed given the loaded renderer entry.
 * - Inert about: URLs (about:blank, about:srcdoc) are always permitted.
 * - For url entries (dev mode): compare origins (scheme + host + port).
 * - For file entries (packaged mode): allow file: URLs whose path is within the
 *   same directory as the loaded htmlPath. Uses fileURLToPath for cross-platform
 *   correctness (Windows drive-letter paths).
 */
export const isAllowedNavigation = (
  targetUrl: string,
  rendererEntry: RendererEntry
): boolean => {
  try {
    const target = new URL(targetUrl);

    // Allow inert about: URLs (about:blank, about:srcdoc).
    // These are Electron-internal; they cannot host attacker scripts
    // and have no exploitable origin.
    if (
      target.protocol === 'about:' &&
      ALLOWED_ABOUT_PATHNAMES.has(target.pathname)
    ) {
      return true;
    }

    if (rendererEntry.type === 'url') {
      // H12: Parse entry URL separately to log a distinct config-error message
      // if the renderer entry URL itself is malformed.
      let allowed: URL;
      try {
        allowed = new URL(rendererEntry.url);
      } catch {
        console.error(
          `[Security] Config error: renderer entry URL is malformed: ${sanitizeForLog(rendererEntry.url)}`
        );
        return false;
      }
      return target.origin === allowed.origin;
    }

    if (rendererEntry.type === 'app') {
      // Packaged mode (app:// protocol): only allow app:// with matching
      // protocol, hostname, and port. We cannot use origin comparison because
      // the WHATWG URL spec returns origin "null" for non-standard schemes
      // (app:// is a custom Electron scheme), which would incorrectly match
      // other opaque-origin URLs like data: and javascript:.
      // Port check added per D4: app://renderer:1234/evil would otherwise
      // pass with only protocol+hostname comparison.
      // H12: Parse entry URL separately to distinguish config errors from
      // target URL issues.
      let allowed: URL;
      try {
        allowed = new URL(rendererEntry.url);
      } catch {
        console.error(
          `[Security] Config error: renderer entry URL is malformed: ${sanitizeForLog(rendererEntry.url)}`
        );
        return false;
      }
      return (
        target.protocol === allowed.protocol &&
        target.hostname === allowed.hostname &&
        target.port === allowed.port
      );
    }

    if (rendererEntry.type === 'file') {
      // Legacy file entry path (fallback — retained for dev-mode file entries
      // when no dev server URL is provided): only allow file: protocol.
      // file:// origins are opaque ("null"), so compare directory paths instead.
      if (target.protocol !== 'file:') {
        return false;
      }

      const allowedDir = path.dirname(rendererEntry.path);
      // Use fileURLToPath for correct conversion on all platforms (Windows drive letters, etc.)
      const targetPath = fileURLToPath(target);
      // Target path must be within the allowed directory
      const resolved = path.resolve(targetPath);
      const resolvedAllowed = path.resolve(allowedDir);
      return resolved.startsWith(resolvedAllowed + path.sep);
    }

    // J14: Exhaustive check — a future 4th RendererEntry variant would be
    // caught at compile time rather than silently falling through.
    const _exhaustive: never = rendererEntry;
    console.error(`[Security] Unknown renderer entry type: ${(_exhaustive as RendererEntry).type}`);
    return false;
  } catch {
    // Malformed URL — deny
    return false;
  }
};

const buildDevUrl = (rendererUrl: string, pathname: string): URL => {
  const baseUrl = rendererUrl.endsWith('/') ? rendererUrl : `${rendererUrl}/`;
  return new URL(pathname, baseUrl);
};

export const getMainRendererEntry = ({
  backendPort,
  htmlPath,
  isPackaged,
  rendererUrl,
}: MainRendererEntryParams): RendererEntry => {
  if (!isPackaged && rendererUrl) {
    if (isAllowedDevRendererUrl(rendererUrl)) {
      const url = buildDevUrl(rendererUrl, '/');
      url.searchParams.set('backend_port', String(backendPort));
      return {
        type: 'url',
        url: url.toString(),
      };
    }
    console.warn(
      `[Security] Rejected dev renderer URL: ${rendererUrl}. Falling back to file entry.`
    );
  }

  if (isPackaged) {
    // Packaged mode: use app:// custom protocol for defense-in-depth.
    // The app-protocol handler serves files from the dist/ directory with
    // path traversal protection and CSP headers.
    return {
      type: 'app',
      url: buildAppUrl('index.html', { backend_port: String(backendPort) }),
    };
  }

  // The app:// scheme is only registered when app.isPackaged, so using
  // it here would load an unregistered scheme and produce a blank window (B1).
  return {
    type: 'file',
    path: htmlPath,
    query: { backend_port: String(backendPort) },
  };
};

export const getSplashRendererEntry = ({
  htmlPath,
  isPackaged,
  rendererUrl,
}: SplashRendererEntryParams): RendererEntry => {
  if (!isPackaged && rendererUrl) {
    if (isAllowedDevRendererUrl(rendererUrl)) {
      return {
        type: 'url',
        url: buildDevUrl(rendererUrl, 'splash.html').toString(),
      };
    }
    console.warn(
      `[Security] Rejected dev renderer URL: ${rendererUrl}. Falling back to file entry.`
    );
  }

  if (isPackaged) {
    // Packaged mode: use app:// custom protocol (consistent with main window).
    return {
      type: 'app',
      url: buildAppUrl('splash.html'),
    };
  }

  return {
    type: 'file',
    path: htmlPath,
  };
};

/**
 * Minimal interface for the webContents methods used by navigation guards.
 * Allows test fakes to satisfy the contract without importing Electron types.
 */
export type NavigationGuardTarget = {
  on(
    event: 'will-navigate' | 'will-redirect',
    listener: (
      event: { preventDefault(): void; url?: string },
      url: string
    ) => void
  ): void;
  setWindowOpenHandler(
    handler: (details: { url: string }) => { action: 'deny' | 'allow' }
  ): void;
};

/**
 * Install will-navigate, will-redirect, and window-open guards on a webContents.
 * Extracted so both window managers share the same logic and it can be tested
 * without a real Electron BrowserWindow.
 */
export const installNavigationGuards = (
  webContents: NavigationGuardTarget,
  rendererEntry: RendererEntry
): void => {
  webContents.on('will-navigate', (event, url) => {
    if (!isAllowedNavigation(url, rendererEntry)) {
      event.preventDefault();
      // H3: sanitize attacker-controlled URL before logging (truncate + strip control chars)
      console.warn(`[Security] Blocked navigation to: ${sanitizeForLog(url)}`);
    }
  });

  // will-navigate does NOT fire for HTTP 3xx redirects; will-redirect covers those.
  webContents.on('will-redirect', (event, url) => {
    if (!isAllowedNavigation(url, rendererEntry)) {
      event.preventDefault();
      // H3: sanitize attacker-controlled URL before logging
      console.warn(`[Security] Blocked redirect to: ${sanitizeForLog(url)}`);
    }
  });

  // Deny all new window creation unconditionally.
  // External links are opened via the explicit 'open-external-url' IPC handler
  // (see ipc.ts), which applies validateExternalUrl and sender verification.
  // Routing window.open through openExternal would bypass those checks.
  webContents.setWindowOpenHandler(({ url }) => {
    // H3: sanitize attacker-controlled URL before logging
    console.warn(`[Security] Denied window.open for URL: ${sanitizeForLog(url)}`);
    return { action: 'deny' as const };
  });
};
