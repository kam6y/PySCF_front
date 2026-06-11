import path from 'node:path';
import { fileURLToPath } from 'node:url';

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

export type RendererEntry = FileRendererEntry | UrlRendererEntry;

// --- SEC-002: Dev renderer URL validation ---

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

    // Only plain http allowed (not https — dev server is local)
    if (parsed.protocol !== 'http:') {
      return false;
    }

    // Reject credentials in URL (user:pass@host)
    if (parsed.username || parsed.password) {
      return false;
    }

    // Only loopback hostnames
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
    // Malformed URL → reject
    return false;
  }
};

// --- SEC-001: Navigation origin validation ---

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
      // Dev mode: origins must match exactly
      const allowed = new URL(rendererEntry.url);
      return target.origin === allowed.origin;
    }

    // Packaged mode (file entry): only allow file: protocol
    // file:// origins are opaque ("null"), so compare directory paths instead
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
  } catch {
    // Malformed URL → deny
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

  return {
    type: 'file',
    path: htmlPath,
    query: {
      backend_port: String(backendPort),
    },
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

  return {
    type: 'file',
    path: htmlPath,
  };
};

// --- SEC-001: Shared navigation guard wiring ---

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
  // Deny same-window navigation to untrusted origins
  webContents.on('will-navigate', (event, url) => {
    if (!isAllowedNavigation(url, rendererEntry)) {
      event.preventDefault();
      console.warn(`[Security] Blocked navigation to: ${url}`);
    }
  });

  // will-navigate does NOT fire for HTTP 3xx redirects; will-redirect covers those.
  webContents.on('will-redirect', (event, url) => {
    if (!isAllowedNavigation(url, rendererEntry)) {
      event.preventDefault();
      console.warn(`[Security] Blocked redirect to: ${url}`);
    }
  });

  // Deny all new window creation unconditionally.
  // External links are opened via the explicit 'open-external-url' IPC handler
  // (see ipc.ts), which applies validateExternalUrl and sender verification.
  // Routing window.open through openExternal would bypass those checks.
  webContents.setWindowOpenHandler(({ url }) => {
    console.warn(`[Security] Denied window.open for URL: ${url}`);
    return { action: 'deny' as const };
  });
};
