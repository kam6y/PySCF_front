/**
 * Custom app:// protocol for serving packaged renderer assets.
 *
 * Electron's guidance recommends custom protocols over file:// because file://
 * has special local-file behavior. If renderer JavaScript is ever compromised,
 * file:// increases the blast radius compared with an app-scoped custom protocol.
 *
 * The protocol is registered as privileged (standard, secure, supportFetchAPI)
 * before app.ready via registerSchemesAsPrivileged, then the actual handler is
 * installed after app.ready via protocol.handle.
 *
 * Security measures:
 * - Path traversal prevention (reject .., absolute escapes, %2e%2e)
 * - Only serves files under the packaged renderer directory
 * - Only GET requests are allowed
 * - Correct Content-Type headers based on file extension
 * - CSP delivered as response header from the existing buildCsp()
 */

import { protocol, net } from 'electron';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { buildCsp } from './session-hardening';

/**
 * The custom scheme name. Must match what is registered in
 * registerSchemesAsPrivileged and used in renderer entry URLs.
 */
export const APP_SCHEME = 'app';

/**
 * The hostname used in app:// URLs. A standard+secure scheme requires a host.
 * Using a single fixed hostname keeps URL handling predictable.
 */
export const APP_HOST = 'renderer';

/**
 * Construct a full app:// URL for a given asset path relative to the
 * renderer directory (e.g. 'index.html', 'assets/foo.js').
 */
export const buildAppUrl = (assetPath: string, query?: Record<string, string>): string => {
  const url = new URL(`${APP_SCHEME}://${APP_HOST}/${assetPath}`);
  if (query) {
    for (const [key, value] of Object.entries(query)) {
      url.searchParams.set(key, value);
    }
  }
  return url.toString();
};

/**
 * MIME type map for common file types served by the renderer.
 * Falls back to 'application/octet-stream' for unknown extensions.
 */
const MIME_TYPES: Readonly<Record<string, string>> = {
  '.html': 'text/html',
  '.htm': 'text/html',
  '.js': 'application/javascript',
  '.mjs': 'application/javascript',
  '.css': 'text/css',
  '.json': 'application/json',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.jpeg': 'image/jpeg',
  '.gif': 'image/gif',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
  '.webp': 'image/webp',
  '.woff': 'font/woff',
  '.woff2': 'font/woff2',
  '.ttf': 'font/ttf',
  '.otf': 'font/otf',
  '.eot': 'application/vnd.ms-fontobject',
  '.wasm': 'application/wasm',
  '.map': 'application/json',
  '.txt': 'text/plain',
  '.xml': 'application/xml',
};

const DEFAULT_MIME_TYPE = 'application/octet-stream';

export const getMimeType = (filePath: string): string => {
  const ext = path.extname(filePath).toLowerCase();
  return MIME_TYPES[ext] ?? DEFAULT_MIME_TYPE;
};

// Patterns checked against the raw URL path. Double-encoding (e.g. %252e)
// decodes to literal '%2e' which is a harmless directory name on disk — the
// final containment check (resolved.startsWith(baseDir)) catches any actual
// escape. These patterns are a fast-reject optimization, not the sole defense.
// Backslash variants (%5c / %5C) are included for defense-in-depth on
// Windows-like path interpretations, even though Node path.resolve on POSIX
// treats backslash as a literal character.
const TRAVERSAL_PATTERNS = [
  '..', '%2e%2e', '%2E%2E', '%2e.', '%2E.', '.%2e', '.%2E',
  '..%5c', '..%5C', '%5c..', '%5C..',
];

/** Maximum length of a pathname to include in log messages (D3). */
const LOG_PATH_MAX_LENGTH = 200;

/**
 * Truncate a string and strip control characters for safe logging.
 * Prevents log injection/flooding from attacker-controlled input (D3).
 */
export const sanitizeForLog = (input: string): string => {
  // Strip control characters (U+0000-U+001F, U+007F, U+0080-U+009F)
  const cleaned = input.replace(/[\x00-\x1f\x7f-\x9f]/g, '');
  if (cleaned.length <= LOG_PATH_MAX_LENGTH) {
    return cleaned;
  }
  return cleaned.slice(0, LOG_PATH_MAX_LENGTH) + '...(truncated)';
};

/**
 * Validate that a request path is safe (no traversal, stays within baseDir).
 * Returns the resolved absolute file path if valid, or null if rejected.
 */
export const validateRequestPath = (
  requestPath: string,
  baseDir: string
): string | null => {
  for (const pattern of TRAVERSAL_PATTERNS) {
    if (requestPath.includes(pattern)) {
      return null;
    }
  }

  let decodedPath: string;
  try {
    decodedPath = decodeURIComponent(requestPath);
  } catch {
    // Malformed percent-encoding — reject
    return null;
  }

  // D2: null bytes cause path truncation in native file APIs
  if (decodedPath.includes('\0')) {
    return null;
  }

  if (decodedPath.includes('..')) {
    return null;
  }

  // Strip leading slashes to make the path relative to baseDir
  const relativePath = decodedPath.replace(/^\/+/, '');

  const targetRelative = relativePath || 'index.html';

  if (path.isAbsolute(targetRelative)) {
    return null;
  }

  const resolved = path.resolve(baseDir, targetRelative);
  const normalizedBase = path.resolve(baseDir);

  // Final containment check: resolved path must be within baseDir.
  // J4: This comparison is case-sensitive, which is safe even on case-insensitive
  // macOS FS: path.resolve(baseDir, rel) always preserves baseDir's exact casing
  // as the prefix, so a case-variant URL path resolves as a literal child under
  // baseDir (e.g. /APP/index.html → baseDir/APP/index.html). No directory escape
  // is possible — at most, the FS serves the same in-dir file via case aliasing.
  if (!resolved.startsWith(normalizedBase + path.sep) && resolved !== normalizedBase) {
    return null;
  }

  return resolved;
};

/**
 * Register the app:// scheme as privileged. MUST be called before app.ready.
 *
 * - standard: true — URLs follow standard URL parsing rules (origin = app://renderer)
 * - secure: true — treated as secure context (like https), enables full web API
 * - supportFetchAPI: true — allows the Fetch API to work with app:// URLs
 */
export const registerAppScheme = (): void => {
  protocol.registerSchemesAsPrivileged([
    {
      scheme: APP_SCHEME,
      privileges: {
        standard: true,
        secure: true,
        supportFetchAPI: true,
        bypassCSP: false,
        allowServiceWorkers: false,
      },
    },
  ]);
};

/**
 * Dependencies for the app:// protocol request handler.
 * Extracted to allow unit testing with stubs (J5).
 */
export type AppProtocolHandlerDeps = {
  baseDir: string;
  isPackaged: boolean;
  getBackendPort: () => number | undefined;
  fetchImpl: (url: string) => Promise<Response>;
};

/**
 * Create the app:// protocol request handler (J5 refactor).
 *
 * The handler serves files from the given base directory with:
 * - Path traversal protection
 * - Correct MIME types
 * - CSP response headers (reusing buildCsp from session-hardening.ts)
 * - Only GET method allowed
 *
 * Exported as a factory so all branches can be unit-tested with a stubbed
 * fetchImpl and Request-like objects, without requiring a live Electron process.
 */
export const createAppProtocolRequestHandler = (
  deps: AppProtocolHandlerDeps
): ((request: { method: string; url: string }) => Response | Promise<Response>) => {
  const normalizedBaseDir = path.resolve(deps.baseDir);

  return (request: { method: string; url: string }): Response | Promise<Response> => {
    // Compute the CSP on every request so we pick up backendPort once it
    // becomes available (the handler is installed before the port is known).
    const csp = buildCsp(deps.isPackaged, deps.getBackendPort());

    if (request.method !== 'GET') {
      return new Response('Method not allowed', {
        status: 405,
        headers: { 'Content-Security-Policy': csp },
      });
    }

    let requestUrl: URL;
    try {
      requestUrl = new URL(request.url);
    } catch {
      return new Response('Bad request', {
        status: 400,
        headers: { 'Content-Security-Policy': csp },
      });
    }

    const filePath = validateRequestPath(requestUrl.pathname, normalizedBaseDir);
    if (!filePath) {
      // Truncate and sanitize the pathname for logging to prevent log
      // injection/flooding from attacker-controlled input (D3).
      const safePath = sanitizeForLog(requestUrl.pathname);
      console.warn(`[Security] Blocked app:// path traversal attempt: ${safePath}`);
      return new Response('Forbidden', {
        status: 403,
        headers: { 'Content-Security-Policy': csp },
      });
    }

    const mimeType = getMimeType(filePath);

    // Use fetchImpl to serve the file. In production this is net.fetch with
    // a file:// URL; in tests it is a stub.
    const fileUrl = pathToFileURL(filePath).href;
    return deps.fetchImpl(fileUrl).then(
      (response) => {
        if (!response.ok) {
          // J2: Log non-OK responses for debuggability in packaged builds
          console.warn(
            `[app-protocol] Non-OK response (${response.status}) for: ${sanitizeForLog(requestUrl.pathname)}`
          );
          return new Response('Not found', {
            status: 404,
            headers: { 'Content-Security-Policy': csp },
          });
        }

        // Build response with proper headers
        const headers = new Headers();
        headers.set('Content-Type', mimeType);
        headers.set('Content-Security-Policy', csp);
        headers.set('X-Content-Type-Options', 'nosniff');

        return new Response(response.body, {
          status: 200,
          headers,
        });
      },
      (error: unknown) => {
        // J2: Log sanitized error details for debuggability in packaged builds.
        // Permission errors (EACCES/EMFILE) would otherwise become silent 404s.
        // The generic 404 response to the renderer is kept to avoid leaking info.
        const errName = error instanceof Error ? error.name : 'UnknownError';
        const errMsg = error instanceof Error ? error.message : String(error);
        console.error(
          `[app-protocol] Fetch failed for ${sanitizeForLog(requestUrl.pathname)}: ${errName}: ${sanitizeForLog(errMsg)}`
        );
        return new Response('Not found', {
          status: 404,
          headers: { 'Content-Security-Policy': csp },
        });
      }
    );
  };
};

/**
 * Install the app:// protocol handler. MUST be called after app.ready.
 *
 * Thin wrapper around {@link createAppProtocolRequestHandler} that wires
 * it into Electron's protocol.handle with net.fetch as the fetchImpl.
 *
 * @param baseDir - The directory containing the packaged renderer assets (dist/)
 * @param isPackaged - Whether the app is running in packaged mode
 * @param getBackendPort - A getter that returns the current backend port, or
 *   undefined if the port is not yet known. The CSP is computed on every
 *   request using this getter so that the protocol handler automatically
 *   picks up the port once it becomes available (fixes A1 — the handler is
 *   installed before the port is known).
 */
export const installAppProtocolHandler = (
  baseDir: string,
  isPackaged: boolean,
  getBackendPort: () => number | undefined = () => undefined
): void => {
  // J17(a): Fail fast if baseDir does not contain index.html — catches
  // main/renderer outDir divergence early. Electron's patched fs.existsSync
  // works transparently with asar archives.
  const indexPath = path.join(path.resolve(baseDir), 'index.html');
  if (!fs.existsSync(indexPath)) {
    throw new Error(
      `[app-protocol] baseDir does not contain index.html: ${baseDir}. ` +
      'Check that the renderer build output directory matches the protocol handler baseDir.'
    );
  }

  const handler = createAppProtocolRequestHandler({
    baseDir,
    isPackaged,
    getBackendPort,
    fetchImpl: (url: string) => net.fetch(url),
  });

  protocol.handle(APP_SCHEME, handler as (request: GlobalRequest) => Response | Promise<Response>);
};
