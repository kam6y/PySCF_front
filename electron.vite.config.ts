import { defineConfig } from 'electron-vite';
import react from '@vitejs/plugin-react';
import fs from 'node:fs/promises';
import path from 'node:path';
import type { Plugin as EsbuildPlugin } from 'esbuild';
import type { Plugin } from 'vite';

const projectRoot = process.cwd();
const resolveProject = (...segments: string[]): string =>
  path.resolve(projectRoot, ...segments);

const isKetcherReactEntry = (id: string): boolean => {
  const normalizedId = id.replace(/\\/g, '/');
  return normalizedId.endsWith('/node_modules/ketcher-react/dist/index.js');
};

const replaceKetcherMacromoleculesEditorImport = (code: string): string => {
  return code.replace(
    /return import\('\.\/index\.modern-[^/']+\.js'\);/,
    'return Promise.resolve({ default: function DisabledMacromoleculesEditor() { return null; } });'
  );
};

/**
 * Dev-only CSP relaxation plugin.
 *
 * In production, index.html ships with a hardened CSP:
 *   script-src 'self' 'wasm-unsafe-eval'
 *
 * During development, Vite HMR and @vitejs/plugin-react (React Fast Refresh)
 * require eval() AND an inline <script> preamble. This plugin appends
 * 'unsafe-eval' and 'unsafe-inline' to script-src ONLY when the dev server is
 * running (apply: 'serve' ensures this plugin is excluded from the build
 * pipeline entirely, so production stays strict). Kept in sync with the
 * dev response-header CSP in src/main/session-hardening.ts.
 */
const devCspRelaxPlugin = (): Plugin => ({
  name: 'dev-csp-relax',
  apply: 'serve',
  transformIndexHtml(html) {
    let result = html.replace(
      "script-src 'self' 'wasm-unsafe-eval'",
      "script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval' 'unsafe-inline'"
    );
    // H5: Fail loud if the marker was not found (CSP format changed)
    if (result === html) {
      throw new Error(
        '[dev-csp-relax] Failed to relax script-src — CSP marker not found in index.html. ' +
        'Has the meta CSP format changed?'
      );
    }

    const beforeConnect = result;
    // Add connect-src wildcards for dev-server and backend (A3: production
    // meta CSP no longer has the wildcard, so dev must inject it here)
    result = result.replace(
      "connect-src 'self'",
      "connect-src 'self' http://127.0.0.1:* ws://127.0.0.1:*"
    );
    // H5: Fail loud if the connect-src marker was not found
    if (result === beforeConnect) {
      throw new Error(
        '[dev-csp-relax] Failed to relax connect-src — CSP marker not found in index.html. ' +
        'Has the meta CSP format changed?'
      );
    }

    return result;
  },
});

/**
 * Production CSP stripping plugin (H2 fix).
 *
 * In packaged mode, the authoritative CSP is delivered as a response header
 * by the app:// protocol handler (app-protocol.ts) and session-hardening.ts.
 * The meta CSP in index.html/splash.html would create a SECOND policy that
 * intersects with the header CSP — per the CSP spec, EVERY policy must allow
 * a request for it to succeed. The meta `connect-src 'self'` (origin
 * app://renderer) does NOT allow http://127.0.0.1:<port>, so backend
 * fetch/SSE would be blocked regardless of the header CSP.
 *
 * This plugin strips the entire <meta http-equiv="Content-Security-Policy">
 * tag from HTML files during the build, leaving only the header CSP as the
 * single authoritative policy in production.
 */
const stripMetaCspPlugin = (): Plugin => ({
  name: 'strip-meta-csp',
  apply: 'build',
  transformIndexHtml(html) {
    // Match the entire <meta http-equiv="Content-Security-Policy" ...> tag
    // across multiple lines. Uses [\s\S]*? instead of [^>]* so the pattern
    // still matches if a formatter wraps the tag's attributes across lines (J10).
    // The tag may be self-closing or not.
    const metaCspPattern = /\s*<meta\s[\s\S]*?http-equiv=["']Content-Security-Policy["'][\s\S]*?\/?>\s*/i;
    const result = html.replace(metaCspPattern, '\n');
    if (result === html) {
      throw new Error(
        '[strip-meta-csp] Failed to strip meta CSP tag — pattern not found in HTML. ' +
        'Has the meta CSP tag been removed or reformatted?'
      );
    }
    return result;
  },
});

const disableKetcherMacromoleculesEditor = (): Plugin => ({
  name: 'disable-ketcher-macromolecules-editor',
  enforce: 'pre',
  transform(code, id) {
    if (!isKetcherReactEntry(id)) {
      return null;
    }

    const nextCode = replaceKetcherMacromoleculesEditorImport(code);

    if (nextCode === code) {
      return null;
    }

    return {
      code: nextCode,
      map: null,
    };
  },
});

const disableKetcherMacromoleculesEditorInOptimizeDeps = (): EsbuildPlugin => ({
  name: 'disable-ketcher-macromolecules-editor-optimize-deps',
  setup(build) {
    build.onLoad(
      { filter: /[/\\]node_modules[/\\]ketcher-react[/\\]dist[/\\]index\.js$/ },
      async args => {
        const code = await fs.readFile(args.path, 'utf8');
        return {
          contents: replaceKetcherMacromoleculesEditorImport(code),
          loader: 'js',
        };
      }
    );
  },
});

const rendererConfig = {
  root: projectRoot,
  publicDir: false,
  resolve: {
    alias: {
      assert: 'assert',
      buffer: 'buffer',
      process: 'process/browser',
      stream: 'stream-browserify',
      util: 'util',
    },
  },
  css: {
    modules: {
      generateScopedName: '[name]__[local]___[hash:base64:5]',
      localsConvention: 'camelCase' as const,
    },
  },
  define: {
    global: 'globalThis',
    'process.env.NODE_ENV': JSON.stringify(
      process.env.NODE_ENV || 'development'
    ),
  },
  optimizeDeps: {
    include: [
      'hoist-non-react-statics',
      'ketcher-core',
      'ketcher-react',
      'ketcher-standalone',
      'lodash',
      'lodash/fp',
    ],
    esbuildOptions: {
      plugins: [disableKetcherMacromoleculesEditorInOptimizeDeps()],
    },
  },
  plugins: [devCspRelaxPlugin(), stripMetaCspPlugin(), disableKetcherMacromoleculesEditor(), react()],
  build: {
    outDir: resolveProject('dist'),
    emptyOutDir: false,
    assetsDir: 'assets',
    assetsInlineLimit: 0,
    rollupOptions: {
      input: {
        index: resolveProject('index.html'),
        splash: resolveProject('splash.html'),
      },
    },
  },
};

export default defineConfig({
  main: {
    build: {
      outDir: resolveProject('dist'),
      emptyOutDir: false,
      rollupOptions: {
        input: {
          main: resolveProject('src/main.ts'),
        },
        output: {
          entryFileNames: '[name].js',
          format: 'cjs',
        },
      },
    },
  },
  preload: {
    build: {
      outDir: resolveProject('dist'),
      emptyOutDir: false,
      rollupOptions: {
        input: {
          preload: resolveProject('src/preload.ts'),
          splashPreload: resolveProject('src/splash/preload.ts'),
        },
        output: {
          entryFileNames: '[name].js',
          format: 'cjs',
        },
      },
    },
  },
  renderer: rendererConfig,
});
