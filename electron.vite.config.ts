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
 * may require eval(). This plugin appends 'unsafe-eval' to script-src ONLY
 * when the dev server is running (apply: 'serve' ensures this plugin is
 * excluded from the build pipeline entirely).
 */
const devCspRelaxPlugin = (): Plugin => ({
  name: 'dev-csp-relax',
  apply: 'serve',
  transformIndexHtml(html) {
    return html.replace(
      "script-src 'self' 'wasm-unsafe-eval'",
      "script-src 'self' 'wasm-unsafe-eval' 'unsafe-eval'"
    );
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
  plugins: [devCspRelaxPlugin(), disableKetcherMacromoleculesEditor(), react()],
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
