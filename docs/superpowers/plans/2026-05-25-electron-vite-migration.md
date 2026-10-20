# Electron Vite Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the webpack/electronmon development and build pipeline with electron-vite while preserving current Electron main, preload, renderer, splash, Python backend startup, and electron-builder packaging behavior.

**Architecture:** Use `electron.vite.config.ts` as the single bundler config for main, preload, and renderer targets. Keep source entry points in their current locations, move HTML files to Vite entry documents at repository root, and update window loaders so development uses `ELECTRON_RENDERER_URL` while packaged builds load built HTML files. Keep OpenAPI codegen explicit through `npm run codegen` and `npm run dev:codegen`.

**Tech Stack:** Electron 37, React 19, TypeScript 5.8, electron-vite, Vite React plugin, electron-builder, npm scripts.

---

### Task 1: Dependency And Script Migration

**Files:**
- Modify: `package.json`
- Modify: `package-lock.json`

- [x] Add `electron-vite` and `@vitejs/plugin-react` as dev dependencies.
- [x] Replace the old webpack/electronmon dev scripts with `electron-vite dev`.
- [x] Make `npm run dev` run `electron-vite dev` without automatic codegen.
- [x] Add `npm run dev:codegen` for `npm run codegen && npm run dev`.
- [x] Replace `build:webpack` with `build:electron`.
- [x] Keep `build`, `build:linux`, `package`, and Docker scripts on the same high-level flow.
- [x] Remove webpack-only dev dependencies once the Vite build is passing.

### Task 2: Vite Entry And Config

**Files:**
- Create: `electron.vite.config.ts`
- Create: `index.html`
- Create: `splash.html`
- Modify: `src/web/index.html`
- Modify: `src/splash/splash.html`

- [x] Configure main entry `src/main.ts` to output `dist/main.js`.
- [x] Configure preload entries `src/preload.ts` and `src/splash/preload.ts` to output `dist/preload.js` and `dist/splashPreload.js`.
- [x] Configure renderer entries `index.html` and `splash.html` to output `dist/index.html` and `dist/splash.html`.
- [x] Preserve CSS modules, asset URLs, React transform, and Ketcher macromolecules replacement.
- [x] Move script tags into Vite HTML entry files.

### Task 3: Window Runtime Loading

**Files:**
- Modify: `src/main/window-manager.ts`
- Modify: `src/main/splash-window-manager.ts`

- [x] In development, load the main renderer with `ELECTRON_RENDERER_URL` plus the existing `backend_port` query.
- [x] In production, load `dist/index.html` as before.
- [x] In development, load the splash renderer from `ELECTRON_RENDERER_URL/splash.html`.
- [x] In production, load `dist/splash.html` as before.
- [x] Keep preload filenames unchanged so existing security settings remain stable.

### Task 4: Packaging And Documentation

**Files:**
- Modify: `package.json`
- Modify: `README.md`
- Modify: `.claude/docs/DESIGN.md`

- [x] Keep `main` pointing to `./dist/main.js`.
- [x] Keep electron-builder packaging inputs compatible with `dist/**/*`.
- [x] Update English and Japanese README setup/codegen guidance.
- [x] Record the electron-vite migration decision in the design document.

### Task 5: Verification

**Commands:**
- `npm run typecheck`
- `npm run build:electron`
- `npm run test:main`
- `npm run test:web`
- `npm run codegen:check` only if OpenAPI artifacts were touched

- [x] Run focused TypeScript checks and Electron main/web tests.
- [x] Run the electron-vite production build.
- [ ] If dev server verification is feasible, start `npm run dev` and confirm the splash and main windows reach the expected built entry points.
