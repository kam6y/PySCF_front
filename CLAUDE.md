# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a React + TypeScript + Electron desktop application boilerplate. The project uses a multi-process architecture with webpack bundling for both main and renderer processes.

## Architecture

### Multi-Process Structure
- **Main Process** (`src/main.ts`): Electron's main process that manages application lifecycle and creates browser windows
- **Preload Script** (`src/preload.ts`): Security bridge between main and renderer processes
- **Renderer Process** (`src/web/`): React application that runs in the browser window

### Build System
The project uses webpack with three separate configurations:
1. **Main Process Bundle**: `electron-main` target for Node.js environment
2. **Preload Bundle**: `electron-preload` target with restricted APIs
3. **Renderer Bundle**: `web` target for browser environment with React

### Key Technologies
- **Electron 37+**: Desktop app framework
- **React 19+**: UI library with modern JSX transform
- **TypeScript**: Strict type checking enabled
- **Webpack 5**: Module bundler with hot reloading
- **CSS Modules**: Extracted CSS files (not inlined for security)

## Development Commands

```bash
# Start development with hot reload for both main and renderer processes
npm run dev

# Build for production
npm run build

# Individual webpack processes (used internally by npm run dev)
npm run dev:webpack  # Webpack watch mode
npm run dev:electron # Electron with file watching
```

## Development Workflow

The `npm run dev` command:
1. Clears the `dist/` directory
2. Runs webpack in watch mode for all three bundles
3. Waits for initial build completion
4. Starts Electron with `electronmon` for automatic restarts

Hot reloading behavior:
- **Main process changes**: Full Electron app restart
- **Renderer process changes**: Browser window reload
- **Preload script changes**: Full Electron app restart

## Security Considerations

- Uses `target: "web"` instead of `electron-renderer` for the renderer process
- CSS is extracted to separate files rather than inlined via `style-loader`
- Content Security Policy is set in the HTML template
- Preload script provides controlled API bridge between processes

## File Structure Notes

- All build output goes to `dist/` (gitignored)
- Main process entry expects `dist/index.html` to exist
- Preload script is referenced as `preload.js` in the main process
- React app mounts to `#root` element in the HTML template
- Assets are placed in `dist/assets/` with original names

## TypeScript Configuration

- Uses modern ES modules with bundler resolution
- React JSX transform (no need to import React in JSX files)
- Strict mode enabled for all TypeScript checking
- Special CommonJS config for `ts-node` to handle webpack config execution