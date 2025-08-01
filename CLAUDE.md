# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a PySCF molecular visualization desktop application built with React + TypeScript + Electron. It provides a 3D molecular structure viewer using 3Dmol.js, allowing users to input XYZ coordinates and visualize molecular structures with various styling options.

## Architecture

### Multi-Process Structure
- **Main Process** (`src/main.ts`): Electron's main process that manages application lifecycle and creates browser windows
- **Preload Script** (`src/preload.ts`): Security bridge between main and renderer processes  
- **Renderer Process** (`src/web/`): React application that runs in the browser window

### Build System
The project uses webpack with three separate configurations defined in `webpack.config.ts`:
1. **Main Process Bundle**: `electron-main` target for Node.js environment
2. **Preload Bundle**: `electron-preload` target with restricted APIs
3. **Renderer Bundle**: `web` target for browser environment with React

### Key Technologies
- **Electron 37+**: Desktop app framework
- **React 19+**: UI library with modern JSX transform
- **TypeScript**: Strict type checking enabled
- **Webpack 5**: Module bundler with hot reloading
- **3Dmol.js**: 3D molecular visualization library
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
1. Clears the `dist/` directory with `rimraf`
2. Runs webpack in watch mode for all three bundles in parallel
3. Waits for initial build completion with `wait-on`
4. Starts Electron with `electronmon` for automatic restarts

Hot reloading behavior:
- **Main process changes**: Full Electron app restart
- **Renderer process changes**: Browser window reload  
- **Preload script changes**: Full Electron app restart

## Application Structure

### Core Components
- **App.tsx** (`src/web/App.tsx`): Main application container managing state and component coordination
- **MoleculeViewer** (`src/web/components/MoleculeViewer.tsx`): 3D molecular visualization component using 3Dmol.js
- **XYZInput** (`src/web/components/XYZInput.tsx`): Input component for XYZ coordinate data
- **StyleControls** (`src/web/components/StyleControls.tsx`): UI controls for molecular styling options

### State Management
The application uses React's built-in state management:
- `hasValidMolecule`: Boolean tracking if valid XYZ data is loaded
- `currentXYZ`: String storing current XYZ coordinate data
- Component refs for imperative API calls to MoleculeViewer

### 3Dmol.js Integration
- Dynamic imports used for proper library loading
- Custom TypeScript definitions in `src/types/3dmol.d.ts` and `src/types/3dmol-build.d.ts`
- Viewer lifecycle managed through React refs and useEffect hooks
- Support for molecular styling, zoom-to-fit, and screenshot functionality

## Security Considerations

- Uses `target: "web"` instead of `electron-renderer` for the renderer process
- CSS is extracted to separate files rather than inlined via `style-loader`
- Content Security Policy is set in the HTML template
- Preload script provides controlled API bridge between processes
- `fsevents` external dependency for macOS compatibility

## File Structure Notes

- All build output goes to `dist/` (gitignored)
- Main process entry expects `dist/index.html` to exist
- Preload script is referenced as `preload.js` in the main process
- React app mounts to `#root` element in the HTML template
- Assets are placed in `dist/assets/` with original names
- Type definitions stored in `src/types/` directory

## TypeScript Configuration

- Uses modern ES modules with bundler resolution
- React JSX transform (no need to import React in JSX files)
- Strict mode enabled for all TypeScript checking
- Special CommonJS config for `ts-node` to handle webpack config execution
- Custom type roots include both `node_modules/@types` and `src/types`

## Molecular Visualization Features

- XYZ coordinate parsing and validation (`src/web/utils/xyzParser.ts`)
- Real-time 3D molecular structure rendering
- Interactive styling controls (stick, ball-and-stick, sphere representations)
- Zoom-to-fit functionality
- Screenshot export capability
- Error handling for invalid molecular data