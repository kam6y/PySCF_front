# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a PySCF Native App - an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates and visualizing 3D molecular structures using 3Dmol.js.

## Development Commands

```bash
# Development mode (builds and runs Electron with hot reload)
npm run dev

# Production build
npm run build

# Clean build directory (automatically done before dev)
rimraf dist

# Individual build commands
npm run dev:webpack    # Build with webpack in development mode
npm run dev:electron   # Start Electron (requires dist files to exist)
```

**Development Workflow:**
- `npm run dev` automatically cleans dist/, builds with webpack in watch mode, and starts Electron
- Webpack compiles 3 separate bundles: main process, preload script, and renderer (React app)
- `electronmon` watches `dist/**/*` and restarts Electron automatically on changes
- `wait-on` ensures build artifacts exist before starting Electron to prevent startup errors
- Development mode includes source maps and opens DevTools in detached mode

## Architecture Overview

### Electron Structure
- **Main Process** (`src/main.ts`): Creates BrowserWindow with custom titlebar styling, loads the React app
- **Preload Script** (`src/preload.ts`): Currently minimal, handles secure communication between main and renderer
- **Renderer Process** (`src/web/`): React application for the UI

### Build System
- Uses Webpack with TypeScript compilation
- Three separate build targets: main process, preload script, and renderer (React app)
- Development mode includes source maps and file watching
- Security-focused configuration (no style-loader, electron-renderer target avoided)
- Uses `electronmon` for automatic Electron restart during development
- CSS is extracted to separate files using `MiniCssExtractPlugin` for security
- Assets (fonts, images) are processed as resources and placed in `dist/assets/`

### Core Components
- **App.tsx**: Main application component with page routing, sidebar/dropdown state management
- **Pages**: Three main app sections (CalculationSettingsPage, CalculationResultsPage, DrawMoleculePage)
- **MoleculeViewer**: React component wrapping 3Dmol.js for 3D molecular visualization (in DrawMoleculePage)
- **XYZInput**: Component for inputting and validating XYZ coordinate data
- **StyleControls**: Controls for adjusting molecular visualization styles
- **Header/Sidebar/DropdownMenu**: Navigation components with page switching functionality

### Key Libraries
- **3Dmol.js**: For 3D molecular visualization and rendering
- **React 19**: UI framework
- **TypeScript**: Type safety throughout the codebase

### Application Flow
1. **Page Management**: App.tsx manages state for sidebar, dropdown navigation, and current page routing
2. **Navigation**: Independent sidebar toggle button and header dropdown provide access to three main pages
3. **Molecular Visualization Workflow** (DrawMoleculePage):
   - User inputs XYZ coordinates via XYZInput component with validation
   - Data is parsed and validated using `xyzParser.ts` utilities
   - Valid coordinates are passed to MoleculeViewer component
   - MoleculeViewer creates a 3Dmol.js viewer instance and loads molecular structure
   - StyleControls allow real-time adjustment of visualization appearance (atoms, bonds, surfaces)
   - Screenshots can be taken directly from the 3D viewer
4. **State Management**: All UI state (sidebar, dropdown, current page) managed in App.tsx with React hooks

### File Structure
```
src/
├── main.ts              # Electron main process
├── preload.ts           # Electron preload script  
├── types/               # TypeScript definitions for 3Dmol.js and Electron
├── web/
│   ├── App.tsx          # Main React application with page routing
│   ├── components/      # React components (Header, Sidebar, MoleculeViewer, etc.)
│   ├── pages/           # Page components (CalculationSettings, Results, DrawMolecule)
│   └── utils/           # Utility functions (XYZ parsing)
```

### Type Definitions
Custom TypeScript definitions are provided for:
- 3Dmol.js library (`3dmol.d.ts`, `3dmol-build.d.ts`)
- Electron APIs (`electron.d.ts`)

## Development Notes

### Electron Security & Configuration
- Custom titlebar with `titleBarStyle: 'hidden'` and transparent `titleBarOverlay` (40px height)
- Security best practices: `nodeIntegration: false`, `contextIsolation: true`
- DevTools automatically opens in detached mode during development

### Code Organization & Standards
- TypeScript strict mode enabled with comprehensive type definitions in `src/types/`
- Custom TypeScript definitions for 3Dmol.js library and Electron APIs
- XYZ parser supports standard molecular coordinate format with validation
- Sample molecules available for testing: water, methane, benzene

### Build System Details
- Webpack configuration includes Japanese comments and prioritizes security
- Three separate build targets prevent security vulnerabilities
- CSS extracted to separate files (no style-loader for security)
- Assets processed as resources and placed in `dist/assets/`
- `fsevents` externalized for macOS compatibility

### UI Architecture
- Independent sidebar toggle button with custom SVG icons
- Dropdown navigation between three main app sections
- Responsive layout that adapts when sidebar is open/closed
- Footer displays technology credits (3Dmol.js, React, TypeScript)

## Troubleshooting

- If build fails on macOS, the `fsevents` external in webpack config provides a workaround
- Electron requires both `dist/index.html` and `dist/main.js` to exist before starting
- Source maps are essential for development mode to avoid "Uncaught EvalError" in DevTools

## Common Tasks

When adding new molecular visualization features, work with the MoleculeViewer component and 3Dmol.js APIs. For new input formats, extend the parsing utilities in `src/web/utils/`. For UI changes, follow the existing component patterns in `src/web/components/`.