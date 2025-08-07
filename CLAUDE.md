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

# Clean build directory
rimraf dist
```

The development command runs webpack in watch mode and launches Electron with automatic restart on changes using electronmon.

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
1. App.tsx manages page routing between three main sections via dropdown navigation
2. DrawMoleculePage contains the molecular visualization workflow:
   - User inputs XYZ coordinates via XYZInput component
   - Data is validated using `xyzParser.ts` utilities
   - Valid coordinates are passed to MoleculeViewer component
   - MoleculeViewer loads the molecular structure into 3Dmol.js viewer
   - User can adjust visualization styles and take screenshots

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

- The app uses a custom titlebar with `titleBarStyle: 'hidden'` and `titleBarOverlay`
- Security best practices are followed (nodeIntegration disabled, contextIsolation enabled)
- XYZ parser supports standard XYZ format with comprehensive validation
- Sample molecules (water, methane, benzene) are available for testing
- Navigation between app sections is implemented via dropdown menu (CalculationSettings, Results, DrawMolecule pages)

## Common Tasks

When adding new molecular visualization features, work with the MoleculeViewer component and 3Dmol.js APIs. For new input formats, extend the parsing utilities in `src/web/utils/`. For UI changes, follow the existing component patterns in `src/web/components/`.