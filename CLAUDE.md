# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a PySCF Native App - an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates, retrieving molecular structures from PubChem database, and visualizing 3D molecular structures using 3Dmol.js.

## Development Commands

```bash
# Development mode (builds and runs Electron with hot reload + Python backend)
npm run dev

# Production build
npm run build

# Clean build directory (automatically done before dev)
rimraf dist

# Individual build commands
npm run dev:webpack    # Build with webpack in development mode
npm run dev:electron   # Start Electron (requires dist files to exist)

# Python backend development (separate terminal)
cd src/python
uv run python app.py   # Start Flask API server for PubChem integration
uv run pytest tests/   # Run Python backend tests
uv sync                # Install/update Python dependencies

# TypeScript compilation and type checking
npx tsc --noEmit       # Type check without emitting files
```

**Development Workflow:**
- `npm run dev` automatically cleans dist/, builds with webpack in watch mode, starts Electron AND Python Flask server
- Webpack compiles 3 separate bundles: main process, preload script, and renderer (React app)
- `electronmon` watches `dist/**/*` and restarts Electron automatically on changes
- `wait-on` ensures build artifacts exist before starting Electron to prevent startup errors
- Development mode includes source maps and opens DevTools in detached mode
- Python Flask server (port 5000) starts automatically and provides PubChem API integration
- SSL verification is disabled in development for PubChem HTTPS requests

## Architecture Overview

### Electron Structure
- **Main Process** (`src/main.ts`): Creates BrowserWindow with custom titlebar styling, manages Python Flask subprocess, loads the React app
- **Preload Script** (`src/preload.ts`): Currently minimal, handles secure communication between main and renderer
- **Renderer Process** (`src/web/`): React application for the UI
- **Python Backend** (`src/python/`): Flask API server providing PubChem integration, SMILES conversion, and molecular data processing

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
- **Flask**: Python web framework for backend API
- **RDKit**: Chemistry toolkit for molecular data processing and SMILES handling
- **uv**: Modern Python package manager for dependency management

### Application Flow
1. **Page Management**: App.tsx manages state for sidebar, dropdown navigation, and current page routing
2. **Navigation**: Independent sidebar toggle button and header dropdown provide access to three main pages
3. **Molecular Visualization Workflow** (CalculationSettingsPage):
   - **PubChem Integration**: User searches molecular compounds by name or CID
   - **Direct API Access**: Flask backend queries PubChem REST API directly (bypassing SSL issues)
   - **Data Retrieval**: Obtains molecular properties (IUPAC name, formula, weight, synonyms) and 3D coordinates
   - **XYZ Conversion**: Backend converts SDF data to compact XYZ format with proper chemical names
   - **UI Display**: Retrieved XYZ data automatically populates Direct XYZ Input/Edit field
   - **Manual Input**: User can directly input XYZ coordinates via XYZInput component with validation
   - **Real-time Validation**: Data is parsed and validated using `xyzParser.ts` utilities
   - **3D Visualization**: Valid coordinates are passed to MoleculeViewer component
   - **Interactive Viewer**: MoleculeViewer creates a 3Dmol.js viewer instance and loads molecular structure
   - **Style Controls**: StyleControls allow real-time adjustment of visualization appearance (atoms, bonds, surfaces)
   - **Screenshots**: Screenshots can be taken directly from the 3D viewer
4. **State Management**: All UI state (sidebar, dropdown, current page, XYZ data) managed in React hooks with centralized control

### File Structure
```
src/
├── main.ts              # Electron main process with Python subprocess management
├── preload.ts           # Electron preload script  
├── types/               # TypeScript definitions for 3Dmol.js and Electron
├── web/
│   ├── App.tsx          # Main React application with page routing
│   ├── components/      # React components (Header, Sidebar, MoleculeViewer, XYZInput, etc.)
│   ├── pages/           # Page components (CalculationSettings, Results, DrawMolecule)
│   ├── utils/           # Utility functions (XYZ parsing)
│   └── index.html       # HTML template with CSP allowing localhost:5000
└── python/              # Python Flask backend
    ├── app.py           # Flask API server main entry point
    ├── pyproject.toml   # uv project configuration with dependencies
    ├── pubchem/         # PubChem integration modules
    │   ├── __init__.py
    │   ├── client.py    # PubChem API client with SSL bypassing
    │   └── parser.py    # XYZ format parsing and conversion utilities
    ├── SMILES/          # SMILES format conversion and processing
    │   ├── __init__.py
    │   └── smiles_converter.py  # RDKit-based SMILES to 3D coordinate conversion
    └── tests/           # Python unit tests and API tests
        ├── __init__.py
        ├── test_pubchem.py      # PubChem client tests
        └── test_flask_api.py    # Flask API endpoint tests
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
- Content Security Policy (CSP) configured to allow connections to localhost:5000 for Flask API
- Python Flask subprocess automatically managed by Electron main process

### Code Organization & Standards
- TypeScript strict mode enabled with comprehensive type definitions in `src/types/`
- Custom TypeScript definitions for 3Dmol.js library and Electron APIs
- XYZ parser supports standard molecular coordinate format with validation
- Sample molecules available for testing: water, methane, benzene
- Python backend follows modern practices with uv package management
- Flask API uses REST conventions with JSON responses and proper error handling
- SSL certificate verification disabled in development for PubChem API access
- Comprehensive test coverage for both Python backend and TypeScript frontend

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
- If PubChem API fails with SSL errors, check that SSL verification is disabled in `pubchem/client.py`
- If Flask server fails to start, ensure port 5000 is not occupied (disable AirPlay Receiver on macOS)
- Python dependencies managed by uv - use `uv sync` in `src/python/` to install/update packages
- If XYZ data doesn't display properly, check that the coordinate format uses compact spacing (8.4f format)

## Common Tasks

### Frontend Development
- **New molecular visualization features**: Work with the MoleculeViewer component and 3Dmol.js APIs
- **New input formats**: Extend the parsing utilities in `src/web/utils/`
- **UI changes**: Follow the existing component patterns in `src/web/components/`
- **XYZ input modifications**: Update XYZInput component with external value control support

### Backend Development
- **PubChem API enhancements**: Modify `src/python/pubchem/client.py` for new compound properties
- **New data formats**: Extend `src/python/pubchem/parser.py` for different molecular formats
- **API endpoints**: Add new routes to `src/python/app.py` following REST conventions
- **Testing**: Add tests to `src/python/tests/` for new functionality

### Full-Stack Integration
- **New molecular databases**: Implement similar to PubChem integration pattern
- **Data flow**: Ensure proper state management between React frontend and Flask backend
- **Error handling**: Implement consistent error reporting across both frontend and backend
- **Security**: Maintain CSP policies when adding new external API connections

## PubChem Integration Usage

### Basic Usage
1. Start the application with `npm run dev` (automatically starts both Electron and Flask server)
2. Navigate to Calculation Settings page
3. Select "Get from PubChem name/ID" option
4. Enter molecular name (e.g., "water", "caffeine", "aspirin") or CID number
5. Click "Convert to XYZ" button
6. Molecular structure appears in both XYZ text field and 3D viewer

### API Endpoints
- **POST** `/api/pubchem/search`: Search compounds and retrieve XYZ coordinates
- **POST** `/api/pubchem/validate`: Validate XYZ format strings
- **POST** `/api/smiles/convert`: Convert SMILES strings to 3D XYZ coordinates using RDKit
- **GET** `/health`: Health check for Flask server

### Data Format
Retrieved XYZ data uses compact formatting:
```
3
oxidane (H2O) - CID:962
O  0.000   0.000   0.000
H  0.277   0.893   0.254
H  0.607  -0.238  -0.717
```