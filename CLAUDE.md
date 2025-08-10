# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a PySCF Native App - an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates, retrieving molecular structures from PubChem, and visualizing 3D molecular structures using 3Dmol.js. The backend is a Python Flask server that handles all chemical computations and data management.

## Development Commands

```bash
# Development mode (builds and runs Electron with hot reload + Python backend)
npm run dev

# Production build (builds frontend and creates Python executable)
npm run build

# Package application for distribution (includes production build)
npm run package

# Clean build directory (automatically done before dev)
rimraf dist

# Individual build commands
npm run dev:webpack    # Build with webpack in development mode
npm run dev:electron   # Start Electron (requires dist files to exist)

# Python backend development (separate terminal)
cd src/python
uv run python app.py   # Start Flask API server
uv run pytest tests/   # Run Python backend tests
uv sync                # Install/update Python dependencies

# Build Python executable only
npm run build:python   # Uses PyInstaller to create standalone executable
```

## Development Workflow:

npm run dev automatically cleans dist/, builds with webpack in watch mode, and starts both the Electron app and the Python Flask server.

Webpack compiles 3 separate bundles: main process, preload script, and renderer (React app).

electronmon watches dist/**/* and restarts Electron automatically on changes.

The Python Flask server (defaults to port 5000, configurable via FLASK_RUN_PORT env var) starts automatically via the Electron main process.

The Main process includes a health check mechanism that pings /health endpoint to ensure the Python server is ready before the UI loads.

## Architecture Overview

### Electron Structure
**Main Process (src/main.ts):** Creates the BrowserWindow, manages the Python Flask subprocess (using spawn with 'uv run python' in dev, PyInstaller executable in production), and handles application lifecycle events. Includes health check mechanism that pings http://127.0.0.1:5000/health before loading UI.

**Preload Script (src/preload.ts):** Securely exposes limited Node.js/Electron APIs to the renderer process (currently minimal).

**Renderer Process (src/web/):** The React application providing the user interface.

### Python Backend (src/python/)
A Flask API server that provides endpoints for:

- PubChem integration (searching compounds, converting to XYZ)
- SMILES string conversion using RDKit  
- Running quantum chemistry calculations via PySCF
- Managing calculation files (listing, renaming, deleting)
- Health check endpoint for startup coordination

**API Port:** Defaults to 5000 (configurable via FLASK_RUN_PORT env var).

### Build System
**Frontend:** Uses Webpack with ts-loader to compile TypeScript and React code into 3 bundles: main process, preload script, and renderer.

**Backend:** Uses PyInstaller to bundle the Flask application and its dependencies into a single executable for production. This is handled by the npm run build:python script.

**Packaging:** electron-builder packages the Electron app and the Python executable (as an extraResource) into distributable formats (DMG for macOS, AppImage for Linux, NSIS for Windows).

### Core Components & State Management
**App.tsx:** Root component managing UI state (sidebar, pages) and orchestrating data flow.

**Custom Hooks (src/web/hooks/):**
- `useCalculations`: Manages the list of all calculations by fetching from the backend API. Handles creation, renaming, and deletion.
- `useActiveCalculation`: Manages which calculation is currently being viewed or edited.

**API Client (src/web/api/apiClient.ts):** Centralized module for all fetch requests to the Python backend.

**State Flow:**
1. UI interaction (e.g., clicking "delete") triggers a function in a React component
2. The function calls a method from the useCalculations hook (e.g., deleteCalculation)
3. The hook's method calls the apiClient to send a request to the Flask backend (e.g., DELETE /api/quantum/calculations/...)
4. The Flask backend performs the operation (e.g., deletes the directory)
5. The API call returns, and the hook updates its local state (e.g., removes the calculation from the list), causing the UI to re-render

## File Structure
src/
├── main.ts              # Electron main process
├── preload.ts           # Electron preload script
├── types/               # TypeScript definitions
├── web/
│   ├── api/
│   │   └── apiClient.ts # Centralized API client
│   ├── App.tsx          # Main React application
│   ├── components/      # Reusable React components
│   ├── hooks/           # Custom hooks for state management
│   ├── pages/           # Page components
│   └── ...
└── python/
    ├── app.py           # Flask API server main entry point
    ├── pyproject.toml   # Python dependencies (uv)
    ├── pubchem/         # PubChem integration modules
    ├── quantum_calc/    # PySCF calculation modules
    └── tests/           # Python unit tests

## Key API Endpoints

- `GET /health` - Health check endpoint used by main process
- `POST /api/pubchem/search` - Search PubChem by name, CID, or formula  
- `POST /api/smiles/convert` - Convert SMILES to XYZ format
- `POST /api/quantum/calculate` - Run DFT calculations using PySCF
- `GET /api/quantum/calculations` - List all calculation directories
- `GET /api/quantum/calculations/<id>` - Get calculation details
- `DELETE /api/quantum/calculations/<id>` - Delete calculation files

## Troubleshooting

**Port Conflict:** If the Flask server fails to start, the port (default 5000) might be in use. You can change it by setting the FLASK_RUN_PORT environment variable.

**Build Failures:** Ensure Python dependencies (uv sync) and Node dependencies (npm install) are up to date. PyInstaller builds can sometimes be tricky; check its logs for errors.

**Python Dependencies:** This project requires PySCF, RDKit, and geometric for quantum chemistry calculations. These have significant dependencies that may need compilation.

**SSL/HTTPS:** pubchem/client.py uses standard HTTPS requests with SSL verification enabled. No special configuration should be needed.