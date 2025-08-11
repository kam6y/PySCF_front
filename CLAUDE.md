CLAUDE.md
This file provides guidance to Claude when working with code in this repository.Your owner is Japanese, so you must use Japanese when reporting your plans.

Project Overview
This is a PySCF_front, an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates, retrieving molecular structures from PubChem and SMILES strings, and visualizing 3D molecular structures using 3Dmol.js. The backend is a Python Flask server that handles all chemical computations and data management using libraries like PySCF and RDKit.

A key feature of this project is its API-first development approach, using an OpenAPI specification as the single source of truth for the API contract between the frontend and backend.

## Development Philosophy

**This is a development-stage application.** Backward compatibility is not a concern, and breaking changes should be made freely in favor of better design and simpler code. When refactoring or improving the codebase:

- **Prioritize simplicity over compatibility** - Remove deprecated patterns and complex fallback logic
- **Make breaking changes confidently** - Don't hesitate to change APIs, data structures, or file formats
- **Clean up legacy code** - Remove old implementations when better alternatives are available
- **Focus on the best solution** - Don't compromise design quality for compatibility with older versions

This approach allows for rapid iteration and prevents technical debt accumulation during the development phase.

Development Commands
Bash

# Install Node.js and Python dependencies
npm install
cd src/python
uv sync
cd ../..

# Development mode (generates code, builds, and runs Electron with hot reload + Python backend)
npm run dev

# Production build (generates code, builds frontend and creates Python executable)
npm run build

# Package application for distribution (includes production build)
npm run package
Individual Commands
Bash

# Clean build directory
npm run clean

# Generate TypeScript types and Python models from OpenAPI spec
npm run codegen

# Build frontend with webpack in development mode
npm run dev:webpack

# Start Electron (requires dist files to exist)
npm run dev:electron

# --- Python Backend Development (in a separate terminal) ---
cd src/python

# Start Flask API server
uv run python app.py

# Run Python backend tests
uv run pytest tests/

# Build Python executable only (uses PyInstaller)
npm run build:python
Development Workflow
The npm run dev script is the primary command for development. It automatically:

Generates Code: Runs npm run codegen to generate TypeScript types and Python Pydantic models from src/api-spec/openapi.yaml. This ensures the frontend and backend are always in sync with the API specification.

Cleans: Cleans the dist/ directory.

Builds Frontend: Builds the frontend code (main process, preload script, and renderer) using Webpack in watch mode.

Starts Backend: Starts the Python Flask server as a subprocess from within the Electron main process.

Starts Electron: Starts the Electron application using electronmon, which watches the dist/ directory for changes and automatically restarts the app.

The Electron main process (src/main.ts) launches the Python Flask server. In development, it uses uv run python app.py. In a packaged application, it runs the PyInstaller executable. The main process includes a health check mechanism; it continuously pings the /health endpoint of the Python server to ensure the backend is fully initialized before loading the UI. The Flask server dynamically finds a free port and reports it to the main process via stdout, ensuring no port conflicts.

Architecture Overview
API-First Development with OpenAPI
The single source of truth for the API is src/api-spec/openapi.yaml.

The npm run codegen command (run automatically with dev and build) uses this file to generate:

Python Pydantic Models (src/python/generated_models.py) via datamodel-code-generator, ensuring type-safe request/response handling in the Flask backend.

TypeScript Type Definitions (src/web/types/generated-api.ts) via openapi-typescript, ensuring the frontend API client is always synchronized with the backend.

This approach prevents mismatches between the frontend and backend interfaces.

Electron Structure
Main Process (src/main.ts): Creates the BrowserWindow, manages the Python Flask subprocess, handles application lifecycle events, and implements the backend health check. It passes the dynamically assigned Flask port number to the renderer process.

Preload Script (src/preload.ts): Securely exposes specific Electron APIs to the renderer process using contextBridge. Currently, it's used to pass the Flask server's port number.

Renderer Process (src/web/): The React (TypeScript) application that provides the entire user interface.

Python Backend (src/python/)
A Flask API server that provides endpoints for:

PubChem Integration: Searching compounds by name, CID, or formula and converting the data to XYZ format.

SMILES Conversion: Converting SMILES strings to 3D XYZ format using RDKit.

Quantum Chemistry Calculations: Running DFT calculations via PySCF. Calculations are executed in background threads to keep the API responsive.

File Management: Listing, renaming, and deleting calculation directories and files.

Health Check: An endpoint (/health) for startup coordination with the Electron main process.

Dynamic Port: The Flask server automatically finds and uses a free port on startup, communicating it back to the Electron process.

Build System
Frontend: Uses Webpack with ts-loader to compile TypeScript and React code into three separate bundles: main.js, preload.js, and the renderer's app.js.

Backend: Uses PyInstaller to bundle the Flask application and its Python dependencies into a single standalone executable for production.

Packaging: electron-builder packages the Electron app and the Python executable (included as an extraResource) into distributable formats (DMG for macOS, NSIS for Windows, AppImage for Linux).

Core Components & State Management
App.tsx: The root React component. It manages the overall UI state, such as the active page and sidebar visibility, and orchestrates data flow between hooks and components.

Custom Hooks (src/web/hooks/):

useCalculations: Manages the list of all calculations. It fetches the list from the backend, handles creating new calculations, and triggers API calls for renaming and deletion.

useActiveCalculation: Manages the currently selected calculation. It fetches detailed data for the active calculation and caches it.

useCalculationPolling: Manages polling for the status of a running calculation. When a calculation is started, this hook periodically requests updates from the backend until the status is completed or error.

API Client (src/web/apiClient.ts): A centralized module containing fetch functions for all communication with the Python backend. It uses the auto-generated TypeScript types for type safety.

Type Wrappers (src/web/types/api-types.ts): This file re-exports types from the auto-generated generated-api.ts to provide convenient, application-wide aliases.

State Flow Example (Starting a Calculation)
This is an asynchronous process involving the frontend, backend, and background threads.

User clicks the "Start Calculation" button in CalculationSettingsPage.

The handleStartCalculation function in App.tsx is called. It calls apiClient.startCalculation with the calculation parameters.

The Flask backend (POST /api/quantum/calculate) receives the request. It immediately:

Creates a new directory for the calculation.

Saves the initial parameters and sets the status to pending.

Starts a new background thread to run the PySCF calculation.

Returns a 202 Accepted response with the newly created (and now persistent) CalculationInstance data.

The useCalculations hook updates the list of calculations, replacing the temporary client-side ID with the new persistent ID from the server. The useActiveCalculation hook updates the active ID.

The CalculationResultsPage detects that the active calculation's status is running and triggers the useCalculationPolling hook.

The useCalculationPolling hook starts making periodic calls to GET /api/quantum/calculations/<id>.

The hook's onUpdate callback updates the central state in App.tsx each time it receives new status information. This causes the UI to reflect the current state (e.g., in the sidebar).

When the backend thread finishes, it updates the status.json on the filesystem to completed or error. The next polling request retrieves this final status, and the hook stops polling.

File Structure
src/
├── api-spec/
│   └── openapi.yaml          # Single source of truth for the API
├── main.ts                   # Electron main process
├── preload.ts                # Electron preload script
├── types/                    # Global TypeScript definitions
├── web/
│   ├── apiClient.ts          # Centralized API client for the frontend
│   ├── App.tsx               # Main React application component
│   ├── components/           # Reusable React components
│   ├── hooks/                # Custom React hooks for state management
│   ├── pages/                # Page components
│   └── types/
│       ├── api-types.ts      # Convenience wrapper for generated types
│       └── generated-api.ts  # (auto-generated) TypeScript types from OpenAPI spec
└── python/
    ├── app.py                # Flask API server main entry point
    ├── pyproject.toml        # Python dependencies (managed by uv)
    ├── generated_models.py   # (auto-generated) Pydantic models from OpenAPI spec
    ├── pubchem/              # PubChem integration modules
    ├── SMILES/               # SMILES conversion modules
    ├── quantum_calc/         # PySCF calculation logic and file management
    └── tests/                # Python unit tests
Key API Endpoints
GET /health: Health check endpoint used by the Electron main process during startup.

POST /api/pubchem/search: Search PubChem by name, CID, or formula.

POST /api/smiles/convert: Convert a SMILES string to XYZ format.

POST /api/pubchem/validate: Validate an XYZ format string.

POST /api/quantum/calculate: Asynchronously starts a DFT calculation in a background thread.

GET /api/quantum/calculations: Lists all saved calculation directories.

GET /api/quantum/calculations/<id>: Gets detailed information and results for a specific calculation.

PUT /api/quantum/calculations/<id>: Updates a calculation's metadata (e.g., renames it).

DELETE /api/quantum/calculations/<id>: Deletes a calculation and all its associated files.

Troubleshooting
Backend Server Fails to Start: The Flask server is designed to find a free port automatically. If it still fails, ensure that no firewall is blocking local network communication and that the Python environment (uv sync) is correctly set up.

Build Failures: Ensure Python dependencies (uv sync in src/python) and Node dependencies (npm install) are up to date. PyInstaller builds can be sensitive; check its logs in the build/pyinstaller directory for errors.

Python Dependencies: This project requires PySCF, RDKit, and geometric. These packages have significant scientific dependencies that may require system-level libraries or compilation. Ensure your Python environment can build them.

SSL/HTTPS: The pubchem/client.py module uses standard HTTPS requests via the requests library with SSL verification enabled. No special configuration should be needed for it to work.