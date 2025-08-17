This file provides guidance to Claude when working with code in this repository.Your owner is Japanese, so you must use Japanese when reporting your plans.

Project Overview
This is a PySCF_front, an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates, retrieving molecular structures from PubChem and SMILES strings, and visualizing 3D molecular structures using 3Dmol.js. The backend is a Python Flask server that handles all chemical computations and data management using libraries like PySCF and RDKit.

The application supports various quantum chemistry calculation methods, including DFT, Hartree-Fock (HF), MP2, and TDDFT, with options for solvent effects and advanced analysis like Natural Transition Orbitals (NTO) for TDDFT.

A key feature of this project is its API-first development approach, using an OpenAPI specification as the single source of truth for the API contract between the frontend and backend. The application now uses WebSockets for real-time status updates of running calculations, providing a more efficient and responsive user experience than the previous polling-based system.

Development Philosophy
This is a development-stage application. Backward compatibility is not a concern, and breaking changes should be made freely in favor of better design and simpler code. When refactoring or improving the codebase:

Prioritize simplicity over compatibility - Remove deprecated patterns and complex fallback logic

Make breaking changes confidently - Don't hesitate to change APIs, data structures, or file formats

Clean up legacy code - Remove old implementations when better alternatives are available

Focus on the best solution - Don't compromise design quality for compatibility with older versions

This approach allows for rapid iteration and prevents technical debt accumulation during the development phase.

Development Commands

## Initial Setup

### Conda Environment (Required)
This project requires conda environment for development. The application no longer supports automatic fallback to other Python environments.

```bash
# Install Node.js dependencies
npm install

# Install Miniforge (if not already installed)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3

# Create and activate conda environment
source $HOME/miniforge3/etc/profile.d/conda.sh
conda create -y -n pyscf-env python=3.12
conda activate pyscf-env

# Install conda packages
conda install -y -c conda-forge pyscf rdkit flask geometric requests flask-cors pydantic gevent threadpoolctl

# Install pip packages
pip install flask-sock flask-pydantic datamodel-code-generator pyinstaller gevent-websocket certifi
```

**Note**: The conda environment setup is mandatory. The application will show an error dialog if the conda environment is not properly configured.

## Development Commands

```bash
# Activate conda environment (if using conda)
conda activate pyscf-env

# Development mode (generates code, builds, and runs Electron with hot reload + Python backend)
npm run dev

# Production build (generates code, builds frontend and creates Python executable)
npm run build

# Package application for distribution (includes production build)
npm run package
```
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

# Activate conda environment (if using conda)
conda activate pyscf-env

# Start Flask API server
python app.py

# Run Python backend tests
pytest tests/

# Build Python executable only (uses PyInstaller)
npm run build:python
Development Workflow
The npm run dev script is the primary command for development. It automatically:

Generates Code: Runs npm run codegen to generate TypeScript types and Python Pydantic models from src/api-spec/openapi.yaml. This ensures the frontend and backend are always in sync with the API specification.

Cleans: Cleans the dist/ directory.

Builds Frontend: Builds the frontend code (main process, preload script, and renderer) using Webpack in watch mode.

Starts Backend: Starts the Python Flask server as a subprocess from within the Electron main process.

Starts Electron: Starts the Electron application using electronmon, which watches the dist/ directory for changes and automatically restarts the app.

The Electron main process (src/main.ts) launches the Python Flask server. In development, it automatically detects and uses the conda environment (pyscf-env) if available, otherwise falls back to uv. In a packaged application, it runs the PyInstaller executable. The main process includes a health check mechanism; it continuously pings the /health endpoint of the Python server to ensure the backend is fully initialized before loading the UI. The Flask server dynamically finds a free port and reports it to the main process via stdout, ensuring no port conflicts.

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
A Flask API server that provides REST endpoints and a WebSocket interface for:

PubChem Integration: Searching compounds by name, CID, or formula and converting the data to XYZ format.

SMILES Conversion: Converting SMILES strings to 3D XYZ format using RDKit.

Quantum Chemistry Calculations: Running various calculations via PySCF. This includes DFT, Hartree-Fock (HF), MP2, and TDDFT. Calculations are executed using a ProcessPoolExecutor-based system that provides true multiprocessing parallelism, allowing multiple calculations to run simultaneously across different CPU cores.

Real-time Status Updates: A WebSocket endpoint (/ws/calculations/<id>) pushes status updates to the frontend as soon as they change on the file system.

File Management: Listing, renaming, and deleting calculation directories and files.

Health Check: An endpoint (/health) for startup coordination with the Electron main process.

Dynamic Port: The Flask server automatically finds and uses a free port on startup, communicating it back to the Electron process.

Build System
Frontend: Uses Webpack with ts-loader to compile TypeScript and React code into three separate bundles: main.js, preload.js, and the renderer's app.js.

Backend: Uses PyInstaller to bundle the Flask application and its Python dependencies into a single standalone executable for production.

Packaging: electron-builder packages the Electron app and the Python executable (included as an extraResource) into distributable formats (DMG for macOS, NSIS for Windows, AppImage for Linux).

Core Components & State Management

The application uses a modern state management architecture with clear separation of concerns:

State Management Architecture
TanStack Query: Manages all server-side state including API data fetching, caching, synchronization, and error handling. Provides automatic request deduplication, background updates, and optimistic updates.

Zustand: Manages global UI state such as the active calculation ID and staged (temporary) calculations for new calculation creation workflows.

App.tsx: The root React component orchestrates the application flow. It's significantly simplified from the previous architecture, focusing on UI logic while delegating data management to specialized hooks and stores.

Global State Store (src/web/store/):

calculationStore.ts: Zustand-based store managing UI state including activeCalculationId and stagedCalculation for handling "create from completed" workflows.

Query Hooks (src/web/hooks/):

useCalculationQueries.ts: TanStack Query-based hooks for server state management:
- useGetCalculations: Fetches and caches the calculation list
- useGetCalculationDetails: Fetches detailed calculation data with smart caching
- useStartCalculation: Mutation for starting new calculations with cache invalidation
- useDeleteCalculation: Mutation for calculation deletion with automatic cache updates
- useUpdateCalculationName: Mutation for renaming calculations
- useSearchPubChem: Mutation for PubChem API integration
- useConvertSmilesToXyz: Mutation for SMILES conversion

useCalculationSubscription: Establishes WebSocket connections for real-time status updates and directly updates TanStack Query cache for seamless state synchronization.

API Client (src/web/apiClient.ts): A centralized module containing fetch functions for all REST communication with the Python backend. It also provides a helper function to generate the correct WebSocket URL. It uses the auto-generated TypeScript types for type safety.

Type Wrappers (src/web/types/api-types.ts): This file re-exports types from the auto-generated generated-api.ts to provide convenient, application-wide aliases.

State Flow Example (Starting a Calculation)
This is an asynchronous process involving the frontend, backend, background threads, TanStack Query cache, and WebSockets.

User clicks the "Start Calculation" button in CalculationSettingsPage.

The handleStartCalculation function in App.tsx is called. It uses the useStartCalculation mutation from TanStack Query to call apiClient.startCalculation with the calculation parameters.

The Flask backend (POST /api/quantum/calculate) receives the request. It immediately:

Creates a new directory for the calculation.

Saves the initial parameters and sets the status to pending.

Submits the calculation to a ProcessPoolExecutor-managed worker process to run the appropriate PySCF calculation (DFT, HF, MP2, or TDDFT) in true parallel execution.

Returns a 202 Accepted response with the newly created (and now persistent) CalculationInstance data.

The useStartCalculation mutation's onSuccess callback automatically invalidates the calculations query cache, triggering a background refetch of the calculation list. The Zustand store clears any staged calculation and updates the activeCalculationId to the new persistent ID.

The CalculationResultsPage (or any component observing the active calculation) detects that the active calculation's status is running through the TanStack Query cache.

The useCalculationSubscription hook is triggered by the running status. It opens a WebSocket connection to ws://.../ws/calculations/<id>.

The Flask backend's WebSocket handler for that calculationId starts monitoring the status.json file on the filesystem.

As the worker process updates the status.json file (e.g., to completed or error), the WebSocket handler detects the change and pushes the complete, updated CalculationInstance data to the connected client.

The useCalculationSubscription hook's onUpdate callback receives the new data and directly updates the TanStack Query cache using queryClient.setQueryData(). This triggers automatic re-renders in all components that depend on this calculation data, providing real-time UI updates.

Once the status is completed or error, the backend WebSocket handler sends the final update and closes the connection. The TanStack Query cache retains the final calculation state for instant access.

File Structure
src/
├── api-spec/
│   └── openapi.yaml          # Single source of truth for the API
├── assets/                   # Fonts and icons
├── main.ts                   # Electron main process
├── preload.ts                # Electron preload script
├── types/                    # Global TypeScript definitions
├── web/
│   ├── apiClient.ts          # Centralized API client for the frontend
│   ├── App.tsx               # Main React application component
│   ├── App.css               # Global styles for the React app
│   ├── components/           # Reusable React components
│   ├── data/                 # Static data (e.g., atomic radii)
│   ├── hooks/                # Custom React hooks for state management
│   │   ├── useCalculationQueries.ts    # TanStack Query hooks for server state
│   │   └── useCalculationSubscription.ts # Real-time updates via WebSocket
│   ├── pages/                # Page components (one per view)
│   ├── store/                # Global state management
│   │   └── calculationStore.ts         # Zustand store for UI state
│   └── types/
│       ├── api-types.ts      # Convenience wrapper for generated types
│       └── generated-api.ts  # (auto-generated) TypeScript types from OpenAPI spec
└── python/
    ├── app.py                # Flask API server main entry point (includes WebSocket handler)
    ├── pyproject.toml        # Python dependencies (conda environment recommended)
    ├── generated_models.py   # (auto-generated) Pydantic models from OpenAPI spec
    ├── pubchem/              # PubChem integration modules
    ├── SMILES/               # SMILES conversion modules
    ├── quantum_calc/         # PySCF calculation logic and file management
    │   ├── base_calculator.py    # Abstract base class for all calculators
    │   ├── dft_calculator.py     # DFT calculation logic
    │   ├── hf_calculator.py      # Hartree-Fock calculation logic
    │   ├── mp2_calculator.py     # MP2 calculation logic
    │   ├── tddft_calculator.py   # TDDFT calculation logic and NTO analysis
    │   ├── exceptions.py         # Custom exceptions for calculations
    │   ├── file_manager.py       # Handles reading/writing calculation files
    │   ├── process_manager.py    # ProcessPoolExecutor-based parallel calculation manager
    │   └── solvent_effects.py    # Logic for applying solvent models (PCM, ddCOSMO)
    └── tests/                # Python unit tests
Key API Endpoints
GET /health: Health check endpoint used by the Electron main process during startup.

POST /api/pubchem/search: Search PubChem by name, CID, or formula.

POST /api/smiles/convert: Convert a SMILES string to XYZ format.

POST /api/pubchem/validate: Validate an XYZ format string.

POST /api/quantum/calculate: Asynchronously starts a quantum chemistry calculation (DFT, HF, MP2, TDDFT) using ProcessPoolExecutor for true parallel execution.

GET /api/quantum/calculations: Lists all saved calculation directories and provides information about active calculations.

GET /api/quantum/calculations/<id>: Gets detailed information and results for a specific calculation.

PUT /api/quantum/calculations/<id>: Updates a calculation's metadata (e.g., renames it).

POST /api/quantum/calculations/<id>/cancel: Cancels a running calculation if possible.

DELETE /api/quantum/calculations/<id>: Deletes a calculation and all its associated files, canceling it first if running.

GET /api/quantum/status: Gets status information about the calculation system including process pool usage and active calculations.

WS /ws/calculations/<id>: WebSocket endpoint for real-time status updates of a running calculation.

Parallel Processing Architecture
The application uses a ProcessPoolExecutor-based system for quantum chemistry calculations to achieve true multiprocessing parallelism:

Process Pool Management: The CalculationProcessManager class (src/python/quantum_calc/process_manager.py) manages a pool of worker processes equal to the system's CPU count. This allows multiple calculations to run simultaneously, fully utilizing available CPU cores.

Worker Process Isolation: Each calculation runs in a completely separate Python process, eliminating the Global Interpreter Lock (GIL) limitation that would prevent true parallelism with threading-based approaches.

Resource Management: The process pool automatically starts when the first calculation is submitted and shuts down cleanly when the Flask server terminates. Proper signal handling ensures no zombie processes are left behind.

Calculation Lifecycle: When a calculation is submitted via POST /api/quantum/calculate, it's queued to the process pool. The worker process updates the status.json file directly, which the WebSocket handler monitors for real-time status updates to connected clients.

Cancellation Support: Running calculations can be cancelled via POST /api/quantum/calculations/<id>/cancel. The process manager attempts to cancel the Future and updates the calculation status accordingly.

Error Handling: Process-level errors are properly captured and reported back through the filesystem-based status system, maintaining consistency with the existing WebSocket notification mechanism.

Troubleshooting
Backend Server Fails to Start: The Flask server is designed to find a free port automatically. If it still fails, ensure that no firewall is blocking local network communication and that the Python environment (conda activate pyscf-env) is correctly set up.

Build Failures: Ensure Python dependencies (conda environment or fallback uv environment) and Node dependencies (npm install) are up to date. PyInstaller builds can be sensitive; check its logs in the build/pyinstaller directory for errors.

Python Dependencies: This project requires PySCF, RDKit, and geometric. These packages have significant scientific dependencies that may require system-level libraries or compilation. Ensure your Python environment can build them.

Process Pool Issues: If calculations appear to hang or the process pool becomes unresponsive, check the /api/quantum/status endpoint for pool health. The process manager automatically handles cleanup, but in extreme cases, restarting the Flask server will reset the process pool.

SSL/HTTPS: The pubchem/client.py module uses standard HTTPS requests via the requests library with SSL verification enabled. No special configuration should be needed for it to work.