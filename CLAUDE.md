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
This project requires conda environment for development. The application now features automated environment setup and verification tools.

#### Quick Setup (Recommended)
```bash
# Install Node.js dependencies
npm install

# Automated environment setup (handles all conda setup)
npm run setup-env

# Verify environment health
npm run verify-env
```

#### Manual Setup
```bash
# Install Node.js dependencies
npm install

# Install Miniforge (if not already installed)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3

# Create conda environment from environment.yml (includes all dependencies)
source $HOME/miniforge3/etc/profile.d/conda.sh
conda env create -f .github/environment.yml

# Activate the environment
conda activate pyscf-env

# Verify the setup
npm run verify-env
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

## Individual Commands

### Environment Management
```bash
# Automated environment setup (conda + dependencies)
npm run setup-env

# Verify environment health and dependencies
npm run verify-env
```

### Build Commands
```bash
# Clean build directory
npm run clean

# Generate TypeScript types and Python models from OpenAPI spec
npm run codegen

# Build frontend with webpack in development mode
npm run dev:webpack

# Start Electron (requires dist files to exist)
npm run dev:electron

# Package conda environment for distribution
npm run build:conda-pack

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

# Code formatting (using Prettier)
npm run format      # Format all source code
npm run format:check # Check if code is properly formatted

Development Workflow
The npm run dev script is the primary command for development. It automatically:

Generates Code: Runs npm run codegen to generate TypeScript types and Python Pydantic models from src/api-spec/openapi.yaml. This ensures the frontend and backend are always in sync with the API specification.

Cleans: Cleans the dist/ directory.

Builds Frontend: Builds the frontend code (main process, preload script, and renderer) using Webpack in watch mode.

Starts Backend: Starts the Python Flask server as a subprocess from within the Electron main process.

Starts Electron: Starts the Electron application using electronmon, which watches the dist/ directory for changes and automatically restarts the app.

The Electron main process (src/main.ts) launches the Python Flask server using an enhanced environment detection system. The `detectPythonEnvironmentPath` function implements a hierarchical detection strategy:

**Environment Detection Priority:**
1. **Bundled conda environment** (packaged apps): Uses conda-pack packaged environment at `process.resourcesPath/conda_env/bin/python`
2. **Development conda environment**: Detects pyscf-env using multiple strategies (environment variables, conda commands, fallback paths)
3. **PyInstaller executable** (fallback): Uses the standalone executable for packaged apps

**Health Check System:** The main process includes a robust health check mechanism that continuously pings the `/health` endpoint with detailed diagnostic information. If the server fails to start, it provides context-specific error messages with troubleshooting guidance.

**Dynamic Port Management:** The Flask server automatically finds a free port and reports it to the main process via stdout, ensuring no port conflicts.

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

Code Quality & Formatting
The project uses Prettier for consistent code formatting across all TypeScript, JavaScript, and JSON files:

Configuration: Code formatting rules are defined in .prettierrc with sensible defaults for the project's style.

Automatic Formatting: Use npm run format to automatically format all source files, or npm run format:check to verify formatting compliance.

CI Integration: GitHub Actions CI enforces code formatting - PRs with formatting issues will fail the build, ensuring consistent code style across all contributions.

Exclusions: Auto-generated files (like generated-api.ts and generated_models.py) are excluded via .prettierignore to prevent conflicts with code generation.

Build System
Frontend: Uses Webpack with ts-loader to compile TypeScript and React code into three separate bundles: main.js, preload.js, and the renderer's app.js.

Backend: Implements a dual-packaging strategy for maximum compatibility:
1. **conda-pack Integration**: Packages the complete conda environment (including PySCF, RDKit, and all dependencies) into a portable format using `npm run build:conda-pack`. This creates a `conda_env/` directory that can be bundled with the application.
2. **PyInstaller Executable**: Bundles the Flask application into a standalone executable as a fallback option.

Packaging: electron-builder packages the Electron app with both Python environments (included as extraResources):
- `conda_env/` - Complete portable conda environment (primary)
- `python_dist/` - PyInstaller executable (fallback)
- Distributable formats: DMG for macOS, NSIS for Windows, AppImage for Linux

**Build Process Flow:**
1. `npm run build:conda-pack` - Packages conda environment
2. `npm run build:python` - Creates PyInstaller executable  
3. `npm run build:webpack` - Builds frontend
4. `electron-builder` - Creates final distributable with both environments

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
├── .github
    ├── environment.yml
    └── workflows
    │   ├── ci.yml
    │   └── release.yml
├── .gitignore
├── .prettierignore
├── .prettierrc
├── CLAUDE.md
├── LICENSE
├── README.md
├── package-lock.json
├── package.json
├── scripts
    ├── setup-environment.sh
    └── verify-environment.py
├── src
    ├── api-spec
    │   └── openapi.yaml
    ├── assets
    │   ├── fonts
    │   │   └── ADLaMDisplay-Regular.ttf
    │   └── icon
    │   │   └── mac
    │   │       └── Pyscf_front.icns
    ├── main.ts
    ├── preload.ts
    ├── python
    │   ├── SMILES
    │   │   ├── __init__.py
    │   │   └── smiles_converter.py
    │   ├── app.py
    │   ├── data
    │   │   ├── __init__.py
    │   │   └── solvent_properties.py
    │   ├── generated_models.py
    │   ├── pubchem
    │   │   ├── __init__.py
    │   │   ├── client.py
    │   │   └── parser.py
    │   ├── quantum_calc
    │   │   ├── __init__.py
    │   │   ├── base_calculator.py
    │   │   ├── ccsd_calculator.py
    │   │   ├── dft_calculator.py
    │   │   ├── exceptions.py
    │   │   ├── file_manager.py
    │   │   ├── file_watcher.py
    │   │   ├── hf_calculator.py
    │   │   ├── mp2_calculator.py
    │   │   ├── orbital_generator.py
    │   │   ├── process_manager.py
    │   │   ├── solvent_effects.py
    │   │   └── tddft_calculator.py
    │   └── tests
    │   │   ├── __init__.py
    │   │   ├── test_concurrent_websockets.py
    │   │   ├── test_error_scenarios.py
    │   │   ├── test_file_watcher.py
    │   │   ├── test_flask_api.py
    │   │   ├── test_pubchem.py
    │   │   ├── test_quantum_calc.py
    │   │   └── test_websocket_integration.py
    ├── types
    │   ├── 3dmol.d.ts
    │   └── electron.d.ts
    └── web
    │   ├── App.css
    │   ├── App.tsx
    │   ├── apiClient.ts
    │   ├── components
    │       ├── DropdownMenu.tsx
    │       ├── Header.tsx
    │       ├── MolecularOrbitalEnergyDiagram.tsx
    │       ├── MolecularOrbitalViewer.tsx
    │       ├── MoleculeViewer.tsx
    │       ├── MoleculeViewerSection.tsx
    │       ├── Sidebar.tsx
    │       ├── StyleControls.tsx
    │       └── XYZInput.tsx
    │   ├── data
    │       └── atomicRadii.ts
    │   ├── hooks
    │       ├── index.ts
    │       ├── useActiveCalculation.ts
    │       ├── useCalculationOperations.ts
    │       ├── useCalculationQueries.ts
    │       ├── useCalculationSubscription.ts
    │       ├── usePageNavigation.ts
    │       └── useSidebarState.ts
    │   ├── index.html
    │   ├── index.tsx
    │   ├── pages
    │       ├── CalculationResultsPage.tsx
    │       ├── CalculationSettingsPage.tsx
    │       └── DrawMoleculePage.tsx
    │   ├── store
    │       └── calculationStore.ts
    │   ├── types
    │       ├── api-types.ts
    │       └── generated-api.ts
    │   └── utils
    │       └── xyzParser.ts
├── tsconfig.json
└── webpack.config.ts

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

## Environment Setup Issues

**Automated Diagnosis:** Always start with environment verification:
```bash
npm run verify-env
```
This command provides comprehensive diagnostic information and specific troubleshooting recommendations.

**Environment Detection Failures:** If the application reports "Python environment not found", check the detection hierarchy:
1. Verify conda is installed and accessible: `conda --version`
2. Check if pyscf-env exists: `conda env list`
3. Try automated setup: `npm run setup-env`
4. Set custom path if needed: `export CONDA_ENV_PATH=/path/to/your/pyscf-env`

**conda-pack Build Issues:** If `npm run build:conda-pack` fails:
- Ensure conda-pack is installed: `conda install conda-pack`
- Verify environment is active: `conda activate pyscf-env`
- Check for disk space and permissions

## Development Server Issues

**Backend Server Fails to Start:** The enhanced error system now provides detailed diagnostic information:
- **Development mode**: Check conda environment and dependencies
- **Packaged mode**: Verify bundled environment integrity
- The Flask server automatically finds free ports, but firewall settings may interfere

**Environment Health Check Failures:** The application performs extensive health checks:
- Python version compatibility (3.10+)
- Required packages (PySCF, RDKit, Flask, etc.)
- Functional tests for critical dependencies
- Port availability and network configuration

## Build and Packaging Issues

**Build Failures:** The dual-packaging system requires:
1. Working conda environment for conda-pack
2. All dependencies installed for PyInstaller
3. Check `build/pyinstaller/` directory for detailed logs

**PyInstaller Specific Issues:** 
- Sensitive to Python environment changes
- Requires all dependencies to be importable
- Check for missing hidden imports in `pyscf_front_api.spec`

**Packaging Verification:** After packaging, verify both environments:
- Test bundled conda environment path
- Confirm PyInstaller executable functionality
- Check extraResources in electron-builder output

## Runtime Issues

**Process Pool Problems:** Enhanced process management with better error reporting:
- Check `/api/quantum/status` for pool health
- Process manager handles cleanup automatically  
- Restart Flask server if pool becomes unresponsive

**SSL/HTTPS Issues:** PubChem integration uses standard HTTPS with certificate verification enabled.

## Development Tools

**Environment Validation Tools:**
- `npm run verify-env` - Comprehensive environment testing
- `npm run setup-env` - Automated environment creation and repair
- Console output provides detailed diagnostic information
- Health check system reports specific failure points

**Debugging Tips:**
- Enable verbose logging in development mode
- Check both stdout and stderr from Python processes
- Monitor resource usage during conda-pack operations
- Verify file permissions for packaged environments