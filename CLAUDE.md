This file provides guidance to Claude when working with code in this repository. Your partner is Japanese, so please always report in Japanese.

Project Overview
This is PySCF_front, an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates, retrieving molecular structures from PubChem and SMILES strings, and visualizing 3D molecular structures using 3Dmol.js. The backend is a Python Flask server that handles all chemical computations and data management using libraries like PySCF and RDKit.

The application supports various quantum chemistry calculation methods, including DFT, Hartree-Fock (HF), MP2, CCSD, TDDFT, CASCI, and CASSCF. It also features geometry optimization, vibrational frequency analysis with IR spectrum visualization, and advanced analysis like Molecular Orbitals (MO) and Natural Transition Orbitals (NTO) for TDDFT. A recent addition is an AI-powered molecular agent that can perform tasks based on natural language prompts. It dynamically loads supported parameters (basis sets, functionals, etc.) from the backend.

A key feature of this project is its API-first development approach, using an OpenAPI specification as the single source of truth for the API contract between the frontend and backend. The application uses WebSockets for real-time status updates of running calculations, providing a more efficient and responsive user experience than a polling-based system.

Development Philosophy
This is a development-stage application. Backward compatibility is not a concern, and breaking changes should be made freely in favor of better design and simpler code. When refactoring or improving the codebase:

Prioritize simplicity over compatibility - Remove deprecated patterns and complex fallback logic.

Make breaking changes confidently - Don't hesitate to change APIs, data structures, or file formats.

Clean up legacy code - Remove old implementations when better alternatives are available.

Focus on the best solution - Don't compromise design quality for compatibility with older versions.

This approach allows for rapid iteration and prevents technical debt accumulation during the development phase.

AI Development Guidelines
When working with external libraries, frameworks, or implementing new features, Claude should ALWAYS:

**Verify Current API Documentation First**

Before writing any code that uses external libraries (LangChain, LangGraph, Flask, React, etc.), ALWAYS use Web Search and Context7 tools to verify the current API:
```
1. Use WebSearch to find recent updates and breaking changes (e.g., "langgraph 2025 API changes")
2. Use Context7 to fetch the latest official documentation for the specific library
3. Check for deprecated parameters or methods in the current version
```

**Mandatory Checks**
- API parameter names and signatures (they change frequently!)
- Deprecated methods or parameters
- New recommended patterns or best practices
- Breaking changes in recent versions

**Example Workflow**

When implementing a feature with LangGraph:
```
1. WebSearch: "langgraph stream messages 2025" to find recent changes
2. Context7: Fetch /langchain-ai/langgraph docs for "create_react_agent parameters"
3. Verify the exact parameter names (e.g., prompt vs messages_modifier vs stateModifier)
4. Write code using the verified, current API
```

**Why This Matters**

Libraries like LangChain and LangGraph update frequently, and parameters/methods get renamed or deprecated. Using outdated APIs leads to runtime errors that could have been avoided. Always verify before coding.

Development Commands
Initial Setup
Conda Environment (Required)
This project requires a conda environment for development. The application features automated environment setup and verification tools, and uses a unified server configuration system for consistent behavior across development and production environments.

Quick Setup (Recommended)

Bash

# Install Node.js dependencies
npm install

# Automated environment setup (handles all conda setup)
npm run setup-env

# Verify environment health
npm run verify-env

# Verify build tools
npm run verify-build-env

# Check server configuration
npm run debug:config
Manual Setup

Bash

# Install Node.js dependencies
npm install

# Install Miniforge (if not already installed)
# Example for macOS ARM:
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3

# Create conda environment from environment.yml (includes all dependencies)
source $HOME/miniforge3/etc/profile.d/conda.sh
conda env create -f .github/environment.yml

# Activate the environment
conda activate pyscf-env

# Verify the setup
npm run verify-env
Note: The conda environment setup is mandatory. The application will show an error dialog if the conda environment is not properly configured.

Development Commands
Bash

# Activate conda environment
conda activate pyscf-env

# Development mode (generates code, builds, and runs Electron with hot reload + Python backend)
npm run dev

# Production build (generates code, builds frontend and creates Python executable)
npm run build

# Package application for distribution (includes production build)
npm run package
Individual Commands
Environment Management

Bash

# Automated environment setup (conda + dependencies)
npm run setup-env

# Verify environment health and dependencies
npm run verify-env

# Verify build tools (conda-pack, Gunicorn)
npm run verify-build-env

# Debug server configuration
npm run debug:config
Build Commands

Bash

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

# Validate build completeness (after build)
npm run validate-build
Code Formatting

Bash

# Format all source code (using Prettier)
npm run format

# Check if code is properly formatted
npm run format:check
Testing and Validation Commands

Bash

# Complete build test (frontend + backend)
npm run test:build

# Test Python imports and dependencies
npm run test:python-build

# Test Gunicorn server locally with unified configuration
npm run test:gunicorn-local

# Full packaging test (build + package)
npm run test:run-packaged

# --- Manual Python Backend Testing (in a separate terminal) ---
cd src/python

# Activate conda environment
conda activate pyscf-env

# Start Flask API server directly
python app.py

# Run Python backend tests
pytest tests/
Development Workflow
The npm run dev script is the primary command for development. It automatically:

Generates Code: Runs npm run codegen to generate TypeScript types and Python Pydantic models from src/api-spec/openapi.yaml. This ensures the frontend and backend are always in sync with the API specification.

Cleans: Cleans the dist/ directory.

Builds Frontend: Builds the frontend code (main process, preload script, and renderer) using Webpack in watch mode.

Starts Backend: Starts the Python Flask server as a subprocess from within the Electron main process using the unified Gunicorn-based execution system.

Starts Electron: Starts the Electron application using electronmon, which watches the dist/ directory for changes and automatically restarts the app.

Unified Server Configuration System
The application uses a configuration-driven approach with config/server-config.json to ensure consistent behavior across development and production environments.

Server settings: Host, port, and auto-detection.

Gunicorn configuration: Workers, threads, and timeout settings.

SocketIO settings: CORS, async mode, and timeouts.

Logging configuration.

Unified Execution Environment: Both development and production environments now use Gunicorn as the WSGI server, eliminating "works in dev but not in production" problems.

The Electron main process (src/main.ts) launches the Python Flask server using a simplified environment detection strategy:

Bundled conda environment (packaged apps): Uses the conda-pack packaged environment at process.resourcesPath/conda_env/bin/python.

Development conda environment: Detects the pyscf-env using environment variables or conda commands.

A robust health check system pings the /health endpoint. If the server fails to start, it provides context-specific error messages.

Architecture Overview
API-First Development with OpenAPI
The single source of truth for the API is src/api-spec/openapi.yaml. The npm run codegen command uses this file to generate:

Python Pydantic Models (src/python/generated_models.py) for type-safe request/response handling in the Flask backend.

TypeScript Type Definitions (src/web/types/generated-api.ts) to keep the frontend API client synchronized with the backend.

Electron Structure
Main Process (src/main.ts): Creates the BrowserWindow, manages the Python Flask subprocess, and handles application lifecycle events.

Preload Script (src/preload.ts): Securely exposes specific Electron APIs to the renderer process.

Renderer Process (src/web/): The React (TypeScript) application that provides the user interface.

Python Backend (src/python/)
A Flask API server with REST endpoints and a WebSocket interface for:

PubChem & SMILES Integration: Searching and converting molecular structures.

Quantum Chemistry Calculations: Running various calculations (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF) via PySCF. This now includes geometry optimization and vibrational frequency analysis. Calculations are executed in parallel using a ProcessPoolExecutor.

Molecular Agent: An AI agent that interprets natural language prompts to perform molecular calculations and analysis using a set of defined tools.

Molecular Orbital Analysis: Generating data for visualizing molecular orbitals, including CUBE files and energy level diagrams.

IR Spectrum Generation: Creating theoretical IR spectra from frequency analysis data.

Dynamic Parameter Loading: Providing lists of supported basis sets, functionals, and solvents via the /api/quantum/supported-parameters endpoint.

Real-time Status Updates: A WebSocket endpoint pushes status updates to the frontend.

File Management: Listing, renaming, and deleting calculation data.

Health Check: An endpoint (/health) for startup coordination.

Core Components & State Management
The application uses a modern, hook-based state management architecture with a clear separation of concerns, moving complex logic out of App.tsx and into reusable hooks.

State Management Architecture:

TanStack Query: Manages all server-side state, including API data fetching, caching, and synchronization.

Zustand: Manages global UI state.

calculationStore.ts: Manages activeCalculationId and stagedCalculation (for new calculation workflows).

notificationStore.ts: Manages global toast notifications for user feedback.

uiStore.ts: Manages UI state like sidebar visibility, current page, and modals.

Key Hooks (src/web/hooks/)
This architecture relies on a collection of custom hooks to encapsulate logic:

useAppState: A central hook that combines UI and calculation state from Zustand stores for simplified access in App.tsx.

useCalculationData: Derives the currently active calculation state. It intelligently selects between a staged (new) calculation, detailed data fetched from the server, or a fallback, providing a unified activeCalculation object to the UI.

useCalculationActions: Encapsulates all mutation actions like starting, renaming, and deleting calculations. It handles API calls and onSuccess/onError logic, including notifications.

useUnifiedWebSocket: A dedicated hook that establishes a WebSocket connection for both global updates and the active calculation, updating the TanStack Query cache in real-time.

useCalculationQueries: Contains all TanStack Query definitions (useQuery, useMutation) for interacting with the backend API.

State Flow Example (Starting a Calculation)

User clicks the "Start Calculation" button.

The handleStartCalculation function from the useCalculationActions hook is called. It uses the useStartCalculation mutation.

The Flask backend (POST /api/quantum/calculate) receives the request, creates a directory, saves initial parameters, sets the status to running, and submits the job to a ProcessPoolExecutor. It returns a 202 Accepted response with the new CalculationInstance data.

The useStartCalculation mutation's onSuccess callback invalidates the calculations query cache, triggering a refetch.

The useCalculationStore clears any staged data, and sets the new activeCalculationId.

The useUnifiedWebSocket hook detects the running status of the new active calculation and opens a WebSocket connection.

The Flask backend's WebSocket handler monitors status.json. As the worker process updates the file (e.g., to completed or error), the handler pushes the complete, updated CalculationInstance data to the client.

The onUpdate callback in useUnifiedWebSocket receives the new data and directly updates the TanStack Query cache using queryClient.setQueryData(), triggering re-renders in all components using that data.

File Structure
.
├── .github/
│   ├── environment.yml
│   └── workflows/
│       ├── ci.yml
│       └── release.yml
├── config/
│   └── server-config.json
├── scripts/
│   ├── setup-environment.sh
│   ├── test-python-standalone.js
│   ├── validate-build-completeness.py
│   └── verify-environment.py
├── src/
│   ├── api-spec/
│   │   └── openapi.yaml
│   ├── assets/
│   │   ├── fonts/
│   │   └── icon/
│   ├── main.ts
│   ├── preload.ts
│   ├── python/
│   │   ├── SMILES/
│   │   │   ├── __init__.py
│   │   │   └── smiles_converter.py
│   │   ├── agent/
│   │   │   ├── __init__.py
│   │   │   ├── molecular_agent.py
│   │   │   ├── prompts/
│   │   │   │   └── system_prompt.txt
│   │   │   └── tools.py
│   │   ├── api/
│   │   │   ├── __init__.py
│   │   │   ├── agent.py
│   │   │   ├── health.py
│   │   │   ├── pubchem.py
│   │   │   ├── quantum.py
│   │   │   └── ...
│   │   ├── quantum_calc/
│   │   │   └── ...
│   │   ├── tests/
│   │   │   └── ...
│   │   ├── websocket/
│   │   │   └── ...
│   │   ├── app.py
│   │   └── generated_models.py
│   ├── types/
│   │   └── ...
│   └── web/
│       ├── components/
│       ├── hooks/
│       ├── pages/
│       ├── store/
│       ├── types/
│       ├── utils/
│       ├── App.tsx
│       ├── apiClient.ts
│       ├── index.html
│       └── index.tsx
├── package.json
├── tsconfig.json
└── webpack.config.ts
Key API Endpoints
GET /health: Health check endpoint.

POST /api/pubchem/search: Search PubChem.

POST /api/smiles/convert: Convert a SMILES string to XYZ.

POST /api/agent/invoke: Invoke the molecular agent with a natural language prompt.

GET /api/quantum/supported-parameters: Get lists of supported calculation methods, basis sets, functionals, etc.

POST /api/quantum/calculate: Asynchronously starts a quantum chemistry calculation.

GET /api/quantum/calculations: Lists all saved calculations.

GET /api/quantum/calculations/<id>: Gets detailed results for a specific calculation.

PATCH /api/quantum/calculations/<id>: Updates a calculation's metadata (e.g., renames it).

DELETE /api/quantum/calculations/<id>: Deletes a calculation.

GET /api/quantum/calculations/<id>/orbitals: Get molecular orbital information.

GET /api/quantum/calculations/<id>/orbitals/{orbitalIndex}/cube: Generate and retrieve a CUBE file for a specific orbital.

GET /api/quantum/calculations/<id>/ir-spectrum: Generate an IR spectrum for a calculation.

GET /api/settings: Get application settings.

PUT /api/settings: Update application settings.

GET /api/system/resource-status: Get system resource status.

WS /ws/calculations/<id>: WebSocket endpoint for real-time status updates.

Parallel Processing Architecture
The application uses a ProcessPoolExecutor-based system for quantum chemistry calculations to achieve true multiprocessing parallelism.

Process Pool Management: The CalculationProcessManager class manages a pool of worker processes, allowing multiple calculations to run simultaneously.

Worker Process Isolation: Each calculation runs in a separate Python process, eliminating the Global Interpreter Lock (GIL) limitation.

Resource Management: The process pool starts on demand and shuts down cleanly when the Flask server terminates.

Spin Multiplicity Notes
PySCF uses the spin attribute to specify the number of unpaired electrons (2S), which differs from the traditional quantum chemistry notation of 2S+1.

Singlet state (closed-shell): spin=0

Doublet state (1 unpaired electron): spin=1

Triplet state (2 unpaired electrons): spin=2

Troubleshooting
Environment Setup Issues
Automated Diagnosis: Always start with environment verification: npm run verify-env. This command provides comprehensive diagnostic information.

Environment Detection Failures: If detection fails, use npm run setup-env or set the CONDA_ENV_PATH environment variable.

Build Tool Issues: If npm run verify-build-env fails, ensure conda-pack and gunicorn are installed in the active pyscf-env environment.

Development Server Issues
The unified Gunicorn-based server reduces environment-specific issues. Most behavior is controlled by config/server-config.json. Use npm run debug:config to verify server settings. The server automatically finds free ports.

Build and Packaging Issues
The build system includes pre- and post-build verification steps (verify-env, verify-build-env, validate-build). A failure in these steps indicates an issue with the conda environment or the server configuration.