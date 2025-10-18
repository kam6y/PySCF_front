This file provides guidance to Claude when working with code in this repository. Your partner is Japanese, so please always report in Japanese.

## Project Overview
This is PySCF_front, an Electron-based desktop application for molecular visualization and quantum chemistry calculations. The app provides a React-based UI for inputting XYZ molecular coordinates, retrieving molecular structures from PubChem and SMILES strings, and visualizing 3D molecular structures using 3Dmol.js. The backend is a Python Flask server that handles all chemical computations and data management using libraries like PySCF and RDKit.

The application supports various quantum chemistry calculation methods, including DFT, Hartree-Fock (HF), MP2, CCSD, TDDFT, CASCI, and CASSCF. It also features geometry optimization, vibrational frequency analysis with IR spectrum visualization, and advanced analysis like Molecular Orbitals (MO) and Natural Transition Orbitals (NTO) for TDDFT. A recent addition is an AI-powered molecular agent that can perform tasks based on natural language prompts. It dynamically loads supported parameters (basis sets, functionals, etc.) from the backend.

A key feature of this project is its API-first development approach, using an OpenAPI specification as the single source of truth for the API contract between the frontend and backend. The application uses WebSockets for real-time status updates of running calculations, providing a more efficient and responsive user experience than a polling-based system.

## Development Philosophy
This is a development-stage application. Backward compatibility is not a concern, and breaking changes should be made freely in favor of better design and simpler code. When refactoring or improving the codebase:

* Prioritize simplicity over compatibility - Remove deprecated patterns and complex fallback logic.
* Make breaking changes confidently - Don't hesitate to change APIs, data structures, or file formats.
* Clean up legacy code - Remove old implementations when better alternatives are available.
* Focus on the best solution - Don't compromise design quality for compatibility with older versions.

This approach allows for rapid iteration and prevents technical debt accumulation during the development phase.

## AI Development Guidelines
When working with external libraries, frameworks, or implementing new features, Claude should ALWAYS:

### Verify Current API Documentation First

Before writing any code that uses external libraries (LangChain, LangGraph, Flask, React, etc.), ALWAYS use Web Search and Context7 tools to verify the current API:

    1. Use WebSearch to find recent updates and breaking changes (e.g., "langgraph 2025 API changes")
    2. Use Context7 to fetch the latest official documentation for the specific library
    3. Check for deprecated parameters or methods in the current version

### Mandatory Checks
* API parameter names and signatures (they change frequently!)
* Deprecated methods or parameters
* New recommended patterns or best practices
* Breaking changes in recent versions

### Example Workflow

When implementing a feature with LangGraph:

    1. WebSearch: "langgraph stream messages 2025" to find recent changes
    2. Context7: Fetch /langchain-ai/langgraph docs for "create_react_agent parameters"
    3. Verify the exact parameter names (e.g., prompt vs messages_modifier vs stateModifier)
    4. Write code using the verified, current API

### Why This Matters

Libraries like LangChain and LangGraph update frequently, and parameters/methods get renamed or deprecated. Using outdated APIs leads to runtime errors that could have been avoided. Always verify before coding.

## Development Commands
### Initial Setup
#### Conda Environment (Required)
This project requires a conda environment for development. The application features automated environment setup and verification tools, and uses a unified server configuration system for consistent behavior across development and production environments.

#### Quick Setup (Recommended)

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

#### Manual Setup

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

### Development Commands

    # Activate conda environment
    conda activate pyscf-env
    
    # Development mode (generates code, builds, and runs Electron with hot reload + Python backend)
    npm run dev
    
    # Production build (generates code, builds frontend and creates Python executable)
    npm run build
    
    # Package application for distribution (includes production build)
    npm run package

### Individual Commands
#### Environment Management

    # Automated environment setup (conda + dependencies)
    npm run setup-env
    
    # Verify environment health and dependencies
    npm run verify-env
    
    # Verify build tools (conda-pack, Gunicorn)
    npm run verify-build-env
    
    # Debug server configuration
    npm run debug:config

#### Build Commands

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

#### Code Formatting

    # Format all source code (using Prettier)
    npm run format
    
    # Check if code is properly formatted
    npm run format:check

#### Testing and Validation Commands

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

## Development Workflow
The `npm run dev` script is the primary command for development. It automatically:

* **Generates Code**: Runs `npm run codegen` to generate TypeScript types and Python Pydantic models from `src/api-spec/openapi.yaml`. This ensures the frontend and backend are always in sync with the API specification.
* **Cleans**: Cleans the `dist/` directory.
* **Builds Frontend**: Builds the frontend code (main process, preload script, and renderer) using Webpack in watch mode.
* **Starts Backend**: Starts the Python Flask server as a subprocess from within the Electron main process using the unified Gunicorn-based execution system.
* **Starts Electron**: Starts the Electron application using `electronmon`, which watches the `dist/` directory for changes and automatically restarts the app.

## Unified Server Configuration System
The application uses a configuration-driven approach with `config/server-config.json` to ensure consistent behavior across development and production environments.

* Server settings: Host, port, and auto-detection.
* Gunicorn configuration: Workers, threads, and timeout settings.
* SocketIO settings: CORS, async mode, and timeouts.
* Logging configuration.

**Unified Execution Environment**: Both development and production environments now use Gunicorn as the WSGI server, eliminating "works in dev but not in production" problems.

The Electron main process (`src/main.ts`) launches the Python Flask server using a simplified environment detection strategy:

1.  **Bundled conda environment** (packaged apps): Uses the conda-pack packaged environment at `process.resourcesPath/conda_env/bin/python`.
2.  **Development conda environment**: Detects the `pyscf-env` using environment variables or conda commands.

A robust health check system pings the `/health` endpoint. If the server fails to start, it provides context-specific error messages.

## Architecture Overview
### API-First Development with OpenAPI
The single source of truth for the API is `src/api-spec/openapi.yaml`. The `npm run codegen` command uses this file to generate:

* **Python Pydantic Models** (`src/python/generated_models.py`) for type-safe request/response handling in the Flask backend.
* **TypeScript Type Definitions** (`src/web/types/generated-api.ts`) to keep the frontend API client synchronized with the backend.

### Multi-Agent AI Architecture
The application features a sophisticated multi-agent AI system powered by LangGraph using a **Supervisor pattern**. This replaces the previous simple router-based architecture.

* **Supervisor (Coordinator)**: Analyzes user intent and delegates tasks to specialized worker agents. It coordinates complex, multi-step workflows.
* **Quantum Calculator**: Handles all quantum chemistry calculations, molecular property analysis, orbital visualization, and general chemistry-related queries. Uses a comprehensive set of tools to interact with the PySCF calculation backend.
* **Literature Surveyor**: Specializes in academic literature search. It can iteratively search multiple sources (arXiv, Tavily, PubMed, etc.), analyze findings, and synthesize comprehensive reports with citations.
* **Science Analyst**: Handles scientific report generation, data interpretation, and visualization tasks based on calculation results or research findings.

The system uses LangGraph's stateful graph architecture, enabling intelligent task distribution and complex, multi-turn interactions.

### Electron Structure
* **Main Process** (`src/main.ts`): Creates the `BrowserWindow`, manages the Python Flask subprocess, and handles application lifecycle events.
* **Preload Script** (`src/preload.ts`): Securely exposes specific Electron APIs to the renderer process.
* **Renderer Process** (`src/web/`): The React (TypeScript) application that provides the user interface.

### Python Backend (`src/python/`)
A Flask API server with REST endpoints and a WebSocket interface for:

* **PubChem & SMILES Integration**: Searching and converting molecular structures.
* **Quantum Chemistry Calculations**: Running various calculations (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF) via PySCF. This now includes geometry optimization and vibrational frequency analysis. Calculations are executed in parallel using a `ProcessPoolExecutor`.
* **Multi-Agent AI System**: A LangGraph-based system with intelligent routing between specialized agents:
    * **Supervisor**: Coordinates tasks between specialized agents.
    * **Quantum Calculator**: Interprets natural language prompts to perform molecular calculations and analysis using a comprehensive toolset.
    * **Literature Surveyor**: Searches academic databases (arXiv, Tavily, PubMed, etc.) for papers, providing summaries and formatted citations.
    * **Science Analyst**: Generates reports and interprets data from calculations and research.
* **Molecular Orbital Analysis**: Generating data for visualizing molecular orbitals, including CUBE files and energy level diagrams.
* **IR Spectrum Generation**: Creating theoretical IR spectra from frequency analysis data.
* **Dynamic Parameter Loading**: Providing lists of supported basis sets, functionals, and solvents via the `/api/quantum/supported-parameters` endpoint.
* **Real-time Status Updates**: A WebSocket endpoint pushes status updates to the frontend.
* **File Management**: Listing, renaming, and deleting calculation data.
* **Health Check**: An endpoint (`/health`) for startup coordination.

## Core Components & State Management
The application uses a modern, hook-based state management architecture with a clear separation of concerns, moving complex logic out of `App.tsx` and into reusable hooks.

### State Management Architecture:

* **TanStack Query**: Manages all server-side state, including API data fetching, caching, and synchronization.
* **Zustand**: Manages global UI state.
    * `agentStore.ts`: Manages AI agent state, such as which agent is active.
    * `calculationStore.ts`: Manages `activeCalculationId` and `stagedCalculation` (for new calculation workflows).
    * `chatHistoryStore.ts`: Manages the state for the chat history sidebar and sessions.
    * `notificationStore.ts`: Manages global toast notifications for user feedback.
    * `uiStore.ts`: Manages UI state like sidebar visibility (instances vs. chats), current page, and modals.

### Key Hooks (`src/web/hooks/`)
This architecture relies on a collection of custom hooks to encapsulate logic:

* **useActiveCalculation**: (Replaces `useCalculationData`) Derives the currently active calculation state. It intelligently selects between a staged (new) calculation, detailed data fetched from the server, or a fallback, providing a unified `activeCalculation` object to the UI.
* **useActiveCalculationId**: Manages the currently selected `activeCalculationId`.
* **useAppSettings**: Hook for fetching and mutating application settings.
* **useAppState**: A central hook that combines UI and calculation state from Zustand stores for simplified access in `App.tsx`.
* **useCalculationActions**: Encapsulates mutation actions (start, rename, delete) for calculations.
* **useCalculationOperations**: (New) Encapsulates complex calculation-related logic and state derivations.
* **useCalculationQueries**: Contains all TanStack Query definitions (`useQuery`, `useMutation`) for the quantum calculation API.
* **useChatHistoryQueries**: (New) Contains all TanStack Query definitions for the chat history API.
* **useUnifiedWebSocket**: A dedicated hook that establishes a WebSocket connection for both global updates and the active calculation, updating the TanStack Query cache in real-time.

### State Flow Example (Starting a Calculation)

1.  User clicks the "Start Calculation" button.
2.  The `handleStartCalculation` function from the `useCalculationActions` hook is called. It uses the `useStartCalculation` mutation.
3.  The Flask backend (POST `/api/quantum/calculate`) receives the request, creates a directory, saves initial parameters, sets the status to `running`, and submits the job to a `ProcessPoolExecutor`. It returns a `202 Accepted` response with the new `CalculationInstance` data.
4.  The `useStartCalculation` mutation's `onSuccess` callback invalidates the `calculations` query cache, triggering a refetch.
5.  The `useCalculationStore` clears any staged data, and sets the new `activeCalculationId`.
6.  The `useUnifiedWebSocket` hook detects the `running` status of the new active calculation and opens a WebSocket connection.
7.  The Flask backend's WebSocket handler monitors `status.json`. As the worker process updates the file (e.g., to `completed` or `error`), the handler pushes the complete, updated `CalculationInstance` data to the client.
8.  The `onUpdate` callback in `useUnifiedWebSocket` receives the new data and directly updates the TanStack Query cache using `queryClient.setQueryData()`, triggering re-renders in all components using that data.

## File Structure
    .
    ├── .github/
    │   ├── environment.yml
    │   └── workflows/
    │       ├── ci.yml
    │       └── release.yml
    ├── .gitignore
    ├── .prettierignore
    ├── .prettierrc
    ├── CLAUDE.md
    ├── LICENSE
    ├── PySCF_front_view.png
    ├── README.md
    ├── TESTING_IMPLEMENTATION_SUMMARY.md
    ├── config/
    │   └── server-config.json
    ├── data/
    │   └── .gitkeep
    ├── package-lock.json
    ├── package.json
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
    │   │   │   └── ADLaMDisplay-Regular.ttf
    │   │   └── icon/
    │   │       └── mac/
    │   │           └── Pyscf_front.icns
    │   ├── main.ts
    │   ├── preload.ts
    │   ├── python/
    │   │   ├── SMILES/
    │   │   │   ├── __init__.py
    │   │   │   └── smiles_converter.py
    │   │   ├── __init__.py
    │   │   ├── agent/
    │   │   │   ├── __init__.py
    │   │   │   ├── graph.py
    │   │   │   ├── literature_surveyor/
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── literature_surveyor_agent.py
    │   │   │   │   ├── prompts/
    │   │   │   │   │   ├── analyze_prompt.txt
    │   │   │   │   │   └── synthesize_prompt.txt
    │   │   │   │   └── tools.py
    │   │   │   ├── quantum_calculator/
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── prompts/
    │   │   │   │   │   └── system_prompt.txt
    │   │   │   │   ├── quantum_calculator_agent.py
    │   │   │   │   └── tools.py
    │   │   │   ├── science_analyst/
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── prompts/
    │   │   │   │   │   └── system_prompt.txt
    │   │   │   │   ├── science_analyst_agent.py
    │   │   │   │   └── tools.py
    │   │   │   ├── supervisor/
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── prompts/
    │   │   │   │   │   └── system_prompt.txt
    │   │   │   │   └── supervisor.py
    │   │   │   └── utils.py
    │   │   ├── api/
    │   │   │   ├── __init__.py
    │   │   │   ├── agent.py
    │   │   │   ├── chat_history.py
    │   │   │   ├── health.py
    │   │   │   ├── pubchem.py
    │   │   │   ├── quantum.py
    │   │   │   ├── settings.py
    │   │   │   ├── smiles.py
    │   │   │   └── system.py
    │   │   ├── app.py
    │   │   ├── config.py
    │   │   ├── data/
    │   │   │   ├── __init__.py
    │   │   │   ├── scale_factors.py
    │   │   │   └── solvent_properties.py
    │   │   ├── database/
    │   │   │   ├── __init__.py
    │   │   │   └── chat_history.py
    │   │   ├── generated_models.py
    │   │   ├── pubchem/
    │   │   │   ├── __init__.py
    │   │   │   ├── client.py
    │   │   │   └── parser.py
    │   │   ├── pytest.ini
    │   │   ├── quantum_calc/
    │   │   │   ├── __init__.py
    │   │   │   ├── base_calculator.py
    │   │   │   ├── casci_calculator.py
    │   │   │   ├── casscf_calculator.py
    │   │   │   ├── ccsd_calculator.py
    │   │   │   ├── config_manager.py
    │   │   │   ├── dft_calculator.py
    │   │   │   ├── exceptions.py
    │   │   │   ├── file_manager.py
    │   │   │   ├── file_watcher.py
    │   │   │   ├── hf_calculator.py
    │   │   │   ├── ir_spectrum.py
    │   │   │   ├── mp2_calculator.py
    │   │   │   ├── orbital_generator.py
    │   │   │   ├── process_manager.py
    │   │   │   ├── resource_manager.py
    │   │   │   ├── settings_manager.py
    │   │   │   ├── solvent_effects.py
    │   │   │   ├── supported_parameters.py
    │   │   │   └── tddft_calculator.py
    │   │   ├── services/
    │   │   │   ├── __init__.py
    │   │   │   ├── chat_history_service.py
    │   │   │   ├── exceptions.py
    │   │   │   ├── pubchem_service.py
    │   │   │   ├── quantum_service.py
    │   │   │   ├── settings_service.py
    │   │   │   ├── smiles_service.py
    │   │   │   └── system_service.py
    │   │   ├── tests/
    │   │   │   ├── E2E_TEST_SCENARIOS.md
    │   │   │   ├── README.md
    │   │   │   ├── __init__.py
    │   │   │   ├── conftest.py
    │   │   │   ├── data/
    │   │   │   │   ├── README.md
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── mock_pubchem_response.json
    │   │   │   │   ├── sample_h2.xyz
    │   │   │   │   └── sample_water.xyz
    │   │   │   ├── integration/
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── test_api_endpoints/
    │   │   │   │   │   ├── __init__.py
    │   │   │   │   │   ├── test_agent_api.py
    │   │   │   │   │   ├── test_health_api.py
    │   │   │   │   │   ├── test_pubchem_api.py
    │   │   │   │   │   ├── test_quantum_api.py
    │   │   │   │   │   └── test_smiles_api.py
    │   │   │   │   ├── test_calculation_workflow.py
    │   │   │   │   └── test_websocket_handlers.py
    │   │   │   ├── test_fixtures.py
    │   │   │   └── unit/
    │   │   │       ├── __init__.py
    │   │   │       ├── test_agent/
    │   │   │       │   ├── __init__.py
    │   │   │       │   ├── test_dispatcher_graph.py
    │   │   │       │   ├── test_molecular_agent.py
    │   │   │       │   ├── test_research_agent.py
    │   │   │       │   ├── test_research_tools.py
    │   │   │       │   └── test_tools.py
    │   │   │       ├── test_quantum_calc/
    │   │   │       │   ├── __init__.py
    │   │   │       │   ├── test_dft_calculator.py
    │   │   │       │   └── test_hf_calculator.py
    │   │   │       └── test_services/
    │   │   │           ├── __init__.py
    │   │   │           ├── test_pubchem_service.py
    │   │   │           ├── test_quantum_service.py
    │   │   │           └── test_smiles_service.py
    │   │   └── websocket/
    │   │       ├── __init__.py
    │   │       └── handlers.py
    │   ├── types/
    │   │   ├── 3dmol.d.ts
    │   │   ├── css-modules.d.ts
    │   │   └── electron.d.ts
    │   └── web/
    │       ├── App.css
    │       ├── App.module.css
    │       ├── App.tsx
    │       ├── apiClient.ts
    │       ├── components/
    │       │   ├── AIAgentSwitch.module.css
    │       │   ├── AIAgentSwitch.tsx
    │       │   ├── CIAnalysisViewer.module.css
    │       │   ├── CIAnalysisViewer.tsx
    │       │   ├── ChatHistoryList.module.css
    │       │   ├── ChatHistoryList.tsx
    │       │   ├── ConfirmationModal.module.css
    │       │   ├── ConfirmationModal.tsx
    │       │   ├── DropdownMenu.module.css
    │       │   ├── DropdownMenu.tsx
    │       │   ├── Header.module.css
    │       │   ├── Header.tsx
    │       │   ├── IRSpectrumViewer.module.css
    │       │   ├── IRSpectrumViewer.tsx
    │       │   ├── InlineOrbitalViewer.module.css
    │       │   ├── InlineOrbitalViewer.tsx
    │       │   ├── MolecularOrbitalEnergyDiagram.module.css
    │       │   ├── MolecularOrbitalEnergyDiagram.tsx
    │       │   ├── MolecularOrbitalViewer.module.css
    │       │   ├── MolecularOrbitalViewer.tsx
    │       │   ├── MoleculeViewer.module.css
    │       │   ├── MoleculeViewer.tsx
    │       │   ├── MoleculeViewerSection.module.css
    │       │   ├── MoleculeViewerSection.tsx
    │       │   ├── Sidebar.module.css
    │       │   ├── Sidebar.tsx
    │       │   ├── StyleControls.module.css
    │       │   ├── StyleControls.tsx
    │       │   ├── ToastContainer.module.css
    │       │   ├── ToastContainer.tsx
    │       │   ├── ToastNotification.module.css
    │       │   ├── ToastNotification.tsx
    │       │   ├── XYZInput.module.css
    │       │   └── XYZInput.tsx
    │       ├── data/
    │       │   └── atomicRadii.ts
    │       ├── hooks/
    │       │   ├── index.ts
    │       │   ├── useActiveCalculation.ts
    │       │   ├── useActiveCalculationId.ts
    │       │   ├── useAppSettings.ts
    │       │   ├── useAppState.ts
    │       │   ├── useCalculationActions.ts
    │       │   ├── useCalculationData.ts
    │       │   ├── useCalculationOperations.ts
    │       │   ├── useCalculationQueries.ts
    │       │   ├── useChatHistoryQueries.ts
    │       │   └── useUnifiedWebSocket.ts
    │       ├── index.html
    │       ├── index.tsx
    │       ├── pages/
    │       │   ├── AgentPage.module.css
    │       │   ├── AgentPage.tsx
    │       │   ├── CalculationResultsPage.module.css
    │       │   ├── CalculationResultsPage.tsx
    │       │   ├── CalculationSettingsPage.module.css
    │       │   ├── CalculationSettingsPage.tsx
    │       │   ├── DrawMoleculePage.module.css
    │       │   ├── SettingsPage.module.css
    │       │   └── SettingsPage.tsx
    │       ├── store/
    │       │   ├── agentStore.ts
    │       │   ├── calculationStore.ts
    │       │   ├── chatHistoryStore.ts
    │       │   ├── notificationStore.ts
    │       │   └── uiStore.ts
    │       ├── types/
    │       │   ├── api-types.ts
    │       │   └── generated-api.ts
    │       └── utils/
    │           └── xyzParser.ts
    ├── test_agent_integration.py
    ├── tsconfig.json
    └── webpack.config.ts

## Key API Endpoints
* `GET /health`: Health check endpoint.
* `POST /api/pubchem/search`: Search PubChem.
* `POST /api/smiles/convert`: Convert a SMILES string to XYZ.
* `POST /api/pubchem/validate`: Validate an XYZ format string.
* `POST /api/agent/chat`: Chat with the multi-agent AI system (streams responses via Server-Sent Events). The system automatically routes queries to the appropriate specialist agent (Quantum Calculator, Literature Surveyor, etc.).
* `GET /api/quantum/supported-parameters`: Get lists of supported calculation methods, basis sets, functionals, etc.
* `POST /api/quantum/calculate`: Asynchronously starts a quantum chemistry calculation.
* `GET /api/quantum/calculations`: Lists all saved calculations.
* `GET /api/quantum/calculations/<id>`: Gets detailed results for a specific calculation.
* `PATCH /api/quantum/calculations/<id>`: Updates a calculation's metadata (e.g., renames it).
* `DELETE /api/quantum/calculations/<id>`: Deletes a calculation.
* `POST /api/quantum/calculations/<id>/cancel`: Cancels a running calculation.
* `GET /api/quantum/calculations/<id>/orbitals`: Get molecular orbital information.
* `GET /api/quantum/calculations/<id>/orbitals/{orbitalIndex}/cube`: Generate and retrieve a CUBE file for a specific orbital.
* `GET /api/quantum/calculations/<id>/orbitals/cube-files`: List all generated CUBE files for a calculation.
* `DELETE /api/quantum/calculations/<id>/orbitals/cube-files`: Delete generated CUBE files (all or by index).
* `GET /api/quantum/calculations/<id>/ir-spectrum`: Generate an IR spectrum for a calculation.
* `GET /api/settings`: Get application settings.
* `PUT /api/settings`: Update application settings.
* `GET /api/system/resource-status`: Get system resource status.
* `GET /api/chat-history/sessions`: Get all chat sessions.
* `POST /api/chat-history/sessions`: Create a new chat session.
* `GET /api/chat-history/sessions/<id>`: Get a specific chat session with messages.
* `PATCH /api/chat-history/sessions/<id>`: Update a chat session's metadata (e.g., rename).
* `DELETE /api/chat-history/sessions/<id>`: Delete a chat session.
* `WS /ws/calculations/<id>`: WebSocket endpoint for real-time status updates.

## Parallel Processing Architecture
The application uses a `ProcessPoolExecutor`-based system for quantum chemistry calculations to achieve true multiprocessing parallelism.

* **Process Pool Management**: The `CalculationProcessManager` class manages a pool of worker processes, allowing multiple calculations to run simultaneously.
* **Worker Process Isolation**: Each calculation runs in a separate Python process, eliminating the Global Interpreter Lock (GIL) limitation.
* **Resource Management**: The process pool starts on demand and shuts down cleanly when the Flask server terminates.

## Spin Multiplicity Notes
PySCF uses the `spin` attribute to specify the number of unpaired electrons (2S), which differs from the traditional quantum chemistry notation of 2S+1.

* **Singlet state** (closed-shell): `spin=0`
* **Doublet state** (1 unpaired electron): `spin=1`
* **Triplet state** (2 unpaired electrons): `spin=2`

## Troubleshooting
### Environment Setup Issues
* **Automated Diagnosis**: Always start with environment verification: `npm run verify-env`. This command provides comprehensive diagnostic information.
* **Environment Detection Failures**: If detection fails, use `npm run setup-env` or set the `CONDA_ENV_PATH` environment variable.
* **Build Tool Issues**: If `npm run verify-build-env` fails, ensure `conda-pack` and `gunicorn` are installed in the active `pyscf-env` environment.

### Development Server Issues
The unified Gunicorn-based server reduces environment-specific issues. Most behavior is controlled by `config/server-config.json`. Use `npm run debug:config` to verify server settings. The server automatically finds free ports.

### Build and Packaging Issues
The build system includes pre- and post-build verification steps (`verify-env`, `verify-build-env`, `validate-build`). A failure in these steps indicates an issue with the conda environment or the server configuration.