# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**IMPORTANT**: Your partner is Japanese, so please always report in Japanese.

## Project Overview

PySCF_front is an Electron + React (TypeScript) + Python (Flask) desktop application for molecular visualization and quantum chemistry calculations. It uses PySCF and RDKit for computations, 3Dmol.js for visualization.

**Key Features:**
- Quantum chemistry calculations: DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF
- Geometry optimization and vibrational frequency analysis
- Molecular orbital and NTO visualization
- PubChem/SMILES molecular structure retrieval
- AI chat assistant (Gemini API) for quantum chemistry

**API-First Development**: `src/api-spec/openapi.yaml` is the single source of truth. Run `npm run codegen` to generate TypeScript types and Python Pydantic models.

## Development Philosophy

This is a development-stage application. **Backward compatibility is not a concern.** Prioritize simplicity over compatibility, make breaking changes confidently, and clean up legacy code when better alternatives exist.

## Development Commands

### Quick Start
```bash
npm install                    # Install Node.js dependencies
npm run setup-env              # Automated conda environment setup
npm run verify-env             # Verify environment health
conda activate pyscf-env && npm run dev  # Run development mode
```

### Essential Commands
```bash
npm run dev                    # Development mode (hot reload + Python backend)
npm run build                  # Production build
npm run package                # Package for distribution
npm run codegen                # Generate types from OpenAPI spec
npm run format                 # Format code with Prettier
```

### Python Testing
**IMPORTANT:** Always use the conda environment's Python interpreter.

```bash
cd src/python
conda activate pyscf-env

# Run all tests
python -m pytest tests/ -v

# Run specific test
python -m pytest tests/integration/test_api_endpoints/test_quantum_api.py -v

# Run specific test function
python -m pytest tests/integration/test_api_endpoints/test_quantum_api.py::TestCalculationSubmissionAPI::test_dft_rejects_casci_parameters -xvs

# Run tests matching keyword
python -m pytest tests/ -k "pause_resume" -v
```

## Architecture Overview

### Three-Layer Structure
```
Electron Main Process (src/main.ts)
    ↓ spawns
Python Flask/Gunicorn Backend (src/python/)
    ↓ serves
React Frontend (src/web/) via BrowserWindow
```

### Code Generation (API-First)
- **OpenAPI spec**: `src/api-spec/openapi.yaml`
- **Generated Python models**: `src/python/generated_models.py`
- **Generated TypeScript types**: `src/web/types/generated-api.ts`

### Frontend State Management
- **TanStack Query**: Server-side state (API data, caching)
- **Zustand stores** (`src/web/store/`): UI state
  - `calculationStore.ts`: Active calculation ID, staged calculation
  - `uiStore.ts`: Sidebar visibility, current page, modals
  - `notificationStore.ts`: Toast notifications

### Key Hooks (`src/web/hooks/`)
- `useActiveCalculation`: Unified active calculation state
- `useCalculationQueries`: TanStack Query definitions for quantum API
- `useUnifiedWebSocket`: Real-time updates via WebSocket

### Python Backend Structure (`src/python/`)
- `api/`: Flask route handlers (blueprints)
- `services/`: Business logic layer
- `quantum_calc/`: PySCF calculation implementations
- `websocket/handlers.py`: Real-time status updates

## File Structure
    .
    ├── CLAUDE.md
    ├── Dockerfile
    ├── LICENSE
    ├── PySCF_front_view.png
    ├── README.md
    ├── code.txt
    ├── config
    │   └── server-config.json
    ├── data
    │   ├── chat_history.db
    │   └── settings.json
    ├── package-lock.json
    ├── package.json
    ├── pyscf_front_api.spec
    ├── scripts
    │   ├── build-conda-pack.sh
    │   ├── build-python-linux.sh
    │   ├── bump-version.js
    │   ├── cleanup-artifacts.sh
    │   ├── docker-build-linux.js
    │   ├── docker-run.js
    │   ├── setup-environment.sh
    │   ├── test-python-standalone.js
    │   ├── uitnize.sh
    │   ├── validate-build-completeness.py
    │   └── verify-environment.py
    ├── src
    │   ├── api-spec
    │   │   └── openapi.yaml
    │   ├── assets
    │   │   ├── fonts
    │   │   │   └── ADLaMDisplay-Regular.ttf
    │   │   └── icon
    │   │       ├── linux
    │   │       │   └── icon.png
    │   │       └── mac
    │   │           └── Pyscf_front.icns
    │   ├── main
    │   │   ├── config.ts
    │   │   ├── ipc.ts
    │   │   ├── menu.ts
    │   │   ├── port-manager.ts
    │   │   ├── python-env.ts
    │   │   ├── python-server.ts
    │   │   ├── splash-window-manager.ts
    │   │   └── window-manager.ts
    │   ├── main.ts
    │   ├── preload.ts
    │   ├── python
    │   │   ├── SMILES
    │   │   │   ├── __init__.py
    │   │   │   └── smiles_converter.py
    │   │   ├── __init__.py
    │   │   ├── agent
    │   │   │   └── __init__.py
    │   │   ├── api
    │   │   │   ├── __init__.py
    │   │   │   ├── agent.py
    │   │   │   ├── chat_history.py
    │   │   │   ├── gpu.py
    │   │   │   ├── health.py
    │   │   │   ├── pubchem.py
    │   │   │   ├── quantum.py
    │   │   │   ├── settings.py
    │   │   │   ├── smiles.py
    │   │   │   ├── swagger_ui.py
    │   │   │   └── system.py
    │   │   ├── app.py
    │   │   ├── config.py
    │   │   ├── data
    │   │   │   ├── __init__.py
    │   │   │   ├── scale_factors.py
    │   │   │   └── solvent_properties.py
    │   │   ├── database
    │   │   │   ├── __init__.py
    │   │   │   └── chat_history.py
    │   │   ├── generated_models.py
    │   │   ├── pubchem
    │   │   │   ├── __init__.py
    │   │   │   ├── client.py
    │   │   │   └── parser.py
    │   │   ├── pyscf_front_api.spec
    │   │   ├── pytest.ini
    │   │   ├── quantum_calc
    │   │   │   ├── __init__.py
    │   │   │   ├── base_calculator.py
    │   │   │   ├── casci_calculator.py
    │   │   │   ├── casscf_calculator.py
    │   │   │   ├── ccsd_calculator.py
    │   │   │   ├── config_manager.py
    │   │   │   ├── dft_calculator.py
    │   │   │   ├── exceptions.py
    │   │   │   ├── file_manager.py
    │   │   │   ├── file_watcher.py
    │   │   │   ├── gpu_manager.py
    │   │   │   ├── hf_calculator.py
    │   │   │   ├── ir_spectrum.py
    │   │   │   ├── method_defaults.py
    │   │   │   ├── mp2_calculator.py
    │   │   │   ├── orbital_generator.py
    │   │   │   ├── pause_manager.py
    │   │   │   ├── process_manager.py
    │   │   │   ├── resource_manager.py
    │   │   │   ├── settings_manager.py
    │   │   │   ├── solvent_effects.py
    │   │   │   ├── supported_parameters.py
    │   │   │   └── tddft_calculator.py
    │   │   ├── services
    │   │   │   ├── __init__.py
    │   │   │   ├── chat_history_service.py
    │   │   │   ├── exceptions.py
    │   │   │   ├── gpu_service.py
    │   │   │   ├── notification_service.py
    │   │   │   ├── pubchem_service.py
    │   │   │   ├── quantum_service.py
    │   │   │   ├── settings_service.py
    │   │   │   ├── smiles_service.py
    │   │   │   └── system_service.py
    │   │   ├── tests
    │   │   │   ├── E2E_TEST_SCENARIOS.md
    │   │   │   ├── README.md
    │   │   │   ├── TESTING_IMPLEMENTATION_SUMMARY.md
    │   │   │   ├── __init__.py
    │   │   │   ├── conftest.py
    │   │   │   ├── data
    │   │   │   │   ├── README.md
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── mock_pubchem_response.json
    │   │   │   │   ├── sample_h2.xyz
    │   │   │   │   └── sample_water.xyz
    │   │   │   ├── integration
    │   │   │   │   ├── __init__.py
    │   │   │   │   ├── test_api_endpoints
    │   │   │   │   │   ├── __init__.py
    │   │   │   │   │   ├── test_agent_api.py
    │   │   │   │   │   ├── test_health_api.py
    │   │   │   │   │   ├── test_pubchem_api.py
    │   │   │   │   │   ├── test_quantum_api.py
    │   │   │   │   │   └── test_smiles_api.py
    │   │   │   │   ├── test_auth_production.py
    │   │   │   │   ├── test_auth_security.py
    │   │   │   │   ├── test_calculation_workflow.py
    │   │   │   │   ├── test_pause_resume_workflow.py
    │   │   │   │   └── test_websocket_handlers.py
    │   │   │   ├── test_fixtures.py
    │   │   │   └── unit
    │   │   │       ├── __init__.py
    │   │   │       ├── test_quantum_calc
    │   │   │       │   ├── __init__.py
    │   │   │       │   ├── test_dft_calculator.py
    │   │   │       │   ├── test_hf_calculator.py
    │   │   │       │   └── test_method_defaults.py
    │   │   │       └── test_services
    │   │   │           ├── __init__.py
    │   │   │           ├── test_pubchem_service.py
    │   │   │           ├── test_quantum_service.py
    │   │   │           └── test_smiles_service.py
    │   │   └── websocket
    │   │       ├── __init__.py
    │   │       └── handlers.py
    │   ├── splash
    │   │   ├── preload.ts
    │   │   ├── splash.css
    │   │   ├── splash.html
    │   │   └── splash.ts
    │   ├── types
    │   │   ├── 3dmol.d.ts
    │   │   ├── css-modules.d.ts
    │   │   ├── electron.d.ts
    │   │   ├── ketcher.d.ts
    │   │   └── splash.d.ts
    │   └── web
    │       ├── App.css
    │       ├── App.module.css
    │       ├── App.tsx
    │       ├── apiClient.ts
    │       ├── components
    │       │   ├── AIAgentSwitch.module.css
    │       │   ├── AIAgentSwitch.tsx
    │       │   ├── CIAnalysisViewer.module.css
    │       │   ├── CIAnalysisViewer.tsx
    │       │   ├── ChatHistoryList.module.css
    │       │   ├── ChatHistoryList.tsx
    │       │   ├── ChatMessage.module.css
    │       │   ├── ChatMessage.tsx
    │       │   ├── ConfirmationModal.module.css
    │       │   ├── ConfirmationModal.tsx
    │       │   ├── DropdownMenu.module.css
    │       │   ├── DropdownMenu.tsx
    │       │   ├── Header.module.css
    │       │   ├── Header.tsx
    │       │   ├── IRSpectrumChart.module.css
    │       │   ├── IRSpectrumChart.tsx
    │       │   ├── InitialSetupDialog.module.css
    │       │   ├── InitialSetupDialog.tsx
    │       │   ├── LazyViewer.tsx
    │       │   ├── MolecularOrbitalEnergyDiagram.module.css
    │       │   ├── MolecularOrbitalEnergyDiagram.tsx
    │       │   ├── MolecularOrbitalViewer.module.css
    │       │   ├── MolecularOrbitalViewer.tsx
    │       │   ├── MoleculeViewer.module.css
    │       │   ├── MoleculeViewer.tsx
    │       │   ├── MoleculeViewerSection.module.css
    │       │   ├── MoleculeViewerSection.tsx
    │       │   ├── MullikenChargeViewer.module.css
    │       │   ├── MullikenChargeViewer.tsx
    │       │   ├── Sidebar.module.css
    │       │   ├── Sidebar.tsx
    │       │   ├── StyleControls.module.css
    │       │   ├── StyleControls.tsx
    │       │   ├── ToastContainer.module.css
    │       │   ├── ToastContainer.tsx
    │       │   ├── ToastNotification.module.css
    │       │   ├── ToastNotification.tsx
    │       │   ├── VibrationModeViewer.module.css
    │       │   ├── VibrationModeViewer.tsx
    │       │   ├── XYZInput.module.css
    │       │   └── XYZInput.tsx
    │       ├── data
    │       │   └── atomicRadii.ts
    │       ├── hooks
    │       │   ├── index.ts
    │       │   ├── useActiveCalculation.ts
    │       │   ├── useActiveCalculationId.ts
    │       │   ├── useAppSettings.ts
    │       │   ├── useAppState.ts
    │       │   ├── useCalculationActions.ts
    │       │   ├── useCalculationData.ts
    │       │   ├── useCalculationOperations.ts
    │       │   ├── useCalculationQueries.ts
    │       │   ├── useChatHistoryQueries.ts
    │       │   ├── useMethodDefaults.ts
    │       │   ├── useProcessedCalculationResults.ts
    │       │   └── useUnifiedWebSocket.ts
    │       ├── index.html
    │       ├── index.tsx
    │       ├── pages
    │       │   ├── AgentPage.module.css
    │       │   ├── AgentPage.tsx
    │       │   ├── CalculationResultsPage.module.css
    │       │   ├── CalculationResultsPage.tsx
    │       │   ├── CalculationSettingsPage.module.css
    │       │   ├── CalculationSettingsPage.tsx
    │       │   ├── DrawMoleculePage.module.css
    │       │   ├── DrawMoleculePage.tsx
    │       │   ├── SettingsPage.module.css
    │       │   └── SettingsPage.tsx
    │       ├── store
    │       │   ├── agentStore.ts
    │       │   ├── calculationStore.ts
    │       │   ├── chatHistoryStore.ts
    │       │   ├── notificationStore.ts
    │       │   └── uiStore.ts
    │       ├── types
    │       │   ├── api-types.ts
    │       │   └── generated-api.ts
    │       └── utils
    │           ├── dateFormatter.ts
    │           ├── errorHandler.ts
    │           ├── irSpectrumConstants.ts
    │           └── xyzParser.ts
    ├── tsconfig.json
    └── webpack.config.ts

## Key Concepts

### API Endpoints
See `src/api-spec/openapi.yaml` for complete API documentation. In development mode, Swagger UI is available at `http://127.0.0.1:5000/api-docs/`.

### Parallel Processing
Quantum calculations use `ProcessPoolExecutor` for true multiprocessing (bypasses Python GIL). See `quantum_calc/process_manager.py`.

### Spin Multiplicity (PySCF Convention)
PySCF uses `spin` = number of unpaired electrons (2S), NOT the traditional 2S+1 notation:
- Singlet: `spin=0`
- Doublet: `spin=1`
- Triplet: `spin=2`

### Server Configuration
`config/server-config.json` controls server behavior (Gunicorn workers, threads, timeouts). Both development and production use Gunicorn for consistency.

## Troubleshooting

Run `npm run verify-env` for automated environment diagnosis. If detection fails, use `npm run setup-env` or set the `CONDA_ENV_PATH` environment variable.