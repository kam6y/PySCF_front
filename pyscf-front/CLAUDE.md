# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PySCF_Front is a quantum chemistry calculation application with a **microservices architecture** consisting of:
- **Flutter frontend** (desktop GUI for molecular visualization and calculation setup)
- **Python backend** (quantum chemistry computation engine using PySCF)
- **MySQL database** (calculation data persistence)
- **Optional MCP server** (Claude Desktop integration for natural language interfaces)

## Architecture Patterns

### Core Data Flow
```
Instance (molecule + calculation settings) → Job Queue → PySCF Engine → Results Storage
```

### Key Architectural Decisions
- **Instance-Centric Design**: All calculations are organized around "instances" that combine molecular structure with computation parameters
- **Plugin-Based Calculation Methods**: Quantum chemistry methods (HF, DFT, MP2, etc.) are implemented as plugins in `backend/src/plugins/`
- **Dual API Pattern**: Both gRPC (for Flutter) and FastAPI (for debugging/tools) endpoints
- **Asynchronous Job Processing**: Long-running calculations use MySQL job queue with status tracking

### Database Schema Relationships
```
instances (1) → (N) molecules
instances (1) → (N) calculations → (N) results
calculations (1) → (1) job_queue
```

## Development Commands

### Environment Setup
```bash
make setup              # Complete environment setup (both backend + frontend)
make setup-backend      # Python venv + dependencies only
make setup-frontend     # Flutter dependencies + code generation
```

### Docker Development (Preferred)
```bash
make docker-up-dev      # Start development containers with hot reload
make docker-up          # Start production containers
make docker-down        # Stop all containers
make docker-build       # Rebuild container images
```

### Local Development
```bash
make run-backend        # Start Python API server (localhost:8000)
make run-frontend       # Start Flutter desktop app
```

### Testing and Quality
```bash
make test               # Run full test suite (pytest + flutter test)
make clean              # Clean cache files

# Individual testing
cd backend && pytest tests/              # Backend unit tests
cd backend && pytest -m integration     # Integration tests only
cd frontend && flutter test             # Flutter widget tests
```

### Backend-Specific Commands
```bash
cd backend
uvicorn src.main:app --reload           # Start FastAPI server
python -m pytest tests/test_*.py        # Run specific test file
black src/                              # Format Python code
mypy src/                               # Type checking
```

### Frontend-Specific Commands
```bash
cd frontend
flutter run -d macos                    # Run on macOS
flutter run -d linux                    # Run on Linux
flutter run -d windows                  # Run on Windows
flutter build linux --release           # Build production binary
dart run build_runner build             # Generate code (protobuf, etc.)
```

## Project Structure Logic

### Backend (`backend/src/`)
- `main.py` - FastAPI application entry point with health endpoints
- `core/` - PySCF configuration and computation engine setup
- `plugins/` - Extensible quantum chemistry method implementations
  - `methods/` - Calculation methods (HF, DFT, post-HF)
  - `basis_sets/` - Basis set definitions
  - `analysis/` - Post-calculation analysis tools
- `api/` - gRPC service implementations (to be implemented)

### Frontend (`frontend/lib/`)
- Planned structure (not yet implemented):
  - `models/` - Data models matching backend entities
  - `views/` - UI screens for molecule input, calculation setup, results
  - `services/` - gRPC client communication
  - `widgets/` - Reusable UI components

### Database (`database/`)
- `schema/001_initial_schema.sql` - Complete database structure
- Core entities: instances → molecules/calculations → results
- Sample data included for testing

## Development Workflow

### Adding New Calculation Methods
1. Create plugin in `backend/src/plugins/methods/`
2. Inherit from `CalculationPlugin` base class
3. Implement: `validate_input()`, `run_calculation()`, `parse_results()`
4. Register in plugin loader

### Database Changes
1. Create migration script in `database/migrations/`
2. Update schema file if needed
3. Test with sample data

### Frontend-Backend Integration
1. Define gRPC service in `protos/`
2. Generate code: `protoc --dart_out=frontend/lib/generated`
3. Implement service in `backend/src/api/`
4. Create client in `frontend/lib/services/`

## Configuration Details

### Environment Variables (.env)
```bash
DATABASE_URL=mysql+pymysql://pyscf_user:pyscf_password@localhost:3306/pyscf_dev
PYSCF_USE_GPU=false                 # Enable GPU acceleration
MCP_SERVER_ENABLED=false            # Enable Claude Desktop integration
```

### Docker Service Ports
- MySQL: 3306
- Backend API: 8000 (FastAPI) + 50051 (gRPC)
- phpMyAdmin: 8080
- Frontend: 8080 (development server)

### Python Dependencies Management
- Core dependencies in `backend/requirements.txt`
- Optional dependencies in `pyproject.toml`:
  - `[dev]` - Testing and code quality tools
  - `[gpu]` - CUDA acceleration (gpu4pyscf)
  - `[mcp]` - Claude Desktop integration

## Code Quality Standards

### Python (Backend)
- **Formatting**: Black (88 char line length)
- **Import sorting**: isort (black-compatible)
- **Type checking**: mypy with strict settings
- **Linting**: flake8
- **Testing**: pytest with async support

### Flutter (Frontend)
- **Linting**: flutter_lints package
- **Testing**: flutter_test + mockito for mocking
- **Code generation**: build_runner for protobuf

## Key Implementation Notes

### PySCF Integration
- Configuration in `backend/src/core/pyscf_config.py`
- Memory and threading settings via environment variables
- GPU detection and optional acceleration setup

### Database Design Patterns
- UUID primary keys for distributed scalability
- JSON columns for flexible parameter storage
- Foreign key cascading for data integrity
- Indexed fields for query performance

### MCP Server (Optional Feature)
- Disabled by default (`MCP_SERVER_ENABLED=false`)
- Provides natural language interface to quantum chemistry calculations
- Independent server that exposes PySCF_Front functionality to Claude Desktop
- Tools for molecule generation, calculation recommendations, results interpretation