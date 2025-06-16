# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PySCF_Front is a desktop quantum chemistry calculation application built entirely in Python:
- **GUI Framework**: PySide6 (Qt 6.7+) for native cross-platform desktop interface
- **Computational Engine**: PySCF 2.9+ for quantum chemistry calculations
- **Database**: MySQL 8.0+ with SQLAlchemy ORM for data persistence
- **Visualization**: VTK and Matplotlib integration for 3D molecular visualization
- **Optional MCP Server**: Model Context Protocol server for Claude Desktop integration

## Architecture

### Unified Python Architecture
- **Single Language**: All components written in Python for seamless integration
- **Direct Integration**: No gRPC needed - direct function calls between GUI and computation
- **Native Performance**: Qt-based GUI with OpenGL/VTK acceleration
- **Scientific Ecosystem**: Direct integration with NumPy, SciPy, Matplotlib, RDKit

### Core Components
- **PySide6 GUI**: Main window with dockable panels, 3D viewer, and dialogs
- **Calculation Engine**: PySCF-based quantum chemistry computations with job queue
- **Database Layer**: SQLAlchemy models for molecules, calculations, and results
- **Plugin System**: Extensible architecture for methods, basis sets, and analysis tools
- **3D Visualization**: VTK-based molecular viewer with interactive rendering

### Data Flow
1. User inputs molecule via PySide6 GUI dialogs
2. Direct Python function calls to calculation engine
3. QThreadPool manages background calculations with progress signals
4. Results stored via SQLAlchemy to MySQL database
5. VTK/OpenGL renders 3D molecular structures and properties

## Development Commands

### Environment Setup
```bash
# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
pip install -r requirements-dev.txt

# Development mode installation
pip install -e .

# Test environment setup
python test_environment.py
```

### Application Development
```bash
# Run development version
python run_dev.py

# Run with debug logging
PYSCF_FRONT_DEBUG=1 python pyscf_front/main.py

# Run tests
pytest pyscf_front/tests/ --cov=pyscf_front

# Run specific GUI tests
pytest pyscf_front/tests/ --qt-test-timeout=5000

# Run tests with verbose output
pytest -v --tb=short pyscf_front/tests/

# Code formatting
black pyscf_front/ --line-length=88 --target-version=py313

# Type checking
mypy pyscf_front/

# Linting
flake8 pyscf_front/
```

### UI Development
```bash
# Launch Qt Designer for UI editing
pyside6-designer

# Convert UI files to Python (if needed)
pyside6-uic pyscf_front/gui/resources/ui_files/main_window.ui -o pyscf_front/gui/ui_main_window.py

# Update translations
pyside6-lupdate pyscf_front/main.py -ts pyscf_front/translations/pyscf_front_ja.ts
pyside6-lrelease pyscf_front/translations/pyscf_front_ja.ts
```

### Database Management
```bash
# Create database migration
alembic revision --autogenerate -m "description"

# Apply migrations
alembic upgrade head

# Reset database (development)
python -c "from pyscf_front.database.models import Base; from pyscf_front.database.connection import engine; Base.metadata.drop_all(engine); Base.metadata.create_all(engine)"
```

### Packaging and Distribution
```bash
# Build executable with PyInstaller
pyinstaller pyscf_front.spec

# Create distributable archive
python setup.py sdist bdist_wheel

# Install from source
pip install dist/pyscf_front-*.whl
```

## Key Directories

**For detailed project structure documentation, see [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)**

- `gui/`: PySide6 GUI components and main window
- `gui/widgets/`: Reusable UI widgets (molecule viewer, calculation panel, results display)
- `gui/dialogs/`: Modal dialogs for settings, molecule input, about
- `gui/resources/`: UI resources (icons, styles, Qt Designer files)
- `core/`: Business logic and computational engine
- `core/calculation_engine.py`: PySCF integration and job management
- `database/`: SQLAlchemy models and database operations
- `database/models.py`: ORM definitions for molecules, calculations, results
- `services/`: Service layer for business logic (calculation, molecule, instance services)
- `plugins/`: Plugin system for extending functionality
- `utils/`: Utility functions (config, logging, validation)
- `mcp_server/`: Optional MCP server for Claude Desktop integration
- `tests/`: Test suites including GUI tests with pytest-qt

## Core API Structure

The application uses direct Python API calls between GUI and computation engine:

### Calculation Engine API
- `CalculationEngine.submit_calculation()`: Submit calculation jobs to queue
- `CalculationEngine.get_job_status()`: Check calculation progress
- `CalculationEngine.get_results()`: Retrieve completed results
- `CalculationEngine.cancel_job()`: Cancel running calculations

### Molecule Manager API
- `MoleculeManager.create_from_smiles()`: Create molecule from SMILES string
- `MoleculeManager.create_from_xyz()`: Import from XYZ coordinate file
- `MoleculeManager.optimize_geometry()`: Perform geometry optimization
- `MoleculeManager.validate_structure()`: Validate molecular structure

### Database Repository API
- `Repository.save_molecule()`: Persist molecule to database
- `Repository.save_calculation()`: Store calculation parameters and results
- `Repository.get_calculation_history()`: Retrieve calculation history
- `Repository.search_molecules()`: Search stored molecules

## MCP Server Integration

The optional MCP server enables natural language interaction with PySCF_Front through Claude Desktop:
- Enable with `MCP_SERVER_ENABLED=true` environment variable
- Provides tools for molecule creation, calculation recommendations, and result interpretation
- MCP server configuration available in `.env` file
- Security features include IP restrictions and rate limiting

## Database Schema

SQLAlchemy models define the following key tables:
- `instances`: Project containers for organizing work
- `molecules`: Molecular structure data (atoms, coordinates, properties)
- `calculations`: Calculation job information (method, basis set, parameters)
- `results`: Computed results (energies, orbitals, properties)
- `jobs`: Job status and progress tracking with enum states

## PySide6-Specific Notes

### Threading and Signals
- Use `QThreadPool` for background calculations to keep GUI responsive
- Implement `QRunnable` workers with custom `QObject` signals for progress updates
- Connect calculation progress to GUI progress bars via Qt signals/slots

### 3D Visualization
- `QOpenGLWidget` for high-performance molecular rendering
- VTK integration via `QVTKRenderWindowInteractor` for advanced visualization
- Matplotlib embedding in Qt widgets for plots and graphs

### GUI Architecture
- Main window uses dockable panels (`QDockWidget`) for flexible layout
- Custom widgets inherit from appropriate Qt base classes (`QWidget`, `QDialog`)
- Qt Designer integration for visual UI design with `.ui` files

### Testing
- Use `pytest-qt` for GUI testing with `qtbot` fixtures
- Mock calculation engines for GUI testing without heavy computations
- Test signal/slot connections and widget interactions

## PySCF Integration Notes

- Direct Python integration eliminates communication overhead
- Supports HF, DFT (B3LYP, PBE), and post-HF methods (MP2, CCSD)
- GPU acceleration via GPU4PySCF when CUDA is available
- Real-time progress tracking through Python generator functions
- Results processing with NumPy/SciPy for analysis and visualization

## Environment Configuration

### Development Environment Variables
Create a `.env` file based on `.env.example` for local development:
```bash
# Database configuration
DB_HOST=localhost
DB_PORT=3306
DB_NAME=pyscf_front
DB_USER=your_username
DB_PASSWORD=your_password

# Application settings
PYSCF_FRONT_DEBUG=true
PYSCF_FRONT_LOG_LEVEL=DEBUG

# GUI settings
PYSCF_FRONT_THEME=dark
PYSCF_FRONT_LANGUAGE=en

# Testing
PYTEST_QT_API=pyside6
```

### VS Code Integration
The project includes VS Code configuration for optimal development:
- **Python interpreter**: Automatically uses local `venv/bin/python`
- **Code formatting**: Black with line length 88, Python 3.13 target
- **Linting**: Flake8 and MyPy integration
- **Qt tooling**: PySide6 tools configured (designer, uic, rcc)
- **Launch configurations**: Debug and run configurations for both dev and production modes

### Development Tools Configuration
- **Black**: `--line-length=88 --target-version=py313`
- **MyPy**: Strict type checking enabled
- **Flake8**: Standard PEP 8 compliance
- **Pytest**: Qt GUI testing with 5-second timeout for GUI tests

# important-instruction-reminders
Do what has been asked; nothing more, nothing less.
NEVER create files unless they're absolutely necessary for achieving your goal.
ALWAYS prefer editing an existing file to creating a new one.
NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.