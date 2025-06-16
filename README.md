# PySCF_native_app

A modern desktop quantum chemistry calculation application built entirely in Python with PySide6 and PySCF.

## 🚀 Features

- **Native Cross-Platform GUI**: Built with PySide6 (Qt 6.7+) for Windows, macOS, and Linux
- **Quantum Chemistry Calculations**: Powered by PySCF 2.9+ with support for HF, DFT, and post-HF methods
- **3D Molecular Visualization**: Interactive molecular viewer with VTK and Matplotlib integration
- **Database Persistence**: SQLAlchemy ORM with MySQL/SQLite for calculation history and results
- **Asynchronous Computing**: Background calculations with real-time progress tracking
- **Extensible Plugin System**: Modular architecture for custom methods and analysis tools
- **MCP Server Integration**: Optional Claude Desktop integration for natural language interactions

## 📋 Requirements

### System Requirements
- Python 3.11+ (recommended 3.13)
- 4GB+ RAM (8GB+ recommended for larger calculations)
- OpenGL 3.3+ compatible graphics card
- 2GB+ free disk space

### Dependencies
See `requirements.txt` for the complete list. Key dependencies include:
- PySide6 6.7.0+ (Qt GUI framework)
- PySCF 2.9.0+ (quantum chemistry engine)
- RDKit 2024.03.3+ (molecular processing)
- SQLAlchemy 2.0.30+ (database ORM)
- VTK 9.3.0+ (3D visualization)
- NumPy, SciPy, Matplotlib (scientific computing)

## 🛠️ Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/PySCF_native_app.git
cd PySCF_native_app
```

### 2. Create Virtual Environment
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 3. Install Dependencies
```bash
pip install --upgrade pip
pip install -r requirements.txt
pip install -r requirements-dev.txt  # For development
```

### 4. Development Installation
```bash
pip install -e .
```

### 5. Verify Installation
```bash
python test_environment.py
```

## 🎮 Usage

### Quick Start
```bash
# Run the application
python run_dev.py

# Or with debug logging
PYSCF_FRONT_DEBUG=1 python pyscf_front/main.py
```

### Basic Workflow
1. **Create Molecule**: Input via SMILES, XYZ file, or molecular presets
2. **Setup Calculation**: Choose method (HF, B3LYP, etc.) and basis set
3. **Run Calculation**: Monitor progress in real-time
4. **Analyze Results**: View energies, orbitals, and molecular properties
5. **Save Project**: All data automatically saved to database

### Supported Input Formats
- **SMILES strings**: `CCO` (ethanol), `c1ccccc1` (benzene)
- **XYZ coordinate files**: Standard molecular coordinate format
- **Molecular presets**: Water, methane, benzene, and more

### Calculation Methods
- **Hartree-Fock (HF)**: Basic quantum chemistry method
- **DFT Methods**: B3LYP, PBE, and other density functional theory methods
- **Basis Sets**: STO-3G, 6-31G, 6-31G(d), cc-pVDZ, and more
- **Post-HF**: MP2, CCSD (advanced methods)

## 🏗️ Development

### Environment Setup
```bash
# Code formatting
black pyscf_front/ tests/ --line-length=88 --target-version=py313

# Type checking
mypy pyscf_front/

# Linting
flake8 pyscf_front/ tests/

# Run tests
pytest pyscf_front/tests/ -v --cov=pyscf_front
```

### Project Structure
See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for detailed architecture documentation.

### Key Development Commands
```bash
# Run tests
pytest pyscf_front/tests/ -v

# GUI tests specifically
pytest pyscf_front/tests/ --qt-test-timeout=5000

# Launch Qt Designer
pyside6-designer

# Database migrations (if using Alembic)
alembic upgrade head
```

## 📊 Database Configuration

### SQLite (Default for Development)
No additional setup required. Database file created automatically.

### MySQL (Production)
1. Install MySQL 8.0+
2. Create database and user
3. Configure connection in `.env` file:
```bash
DB_HOST=localhost
DB_PORT=3306
DB_NAME=pyscf_front
DB_USER=your_username
DB_PASSWORD=your_password
```

## 🔧 Configuration

### Environment Variables
Create a `.env` file (see `.env.example`):
```bash
# Application settings
PYSCF_FRONT_DEBUG=false
PYSCF_FRONT_LOG_LEVEL=INFO
PYSCF_FRONT_THEME=dark

# Database
DB_HOST=localhost
DB_NAME=pyscf_front

# GUI settings
PYSCF_FRONT_LANGUAGE=en
```

### MCP Server (Optional)
Enable Claude Desktop integration:
```bash
MCP_SERVER_ENABLED=true
MCP_SERVER_PORT=8000
```

## 📚 Documentation

- [Project Structure](PROJECT_STRUCTURE.md) - Detailed architecture overview
- [Development Guide](CLAUDE.md) - Development instructions and API reference
- [Implementation Status](IMPLEMENTATION_STATUS.md) - Current feature completion
- [MCP Server Setup](mcp-server-setup.md) - Claude Desktop integration guide

## 🧪 Testing

The project includes comprehensive test coverage:

```bash
# Run all tests
pytest pyscf_front/tests/ -v

# Run with coverage report
pytest pyscf_front/tests/ --cov=pyscf_front --cov-report=html

# Run specific test categories
pytest pyscf_front/tests/test_services.py -v
pytest pyscf_front/tests/test_calculation_engine_integration.py -v
```

### Test Categories
- **Unit Tests**: Individual component testing
- **Integration Tests**: Service layer and database integration
- **GUI Tests**: PySide6 interface testing with pytest-qt
- **Calculation Tests**: PySCF integration testing (with mocking)

## 🚀 Deployment

### Building Executables
```bash
# Create standalone executable
pyinstaller pyscf_front.spec

# Create distributable package
python setup.py sdist bdist_wheel
```

### Docker Support
```bash
# Build Docker image
docker build -t pyscf-native-app .

# Run in container
docker run -v $(pwd)/data:/app/data pyscf-native-app
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests (`pytest pyscf_front/tests/`)
5. Format code (`black`, `flake8`, `mypy`)
6. Commit changes (`git commit -m 'Add amazing feature'`)
7. Push to branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

### Development Standards
- Follow PEP 8 style guidelines
- Maintain type hints throughout codebase
- Write comprehensive tests for new features
- Update documentation for API changes
- Use conventional commit messages

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **PySCF Team**: For the excellent quantum chemistry library
- **Qt/PySide6**: For the robust GUI framework
- **RDKit**: For molecular processing capabilities
- **Scientific Python Community**: NumPy, SciPy, Matplotlib ecosystems
- **VTK**: For advanced 3D visualization capabilities

## 📞 Support

- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/yourusername/PySCF_native_app/issues)
- 💡 **Feature Requests**: [GitHub Discussions](https://github.com/yourusername/PySCF_native_app/discussions)
- 📧 **Email**: your.email@example.com
- 📖 **Documentation**: [Wiki](https://github.com/yourusername/PySCF_native_app/wiki)

## 🏆 Project Status

**Current Version**: 0.9.5 (95% Complete)

### ✅ Completed Features
- Molecular input system (SMILES, XYZ, presets)
- PySCF calculation engine integration
- GUI with job management and progress tracking
- Database persistence layer
- Service architecture
- Comprehensive testing suite
- Basic 3D molecular visualization

### 🚧 In Progress
- Advanced 3D visualization features
- Additional analysis tools
- Enhanced user documentation
- Performance optimizations

### 🔮 Planned Features
- GPU acceleration support
- Advanced plotting and analysis
- Batch calculation processing
- Plugin marketplace
- Cloud calculation integration

---

**Made with ❤️ by the quantum chemistry community**