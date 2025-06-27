# PySCF Front

A cross-platform quantum chemistry calculation application that provides an intuitive GUI frontend for PySCF (Python-based Simulations of Chemistry Framework).

## Features

- **Intuitive GUI**: Modern Electron-based interface with React
- **Powerful Backend**: PySCF integration for accurate quantum chemistry calculations
- **Cross-Platform**: Works on Windows, macOS, and Linux
- **Real-time 3D Visualization**: Molecular structure and orbital visualization (3Dpymol integration planned)
- **Multiple Calculation Methods**: Support for HF, DFT, MP2, and more
- **Project Management**: Organize calculations and results
- **Batch Processing**: Queue multiple calculations

## Technology Stack

- **Frontend**: Electron + React + TypeScript
- **Backend**: Python 3.13 + PySCF
- **Visualization**: 3Dpymol (planned)
- **Styling**: Custom CSS with modern design system

## Prerequisites

- Node.js >= 18.0.0
- Python 3.13
- npm >= 8.0.0

## Installation & Setup

### 1. Clone and Setup

```bash
git clone <repository-url>
cd pyscf-front
```

### 2. Install Dependencies

```bash
# Install Node.js dependencies
npm install

# Set up Python virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt
```

### 3. Development

```bash
# Start development mode (runs both frontend and backend)
npm run dev

# Or run components separately:
npm run dev:main     # Electron main process (TypeScript compilation)
npm run dev:renderer # React frontend (Vite dev server)
```

### 4. Building

```bash
# Build for production
npm run build

# Run the built application
npm start

# Create distributable packages
npm run dist
```

## Development Commands

| Command | Description |
|---------|-------------|
| `npm run dev` | Start development mode with hot reload |
| `npm run build` | Build for production |
| `npm start` | Run the built application |
| `npm run dist` | Create distributable packages |
| `npm test` | Run tests |
| `npm run lint` | Run ESLint |
| `npm run typecheck` | Run TypeScript type checking |

## Project Structure

```
pyscf-front/
├── src/
│   ├── main/                 # Electron main process
│   │   ├── index.ts         # Entry point
│   │   └── preload.ts       # Context bridge
│   ├── renderer/            # React frontend
│   │   ├── components/      # UI components
│   │   ├── App.tsx         # Main app component
│   │   └── main.tsx        # React entry point
│   └── python/             # Python backend
│       ├── main.py         # Backend entry point
│       ├── calculations/   # Calculation engine
│       └── utils/          # Utilities
├── assets/                 # Static resources
├── config/                 # Configuration files
├── data/                   # Data storage
└── tests/                  # Test files
```

## Usage

1. **Start the Application**: Run `npm run dev` for development
2. **Load a Molecule**: Select from sample molecules or load your own
3. **Configure Calculation**: Choose method, functional, and basis set
4. **Run Calculation**: Click "Start Calculation" and monitor progress
5. **View Results**: Results appear in the calculation panel with detailed information

## Supported Calculation Methods

- **Hartree-Fock (HF)**: RHF, UHF, ROHF
- **Density Functional Theory (DFT)**: B3LYP, PBE, PBE0, M06-2X, and more
- **Møller-Plesset (MP2)**: Second-order perturbation theory
- **Future**: CCSD(T), CASSCF, TD-DFT

## Supported Basis Sets

- Pople series: STO-3G, 6-31G, 6-31G*, 6-311G*, etc.
- Dunning series: cc-pVDZ, cc-pVTZ, aug-cc-pVDZ, etc.
- Custom basis sets (planned)

## Input Formats

- **Coordinates**: Direct atomic coordinates input
- **XYZ Files**: Standard XYZ format
- **SMILES**: Chemical notation (requires RDKit - planned)
- **PDB**: Protein Data Bank format (planned)

## Architecture

The application uses a multi-process architecture:

1. **Electron Main Process**: Manages the application lifecycle and window
2. **Electron Renderer Process**: React-based UI
3. **Python Backend**: PySCF calculations via child process communication

Communication flows through IPC (Inter-Process Communication) between Electron processes and JSON messaging with the Python backend.

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make your changes
4. Run tests: `npm test`
5. Submit a pull request

## Development Notes

- The Python backend communicates with Electron via stdin/stdout
- All calculation requests are queued and processed sequentially
- Results are cached locally in JSON format
- The application follows the design patterns specified in `claude.md`

## Future Enhancements

- 3Dpymol integration for molecular visualization
- RDKit integration for SMILES processing
- MCP server implementation for Claude Desktop integration
- GPU acceleration support
- Cloud calculation backends
- Machine learning integration

## License

MIT License - see LICENSE file for details

## Support

For issues and questions, please use the GitHub issue tracker.