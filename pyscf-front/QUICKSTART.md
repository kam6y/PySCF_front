# Quick Start Guide

## Development Environment Setup Complete! 🎉

Your PySCF Front development environment has been successfully set up with:

### ✅ Completed Setup
1. **Project Structure**: Complete directory structure matching specifications
2. **Node.js Environment**: Electron + React + TypeScript configured
3. **Python Environment**: Virtual environment with PySCF and dependencies
4. **Build System**: Vite + TypeScript compilation working
5. **Development Tools**: ESLint, Prettier, Jest configured
6. **IPC Communication**: Electron ↔ Python backend ready
7. **Basic UI**: React components with modern design system

### 🚀 Quick Start Commands

```bash
# Start development mode (recommended)
npm run dev

# Run Python backend tests
source venv/bin/activate && python -c "
import sys; sys.path.insert(0, 'src/python')
from utils.molecule_builder import MoleculeBuilder
mb = MoleculeBuilder()
mol = mb.build_from_data(mb.get_test_molecules()['water'])
print(f'✓ Backend working: {mol.natm} atoms loaded')
"

# Build for production
npm run build

# Run type checking
npm run typecheck

# Run linting
npm run lint
```

### 🧪 Test the Application

1. **Start Development Mode**:
   ```bash
   npm run dev
   ```

2. **The app should open with**:
   - Left panel: Sample molecules (Water, Methane, Hydrogen)
   - Center: 3D viewer placeholder
   - Right panel: Calculation controls
   - Bottom: Status bar with logs

3. **Try a Calculation**:
   - Click on "Water (H₂O)" in the left panel
   - Ensure method is set to "HF" and basis to "STO-3G"
   - Click "Start Calculation"
   - Watch the results appear in the right panel

### 📁 Key Files

- `src/main/index.ts` - Electron main process
- `src/python/main.py` - Python backend entry point
- `src/renderer/App.tsx` - React main component
- `package.json` - Node.js dependencies and scripts
- `requirements.txt` - Python dependencies

### 🔧 Development Workflow

1. **Frontend Changes**: Hot reload via Vite (automatic)
2. **Backend Changes**: Restart with `npm run dev`
3. **Main Process Changes**: Restart Electron
4. **Type Checking**: `npm run typecheck`
5. **Code Style**: `npm run lint:fix`

### 🧬 Supported Features

- **Molecules**: Water, Methane, Hydrogen (sample)
- **Methods**: HF, DFT, MP2
- **Functionals**: B3LYP, PBE, PBE0, M06-2X, wB97X-D
- **Basis Sets**: STO-3G, 6-31G*, cc-pVDZ, and more

### 📊 What Works Now

- ✅ Molecule loading and building
- ✅ Quantum chemistry calculations (HF, DFT, MP2)
- ✅ Real-time calculation status
- ✅ Results display with energy and orbital info
- ✅ Cross-platform Electron app
- ✅ Modern React UI with TypeScript

### 🔮 Next Steps

- **3Dpymol Integration**: Add real molecular visualization
- **RDKit Integration**: SMILES input support
- **File I/O**: XYZ, PDB file loading
- **Advanced Calculations**: Geometry optimization, frequency analysis
- **Project Management**: Save/load calculation projects
- **MCP Server**: Claude Desktop integration

### 🐛 Troubleshooting

**Python Backend Issues**:
```bash
source venv/bin/activate
python src/python/main.py
# Should wait for JSON input
```

**Build Issues**:
```bash
npm run typecheck  # Check TypeScript errors
npm run lint       # Check code style
```

**Dependencies**:
```bash
npm install        # Reinstall Node.js deps
pip install -r requirements.txt  # Reinstall Python deps
```

### 🏗️ Architecture Overview

```
┌─────────────────────┐
│  Electron Frontend  │ ← React + TypeScript
│  (Port 3000)       │
└──────────┬──────────┘
           │ IPC
┌──────────┴──────────┐
│ Electron Main       │ ← Node.js + Electron
│ (Process Manager)   │
└──────────┬──────────┘
           │ Child Process + JSON
┌──────────┴──────────┐
│   Python Backend    │ ← PySCF + Scientific Python
│     (Calculations)  │
└─────────────────────┘
```

**Congratulations! Your PySCF Front development environment is ready for quantum chemistry calculations! 🎉**