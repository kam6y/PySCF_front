"""File management utilities for quantum chemistry calculations."""

import os
import shutil
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime


class CalculationFileManager:
    """Manages files generated during quantum chemistry calculations."""
    
    def __init__(self, base_dir: Optional[str] = None):
        """Initialize file manager with base directory."""
        if base_dir is None:
            # Use user's home directory for calculations
            home = Path.home()
            self.base_dir = home / "PySCF_Calculations"
        else:
            self.base_dir = Path(base_dir)
        
        # Create base directory if it doesn't exist
        self.base_dir.mkdir(exist_ok=True)
    
    def create_calculation_dir(self, molecule_name: Optional[str] = None) -> str:
        """Create a new directory for calculation files."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if molecule_name:
            # Clean molecule name for use in directory name
            clean_name = "".join(c for c in molecule_name if c.isalnum() or c in "._-")
            dir_name = f"{clean_name}_{timestamp}"
        else:
            dir_name = f"calculation_{timestamp}"
        
        calc_dir = self.base_dir / dir_name
        calc_dir.mkdir(exist_ok=True)
        
        return str(calc_dir)
    
    def get_checkpoint_path(self, calc_dir: str, filename: str = "calculation.chk") -> str:
        """Get the path to a checkpoint file in the calculation directory."""
        return os.path.join(calc_dir, filename)
    
    def save_calculation_info(self, calc_dir: str, info: Dict[str, Any]) -> None:
        """Save calculation information to a text file."""
        info_file = Path(calc_dir) / "calculation_info.txt"
        
        with open(info_file, 'w') as f:
            f.write("PySCF Quantum Chemistry Calculation\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Calculation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Molecule: {info.get('molecule_name', 'Unknown')}\n")
            f.write(f"Method: {info.get('calculation_method', 'DFT')}\n")
            f.write(f"Basis Set: {info.get('basis', 'N/A')}\n")
            f.write(f"XC Functional: {info.get('xc_functional', 'N/A')}\n")
            f.write(f"Charge: {info.get('charge', 0)}\n")
            f.write(f"Spin Multiplicity: {info.get('spin_multiplicity', 1)}\n")
            f.write(f"SCF Energy: {info.get('scf_energy', 'N/A')} hartree\n")
            f.write(f"Converged: {info.get('converged', False)}\n")
            f.write(f"HOMO Index: {info.get('homo_index', 'N/A')}\n")
            f.write(f"LUMO Index: {info.get('lumo_index', 'N/A')}\n")
            f.write(f"Occupied Orbitals: {info.get('num_occupied_orbitals', 'N/A')}\n")
            f.write(f"Virtual Orbitals: {info.get('num_virtual_orbitals', 'N/A')}\n")
            
            f.write(f"\nFiles in this directory:\n")
            f.write(f"- calculation.chk: PySCF checkpoint file with MO data\n")
            f.write(f"- calculation_info.txt: This summary file\n")
            f.write(f"- optimized_geometry.xyz: Optimized molecular structure\n")
    
    def save_geometry(self, calc_dir: str, xyz_string: str) -> None:
        """Save optimized geometry to XYZ file."""
        xyz_file = Path(calc_dir) / "optimized_geometry.xyz"
        with open(xyz_file, 'w') as f:
            f.write(xyz_string)
    
    def list_calculations(self) -> list:
        """List all calculation directories."""
        if not self.base_dir.exists():
            return []
        
        calculations = []
        for item in self.base_dir.iterdir():
            if item.is_dir():
                info_file = item / "calculation_info.txt"
                if info_file.exists():
                    calculations.append({
                        'name': item.name,
                        'path': str(item),
                        'date': datetime.fromtimestamp(item.stat().st_mtime),
                        'has_checkpoint': (item / "calculation.chk").exists()
                    })
        
        return sorted(calculations, key=lambda x: x['date'], reverse=True)
    
    def get_base_directory(self) -> str:
        """Get the base directory path."""
        return str(self.base_dir)
    
    def read_calculation_info(self, calc_dir: str) -> Optional[Dict[str, Any]]:
        """Read calculation info from the calculation directory."""
        info_file = Path(calc_dir) / "calculation_info.txt"
        
        if not info_file.exists():
            return None
        
        try:
            info = {}
            with open(info_file, 'r') as f:
                content = f.read()
                
            # Parse basic info from the file content
            lines = content.split('\n')
            for line in lines:
                if ':' in line and not line.startswith('=') and not line.startswith('-'):
                    key, value = line.split(':', 1)
                    key = key.strip().lower().replace(' ', '_')
                    value = value.strip()
                    info[key] = value
            
            return info
            
        except Exception as e:
            # If there's an error reading the file, return None
            return None
    
    def file_exists(self, calc_dir: str, filename: str) -> bool:
        """Check if a specific file exists in the calculation directory."""
        file_path = Path(calc_dir) / filename
        return file_path.exists()