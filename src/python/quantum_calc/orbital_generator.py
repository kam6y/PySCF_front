"""Molecular orbital generator for visualization using PySCF CUBE files."""

import os
import tempfile
import logging
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import numpy as np
from pyscf import gto, lib, tools
from .exceptions import CalculationError, FileManagerError

logger = logging.getLogger(__name__)

# Hartree to eV conversion factor
HARTREE_TO_EV = 27.211386245988


class MolecularOrbitalGenerator:
    """Generate molecular orbital CUBE files and orbital information for visualization."""
    
    def __init__(self, working_dir: str):
        """
        Initialize the molecular orbital generator.
        
        Args:
            working_dir: Directory containing the calculation checkpoint file
        """
        self.working_dir = working_dir
        self.checkpoint_file = os.path.join(working_dir, "calculation.chk")
        self.mol: Optional[gto.Mole] = None
        self.mf = None
        self._orbital_info_cache: Optional[List[Dict[str, Any]]] = None
        
    def _load_calculation_data(self) -> None:
        """Load molecular and SCF data from checkpoint file."""
        if not os.path.exists(self.checkpoint_file):
            raise FileManagerError(f"Checkpoint file not found: {self.checkpoint_file}")
        
        try:
            # Load the molecule and SCF data from checkpoint file
            lib.chkfile.load(self.checkpoint_file, 'mol')
            self.mol = lib.chkfile.load_mol(self.checkpoint_file)
            
            # Try to load SCF data - works for DFT, HF, and MP2
            scf_data = lib.chkfile.load(self.checkpoint_file, 'scf')
            
            # Create appropriate mean-field object based on the calculation type
            from pyscf import dft, scf
            
            # Determine if this is a DFT or HF calculation
            xc = scf_data.get('xc', None)
            if xc and xc != 'HF':
                # DFT calculation
                if 'uhf' in scf_data and scf_data['uhf']:
                    self.mf = dft.UKS(self.mol)
                else:
                    self.mf = dft.RKS(self.mol)
                self.mf.xc = xc
            else:
                # HF calculation
                if 'uhf' in scf_data and scf_data['uhf']:
                    self.mf = scf.UHF(self.mol)
                else:
                    self.mf = scf.RHF(self.mol)
            
            # Load the coefficient matrix and other SCF data
            self.mf.mo_coeff = scf_data['mo_coeff']
            self.mf.mo_energy = scf_data['mo_energy']
            self.mf.mo_occ = scf_data['mo_occ']
            
            logger.info(f"Loaded calculation data from {self.checkpoint_file}")
            logger.info(f"Molecule: {self.mol.natm} atoms, {len(self.mf.mo_energy)} orbitals")
            
        except Exception as e:
            logger.error(f"Failed to load calculation data: {e}")
            raise CalculationError(f"Cannot load calculation data from checkpoint: {e}")
    
    def get_orbital_info(self) -> List[Dict[str, Any]]:
        """
        Get information about all molecular orbitals.
        
        Returns:
            List of orbital information dictionaries
        """
        if self._orbital_info_cache is not None:
            return self._orbital_info_cache
        
        if self.mol is None or self.mf is None:
            self._load_calculation_data()
        
        orbitals = []
        mo_energy = self.mf.mo_energy
        mo_occ = self.mf.mo_occ
        
        # Handle both RKS/RHF (1D array) and UKS/UHF (2D array) cases
        if hasattr(mo_energy, 'ndim') and mo_energy.ndim == 2:
            # UKS/UHF case: use alpha orbitals for display
            mo_energy = mo_energy[0]
            mo_occ = mo_occ[0]
        
        # Find HOMO and LUMO indices
        occupied_indices = np.where(mo_occ > 0)[0]
        homo_idx = occupied_indices[-1] if len(occupied_indices) > 0 else -1
        lumo_idx = np.where(mo_occ == 0)[0][0] if len(np.where(mo_occ == 0)[0]) > 0 else len(mo_energy)
        
        for i, (energy, occ) in enumerate(zip(mo_energy, mo_occ)):
            # Determine orbital type and label
            if i < homo_idx - 5:
                orbital_type = "core"
                label = f"Core {i}"
            elif i == homo_idx:
                orbital_type = "homo"
                label = "HOMO"
            elif i == lumo_idx:
                orbital_type = "lumo"
                label = "LUMO"
            elif i < homo_idx:
                diff = homo_idx - i
                orbital_type = "homo"
                label = f"HOMO-{diff}"
            elif i > lumo_idx:
                diff = i - lumo_idx
                orbital_type = "virtual"
                label = f"LUMO+{diff}"
            else:
                orbital_type = "virtual"
                label = f"Virtual {i}"
            
            orbital_info = {
                "index": int(i),
                "energy_hartree": float(energy),
                "energy_ev": float(energy * HARTREE_TO_EV),
                "occupancy": float(occ),
                "orbital_type": orbital_type,
                "label": label
            }
            orbitals.append(orbital_info)
        
        self._orbital_info_cache = orbitals
        return orbitals
    
    def get_orbital_summary(self) -> Dict[str, Any]:
        """
        Get summary information about orbitals.
        
        Returns:
            Dictionary with orbital summary
        """
        orbitals = self.get_orbital_info()
        
        # Find HOMO and LUMO
        homo_idx = None
        lumo_idx = None
        num_occupied = 0
        num_virtual = 0
        
        for orbital in orbitals:
            if orbital["orbital_type"] == "homo":
                homo_idx = orbital["index"]
            elif orbital["orbital_type"] == "lumo":
                lumo_idx = orbital["index"]
            
            if orbital["occupancy"] > 0:
                num_occupied += 1
            else:
                num_virtual += 1
        
        return {
            "orbitals": orbitals,
            "homo_index": homo_idx,
            "lumo_index": lumo_idx,
            "total_orbitals": len(orbitals),
            "num_occupied": num_occupied,
            "num_virtual": num_virtual
        }
    
    def generate_cube_file(
        self, 
        orbital_index: int,
        grid_size: int = 80,
        isovalue_pos: float = 0.02,
        isovalue_neg: float = -0.02,
        return_content: bool = True,
        save_to_disk: bool = False
    ) -> Dict[str, Any]:
        """
        Generate CUBE file for a specific molecular orbital.
        
        Args:
            orbital_index: Index of the orbital to generate
            grid_size: Grid size for CUBE file (default: 80)
            isovalue_pos: Positive isovalue for visualization
            isovalue_neg: Negative isovalue for visualization
            return_content: Whether to return file content as string
            save_to_disk: Whether to save file permanently to disk
        
        Returns:
            Dictionary with CUBE data and metadata
        """
        if self.mol is None or self.mf is None:
            self._load_calculation_data()
        
        orbitals = self.get_orbital_info()
        
        if orbital_index < 0 or orbital_index >= len(orbitals):
            raise CalculationError(f"Invalid orbital index: {orbital_index}. Available range: 0-{len(orbitals)-1}")
        
        orbital_info = orbitals[orbital_index]
        
        # Determine cube file path
        cube_filename = f"orbital_{orbital_index}_grid{grid_size}.cube"
        
        if save_to_disk:
            # Create orbital subdirectory and save there
            orbital_dir = os.path.join(self.working_dir, "orbital")
            os.makedirs(orbital_dir, exist_ok=True)
            cube_file_path = os.path.join(orbital_dir, cube_filename)
            
            # Check if file already exists and is valid
            if os.path.exists(cube_file_path):
                logger.info(f"Using existing CUBE file: {cube_file_path}")
                cube_content = ""
                file_size_kb = 0
                
                if return_content:
                    with open(cube_file_path, 'r') as f:
                        cube_content = f.read()
                    file_size_kb = os.path.getsize(cube_file_path) / 1024.0
                
                return {
                    "cube_data": cube_content,
                    "orbital_info": orbital_info,
                    "generation_params": {
                        "grid_size": grid_size,
                        "isovalue_positive": isovalue_pos,
                        "isovalue_negative": isovalue_neg,
                        "file_size_kb": file_size_kb
                    },
                    "file_path": cube_file_path,
                    "cached": True
                }
        else:
            # Create temporary file for CUBE output
            with tempfile.NamedTemporaryFile(mode='w', suffix='.cube', delete=False) as temp_file:
                cube_file_path = temp_file.name
        
        try:
            # Generate CUBE file using PySCF
            # Handle both RKS/RHF and UKS/UHF cases
            mo_coeff = self.mf.mo_coeff
            if hasattr(mo_coeff, 'ndim') and mo_coeff.ndim == 3:
                # UKS/UHF case: use alpha orbitals
                mo_coeff = mo_coeff[0]
            
            tools.cubegen.orbital(
                self.mol, 
                cube_file_path, 
                mo_coeff[:, orbital_index],
                nx=grid_size, 
                ny=grid_size, 
                nz=grid_size
            )
            
            # Read the generated CUBE file
            cube_content = ""
            file_size_kb = 0
            
            if return_content and os.path.exists(cube_file_path):
                with open(cube_file_path, 'r') as f:
                    cube_content = f.read()
                file_size_kb = os.path.getsize(cube_file_path) / 1024.0
            
            logger.info(f"Generated CUBE file for orbital {orbital_index} ({orbital_info['label']})")
            logger.info(f"Grid size: {grid_size}x{grid_size}x{grid_size}, File size: {file_size_kb:.1f} KB")
            
            return {
                "cube_data": cube_content,
                "orbital_info": orbital_info,
                "generation_params": {
                    "grid_size": grid_size,
                    "isovalue_positive": isovalue_pos,
                    "isovalue_negative": isovalue_neg,
                    "file_size_kb": file_size_kb
                },
                "file_path": cube_file_path if save_to_disk else None,
                "cached": False
            }
            
        except Exception as e:
            logger.error(f"Failed to generate CUBE file for orbital {orbital_index}: {e}")
            raise CalculationError(f"CUBE file generation failed: {e}")
        
        finally:
            # Clean up temporary file only if not saved to disk
            if not save_to_disk and os.path.exists(cube_file_path):
                try:
                    os.unlink(cube_file_path)
                except Exception as e:
                    logger.warning(f"Failed to clean up temporary CUBE file: {e}")
    
    def list_cube_files(self) -> List[Dict[str, Any]]:
        """
        List all CUBE files in the working directory.
        
        Returns:
            List of CUBE file information
        """
        cube_files = []
        
        orbital_dir = os.path.join(self.working_dir, "orbital")
        if not os.path.exists(orbital_dir):
            return cube_files
        
        # Find all cube files matching our naming pattern
        import glob
        pattern = os.path.join(orbital_dir, "orbital_*_grid*.cube")
        
        for file_path in glob.glob(pattern):
            filename = os.path.basename(file_path)
            
            # Extract orbital index and grid size from filename
            try:
                # Parse filename: orbital_{index}_grid{size}.cube
                parts = filename.replace('.cube', '').split('_')
                if len(parts) >= 3 and parts[0] == 'orbital' and parts[2].startswith('grid'):
                    orbital_index = int(parts[1])
                    grid_size = int(parts[2].replace('grid', ''))
                    
                    file_size_kb = os.path.getsize(file_path) / 1024.0
                    modified_time = datetime.fromtimestamp(os.path.getmtime(file_path))
                    
                    cube_files.append({
                        "filename": filename,
                        "file_path": file_path,
                        "orbital_index": orbital_index,
                        "grid_size": grid_size,
                        "file_size_kb": file_size_kb,
                        "modified": modified_time.isoformat()
                    })
            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse CUBE filename {filename}: {e}")
                continue
        
        return sorted(cube_files, key=lambda x: x["orbital_index"])
    
    def delete_cube_file(self, orbital_index: int, grid_size: int = None) -> bool:
        """
        Delete a specific CUBE file or all CUBE files for an orbital.
        
        Args:
            orbital_index: Index of the orbital
            grid_size: Specific grid size to delete, or None to delete all grid sizes
        
        Returns:
            True if file(s) were deleted successfully
        """
        deleted_count = 0
        orbital_dir = os.path.join(self.working_dir, "orbital")
        
        if not os.path.exists(orbital_dir):
            return False
        
        if grid_size is not None:
            # Delete specific file
            filename = f"orbital_{orbital_index}_grid{grid_size}.cube"
            file_path = os.path.join(orbital_dir, filename)
            
            if os.path.exists(file_path):
                try:
                    os.unlink(file_path)
                    deleted_count += 1
                    logger.info(f"Deleted CUBE file: {filename}")
                except Exception as e:
                    logger.error(f"Failed to delete CUBE file {filename}: {e}")
                    return False
        else:
            # Delete all grid sizes for this orbital
            import glob
            pattern = os.path.join(orbital_dir, f"orbital_{orbital_index}_grid*.cube")
            
            for file_path in glob.glob(pattern):
                try:
                    os.unlink(file_path)
                    deleted_count += 1
                    logger.info(f"Deleted CUBE file: {os.path.basename(file_path)}")
                except Exception as e:
                    logger.error(f"Failed to delete CUBE file {file_path}: {e}")
        
        return deleted_count > 0
    
    def delete_all_cube_files(self) -> int:
        """
        Delete all CUBE files in the working directory.
        
        Returns:
            Number of files deleted
        """
        deleted_count = 0
        orbital_dir = os.path.join(self.working_dir, "orbital")
        
        if not os.path.exists(orbital_dir):
            return deleted_count
        
        import glob
        pattern = os.path.join(orbital_dir, "orbital_*_grid*.cube")
        
        for file_path in glob.glob(pattern):
            try:
                os.unlink(file_path)
                deleted_count += 1
                logger.info(f"Deleted CUBE file: {os.path.basename(file_path)}")
            except Exception as e:
                logger.error(f"Failed to delete CUBE file {file_path}: {e}")
        
        logger.info(f"Deleted {deleted_count} CUBE files from {orbital_dir}")
        return deleted_count

    def validate_calculation(self) -> bool:
        """
        Validate that the calculation data is available for orbital generation.
        
        Returns:
            True if calculation data is valid, False otherwise
        """
        try:
            if self.mol is None or self.mf is None:
                self._load_calculation_data()
            
            # Check if we have the required data
            if (self.mf.mo_coeff is None or 
                self.mf.mo_energy is None or 
                self.mf.mo_occ is None):
                return False
            
            return True
            
        except Exception as e:
            logger.warning(f"Calculation validation failed: {e}")
            return False