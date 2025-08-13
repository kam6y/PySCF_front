"""Base calculator class for quantum chemistry calculations."""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
import os
import tempfile


class BaseCalculator(ABC):
    """Abstract base class for quantum chemistry calculations."""
    
    def __init__(self, working_dir: Optional[str] = None):
        """Initialize calculator with optional working directory."""
        self.working_dir = working_dir or tempfile.mkdtemp(prefix="pyscf_calc_")
        self.results: Dict[str, Any] = {}
        
    def parse_xyz(self, xyz_string: str) -> List[List]:
        """Parse XYZ format string into atom list."""
        lines = xyz_string.strip().split('\n')
        if len(lines) < 3:
            raise ValueError("Invalid XYZ format: insufficient lines")
        
        try:
            atom_count = int(lines[0])
        except ValueError:
            raise ValueError("Invalid XYZ format: first line must be atom count")
        
        if len(lines) < atom_count + 2:
            raise ValueError(f"Invalid XYZ format: expected {atom_count + 2} lines, got {len(lines)}")
        
        atoms = []
        for i in range(2, atom_count + 2):
            parts = lines[i].split()
            if len(parts) < 4:
                raise ValueError(f"Invalid XYZ format at line {i + 1}: insufficient columns")
            
            symbol = parts[0]
            try:
                coords = [float(parts[j]) for j in range(1, 4)]
            except ValueError:
                raise ValueError(f"Invalid XYZ format at line {i + 1}: invalid coordinates")
            
            atoms.append([symbol, coords])
        
        return atoms
    
    def get_checkpoint_path(self) -> str:
        """Get the path to the checkpoint file."""
        return os.path.join(self.working_dir, "calculation.chk")
    
    def apply_resource_settings(self, mol, memory_mb: Optional[int] = None, cpu_cores: Optional[int] = None) -> None:
        """Apply resource settings to PySCF molecule object."""
        import logging
        logger = logging.getLogger(__name__)
        
        # メモリ設定
        if memory_mb is not None and memory_mb > 0:
            # PySCF expects memory in MB
            mol.max_memory = int(memory_mb)
            logger.info(f"Set PySCF max_memory to {memory_mb} MB")
        else:
            # デフォルトメモリ設定を適用
            mol.max_memory = 2000
            logger.info("Using default PySCF max_memory: 2000 MB")
        
        # CPU コア数設定 (PySCF lib.num_threads()を使用)
        if cpu_cores is not None and cpu_cores > 0:
            try:
                from pyscf import lib
                # 現在のスレッド数を取得
                current_threads = lib.num_threads()
                logger.info(f"Current PySCF threads: {current_threads}")
                
                # CPUコア数を設定
                lib.num_threads(int(cpu_cores))
                new_threads = lib.num_threads()
                logger.info(f"Set PySCF threads to: {new_threads} (requested: {cpu_cores})")
                
                if new_threads != int(cpu_cores):
                    if new_threads == 1 and int(cpu_cores) > 1:
                        logger.warning(f"PySCF parallelization is not available. This is likely because:")
                        logger.warning(f"1. PySCF was compiled without OpenMP support (common on macOS)")
                        logger.warning(f"2. OpenMP runtime is not available")
                        logger.warning(f"CPU cores requested: {cpu_cores}, but only 1 core will be used")
                        logger.warning(f"To enable multi-core support, reinstall PySCF with OpenMP:")
                        logger.warning(f"  pip uninstall pyscf")
                        logger.warning(f"  pip install pyscf[mkl] --upgrade")
                        logger.warning(f"  or install Intel MKL libraries manually")
                    else:
                        logger.warning(f"PySCF thread count mismatch: requested {cpu_cores}, got {new_threads}")
                else:
                    logger.info(f"✓ Successfully configured PySCF to use {new_threads} CPU cores")
                    
            except ImportError as e:
                logger.error(f"Failed to import PySCF lib module: {e}")
            except Exception as e:
                logger.error(f"Failed to set PySCF thread count: {e}")
        else:
            logger.info("No CPU core count specified, using PySCF default")
    
    @abstractmethod
    def setup_calculation(self, atoms: List[List], **kwargs) -> None:
        """Setup the quantum chemistry calculation."""
        pass
    
    @abstractmethod
    def run_calculation(self) -> Dict[str, Any]:
        """Run the quantum chemistry calculation and return results."""
        pass
    
    def cleanup(self, keep_files: bool = False) -> None:
        """Clean up temporary files."""
        import shutil
        if not keep_files and os.path.exists(self.working_dir) and "pyscf_calc_" in self.working_dir:
            shutil.rmtree(self.working_dir)
        elif keep_files:
            print(f"Calculation files preserved in: {self.working_dir}")