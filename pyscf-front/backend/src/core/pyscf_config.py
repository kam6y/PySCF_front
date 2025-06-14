"""
PySCF Configuration Module
"""
import os
from typing import Optional


class PySCFConfig:
    """PySCF configuration management"""
    
    def __init__(self):
        self.setup_pyscf()
    
    def setup_pyscf(self):
        """Initialize PySCF configuration"""
        try:
            import pyscf
            from pyscf import lib
            
            # Memory and threading configuration
            num_threads = int(os.environ.get('OMP_NUM_THREADS', 4))
            lib.num_threads(num_threads)
            
            # Temporary directory
            pyscf.lib.misc.TMPDIR = os.environ.get('PYSCF_TMPDIR', '/tmp')
            
            # Convergence settings
            self.conv_tol = float(os.environ.get('PYSCF_CONV_TOL', 1e-9))
            self.conv_tol_grad = float(os.environ.get('PYSCF_CONV_TOL_GRAD', 1e-6))
            
            # GPU configuration
            self.use_gpu = os.environ.get('PYSCF_USE_GPU', 'false').lower() == 'true'
            self.gpu_available = self._check_gpu_availability()
            
            print(f"PySCF initialized with {num_threads} threads")
            if self.gpu_available:
                print("GPU acceleration available")
            
        except ImportError as e:
            print(f"Warning: PySCF not available: {e}")
            self.conv_tol = 1e-9
            self.conv_tol_grad = 1e-6
            self.use_gpu = False
            self.gpu_available = False
    
    def _check_gpu_availability(self) -> bool:
        """Check if GPU acceleration is available"""
        if not self.use_gpu:
            return False
        
        try:
            import gpu4pyscf
            return True
        except ImportError:
            print("Warning: GPU4PySCF not available")
            return False
    
    def get_default_memory(self) -> int:
        """Get default memory allocation in MB"""
        return int(os.environ.get('PYSCF_MAX_MEMORY', 4000))
    
    def get_scratch_dir(self) -> str:
        """Get scratch directory path"""
        import pyscf
        return pyscf.lib.misc.TMPDIR


# Global configuration instance
config = PySCFConfig()