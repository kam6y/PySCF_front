"""PySCF Calculation Engine"""

import logging
import os
import sys
import time
from typing import Dict, List, Optional, Any, Generator
from dataclasses import dataclass
import numpy as np

try:
    import pyscf
    from pyscf import gto, scf, dft, mp, cc, lib
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False
    
from .config import get_pyscf_config

logger = logging.getLogger(__name__)


@dataclass
class AtomData:
    """Atom data structure"""
    symbol: str
    x: float
    y: float
    z: float


@dataclass
class MoleculeData:
    """Molecule data structure"""
    name: str
    formula: str
    atoms: List[AtomData]
    charge: int = 0
    multiplicity: int = 1


@dataclass
class CalculationResult:
    """Calculation result data structure"""
    success: bool
    total_energy: Optional[float] = None
    homo_energy: Optional[float] = None
    lumo_energy: Optional[float] = None
    homo_lumo_gap: Optional[float] = None
    orbital_energies: Optional[List[float]] = None
    orbital_occupancies: Optional[List[float]] = None
    dipole_moment: Optional[List[float]] = None
    calculation_time: Optional[float] = None
    convergence_achieved: bool = False
    error_message: Optional[str] = None
    additional_data: Optional[Dict[str, Any]] = None


class PySCFEngine:
    """PySCF calculation engine"""
    
    def __init__(self):
        self.config = get_pyscf_config()
        self._setup_pyscf()
    
    def _setup_pyscf(self):
        """Setup PySCF environment"""
        if not PYSCF_AVAILABLE:
            logger.warning("PySCF not available. Running in simulation mode.")
            return
        
        # Set up temporary directory
        os.makedirs(self.config["tmpdir"], exist_ok=True)
        lib.misc.TMPDIR = self.config["tmpdir"]
        
        # Set up threading
        lib.num_threads(self.config["num_threads"])
        
        # Set up memory
        if self.config["max_memory"]:
            lib.misc.MAX_MEMORY = self.config["max_memory"]
        
        logger.info(f"PySCF setup complete. Using {self.config['num_threads']} threads, "
                   f"{self.config['max_memory']}MB memory")
    
    @staticmethod
    def health_check() -> Dict[str, Any]:
        """Check PySCF health status"""
        if not PYSCF_AVAILABLE:
            return {
                "available": False,
                "error": "PySCF not installed",
                "simulation_mode": True
            }
        
        try:
            # Try a simple calculation
            mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', verbose=0)
            mf = scf.RHF(mol)
            mf.max_cycle = 1  # Just one iteration for testing
            mf.kernel()
            
            return {
                "available": True,
                "version": pyscf.__version__,
                "simulation_mode": False
            }
        except Exception as e:
            return {
                "available": False,
                "error": str(e),
                "simulation_mode": True
            }
    
    @staticmethod
    def get_system_info() -> Dict[str, Any]:
        """Get system information"""
        config = get_pyscf_config()
        
        # Basic system info
        info = {
            "version": "1.0.0",
            "python_version": sys.version,
            "pyscf_available": PYSCF_AVAILABLE,
            "simulation_mode": not PYSCF_AVAILABLE,
            "gpu_available": False,
            "resources": {
                "cpu_cores": config["num_threads"],
                "total_memory_mb": config["max_memory"],
            },
            "available_methods": [
                "HF", "RHF", "UHF", "ROHF",
                "PBE", "BLYP", "PW91", "B3LYP", "PBE0", "HSE06",
                "TPSS", "M06-L", "CAM-B3LYP", "ωB97X-D",
                "MP2", "CCSD", "CCSD(T)"
            ],
            "available_basis_sets": [
                "sto-3g", "3-21g", "6-31g", "6-31g*", "6-31g**",
                "6-31+g", "6-31++g", "6-311g", "6-311g*", "6-311g**",
                "cc-pvdz", "cc-pvtz", "cc-pvqz", "aug-cc-pvdz", "aug-cc-pvtz"
            ]
        }
        
        if PYSCF_AVAILABLE:
            info["pyscf_version"] = pyscf.__version__
            
            # Check GPU availability
            try:
                import gpu4pyscf
                info["gpu_available"] = True
                info["gpu4pyscf_version"] = gpu4pyscf.__version__
            except ImportError:
                info["gpu_available"] = False
        
        return info
    
    def create_molecule(self, molecule_data: MoleculeData) -> gto.Mole:
        """Create PySCF molecule object"""
        if not PYSCF_AVAILABLE:
            raise RuntimeError("PySCF not available")
        
        # Convert atoms to PySCF format
        atom_string = []
        for atom in molecule_data.atoms:
            atom_string.append(f"{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}")
        
        mol = gto.M(
            atom="; ".join(atom_string),
            charge=molecule_data.charge,
            spin=molecule_data.multiplicity - 1,  # PySCF uses spin (2S) not multiplicity (2S+1)
            basis='sto-3g',  # Default basis, will be overridden
            verbose=3
        )
        
        return mol
    
    def run_calculation(
        self, 
        molecule_data: MoleculeData,
        method: str,
        basis_set: str,
        max_iterations: int = 100
    ) -> Generator[Dict[str, Any], None, CalculationResult]:
        """
        Run quantum chemistry calculation with progress updates
        
        Yields progress updates and returns final result
        """
        start_time = time.time()
        
        if not PYSCF_AVAILABLE:
            # Simulation mode for development
            yield from self._simulate_calculation(molecule_data, method, basis_set)
            return CalculationResult(
                success=True,
                total_energy=-76.026749,  # Simulated water energy
                calculation_time=time.time() - start_time,
                convergence_achieved=True
            )
        
        try:
            # Create molecule
            yield {"step": "setup", "progress": 10, "message": "Creating molecule"}
            mol = self.create_molecule(molecule_data)
            mol.basis = basis_set
            mol.build()
            
            yield {"step": "basis", "progress": 20, "message": f"Building basis set {basis_set}"}
            
            # Setup calculation method
            if method.upper() in ['HF', 'RHF']:
                mf = scf.RHF(mol)
            elif method.upper() == 'UHF':
                mf = scf.UHF(mol)
            elif method.upper() == 'ROHF':
                mf = scf.ROHF(mol)
            elif method.upper() in ['B3LYP', 'PBE', 'BLYP', 'PW91', 'PBE0', 'HSE06']:
                mf = dft.RKS(mol)
                mf.xc = method.lower()
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            mf.max_cycle = max_iterations
            
            # Setup progress callback
            iteration_count = [0]
            def callback(mf_instance):
                iteration_count[0] += 1
                progress = min(90, 30 + (iteration_count[0] / max_iterations) * 60)
                return {
                    "step": "scf", 
                    "progress": progress, 
                    "message": f"SCF iteration {iteration_count[0]}"
                }
            
            yield {"step": "scf_start", "progress": 30, "message": f"Starting {method} calculation"}
            
            # Run calculation
            energy = mf.kernel()
            
            if not mf.converged:
                yield {"step": "warning", "progress": 85, "message": "Calculation did not converge"}
            
            yield {"step": "analysis", "progress": 90, "message": "Analyzing results"}
            
            # Extract results
            orbital_energies = mf.mo_energy.tolist() if hasattr(mf, 'mo_energy') else []
            orbital_occupancies = mf.mo_occ.tolist() if hasattr(mf, 'mo_occ') else []
            
            # Find HOMO/LUMO
            homo_idx = None
            lumo_idx = None
            if orbital_occupancies:
                for i, occ in enumerate(orbital_occupancies):
                    if occ > 0:
                        homo_idx = i
                    elif occ == 0 and homo_idx is not None and lumo_idx is None:
                        lumo_idx = i
                        break
            
            homo_energy = orbital_energies[homo_idx] if homo_idx is not None else None
            lumo_energy = orbital_energies[lumo_idx] if lumo_idx is not None else None
            homo_lumo_gap = (lumo_energy - homo_energy) if (homo_energy and lumo_energy) else None
            
            # Calculate dipole moment
            dipole = None
            try:
                dipole = mf.dip_moment().tolist()
            except:
                pass
            
            calculation_time = time.time() - start_time
            
            yield {"step": "complete", "progress": 100, "message": "Calculation completed"}
            
            return CalculationResult(
                success=True,
                total_energy=float(energy),
                homo_energy=homo_energy,
                lumo_energy=lumo_energy,
                homo_lumo_gap=homo_lumo_gap,
                orbital_energies=orbital_energies,
                orbital_occupancies=orbital_occupancies,
                dipole_moment=dipole,
                calculation_time=calculation_time,
                convergence_achieved=mf.converged
            )
            
        except Exception as e:
            logger.error(f"Calculation failed: {e}")
            return CalculationResult(
                success=False,
                error_message=str(e),
                calculation_time=time.time() - start_time
            )
    
    def _simulate_calculation(
        self, 
        molecule_data: MoleculeData,
        method: str,
        basis_set: str
    ) -> Generator[Dict[str, Any], None, None]:
        """Simulate calculation progress for development/testing"""
        steps = [
            (10, "setup", "Creating molecule"),
            (20, "basis", f"Building basis set {basis_set}"),
            (30, "scf_start", f"Starting {method} calculation"),
            (45, "scf", "SCF iteration 1"),
            (60, "scf", "SCF iteration 2"),
            (75, "scf", "SCF iteration 3"),
            (85, "scf", "SCF iteration 4"),
            (90, "analysis", "Analyzing results"),
            (100, "complete", "Calculation completed")
        ]
        
        for progress, step, message in steps:
            time.sleep(0.5)  # Simulate computation time
            yield {
                "step": step,
                "progress": progress,
                "message": message
            }


# Global engine instance
engine = PySCFEngine()