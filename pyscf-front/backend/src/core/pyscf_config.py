"""
PySCF Configuration and Setup
"""
import os
import logging
from typing import Dict, Any, Optional
import pyscf
from pyscf import gto, scf, dft, mp, cc, fci
import numpy as np

logger = logging.getLogger(__name__)

class PySCFConfig:
    """PySCF configuration manager"""
    
    def __init__(self):
        self.config = self._load_config()
        self._setup_pyscf()
    
    def _load_config(self) -> Dict[str, Any]:
        """Load PySCF configuration from environment"""
        return {
            'use_gpu': os.getenv('PYSCF_USE_GPU', 'false').lower() == 'true',
            'omp_threads': int(os.getenv('OMP_NUM_THREADS', '4')),
            'max_memory': int(os.getenv('PYSCF_MAX_MEMORY', '4000')),  # MB
            'scratch_dir': os.getenv('PYSCF_SCRATCH_DIR', '/tmp'),
            'verbose': int(os.getenv('PYSCF_VERBOSE', '3')),
            'convergence_threshold': float(os.getenv('PYSCF_CONV_TOL', '1e-8')),
            'max_cycles': int(os.getenv('PYSCF_MAX_CYCLES', '50'))
        }
    
    def _setup_pyscf(self):
        """Setup PySCF global configuration"""
        # Set number of threads
        pyscf.lib.num_threads(self.config['omp_threads'])
        
        # Set memory limit
        pyscf.lib.param.MAX_MEMORY = self.config['max_memory']
        
        # Setup GPU if available
        if self.config['use_gpu']:
            try:
                import gpu4pyscf
                logger.info("GPU acceleration enabled")
            except ImportError:
                logger.warning("GPU requested but gpu4pyscf not available")
                self.config['use_gpu'] = False
        
        logger.info(f"PySCF configured with {self.config['omp_threads']} threads, "
                   f"{self.config['max_memory']} MB memory")
    
    def get_available_methods(self) -> list:
        """Get list of available quantum chemistry methods"""
        methods = [
            'HF',     # Hartree-Fock
            'DFT',    # Density Functional Theory
            'MP2',    # Møller-Plesset 2nd order
            'CCSD',   # Coupled Cluster Singles Doubles
            'CASSCF', # Complete Active Space SCF
            'FCI'     # Full Configuration Interaction
        ]
        
        # Check for GPU-accelerated methods
        if self.config['use_gpu']:
            methods.extend(['HF-GPU', 'DFT-GPU'])
        
        return methods
    
    def get_available_basis_sets(self) -> list:
        """Get list of available basis sets"""
        return [
            'sto-3g',
            '3-21g',
            '6-31g',
            '6-31g*',
            '6-31g**',
            '6-311g',
            '6-311g*',
            '6-311g**',
            'cc-pvdz',
            'cc-pvtz',
            'cc-pvqz',
            'aug-cc-pvdz',
            'aug-cc-pvtz',
            'def2-sv',
            'def2-svp',
            'def2-tzvp'
        ]
    
    def get_available_functionals(self) -> list:
        """Get list of available DFT functionals"""
        return [
            'lda',
            'pbe',
            'b3lyp',
            'b3lyp5',
            'pbe0',
            'wb97x',
            'm06',
            'm06-2x',
            'tpss',
            'scan'
        ]

# Global configuration instance
pyscf_config = PySCFConfig()


class MoleculeBuilder:
    """Helper class for building PySCF molecule objects"""
    
    @staticmethod
    def from_atoms(atoms: list, charge: int = 0, spin: int = 0, 
                   basis: str = 'sto-3g', symmetry: str = 'c1') -> gto.Mole:
        """
        Build PySCF molecule from atom list
        
        Args:
            atoms: List of (symbol, x, y, z) tuples
            charge: Total charge
            spin: 2*S (multiplicity - 1)
            basis: Basis set name
            symmetry: Point group symmetry
        """
        mol = gto.Mole()
        
        # Build atom string
        atom_string = []
        for atom in atoms:
            if len(atom) == 4:  # (symbol, x, y, z)
                symbol, x, y, z = atom
                atom_string.append(f"{symbol} {x} {y} {z}")
            else:
                raise ValueError(f"Invalid atom format: {atom}")
        
        mol.atom = '; '.join(atom_string)
        mol.charge = charge
        mol.spin = spin
        mol.basis = basis
        mol.symmetry = symmetry
        mol.verbose = pyscf_config.config['verbose']
        
        try:
            mol.build()
            logger.info(f"Built molecule: {mol.formula} with {mol.natm} atoms, "
                       f"charge={charge}, spin={spin}")
            return mol
        except Exception as e:
            logger.error(f"Failed to build molecule: {e}")
            raise
    
    @staticmethod
    def from_xyz_file(xyz_path: str, charge: int = 0, spin: int = 0,
                     basis: str = 'sto-3g') -> gto.Mole:
        """Build molecule from XYZ file"""
        mol = gto.Mole()
        mol.fromfile(xyz_path)
        mol.charge = charge
        mol.spin = spin
        mol.basis = basis
        mol.verbose = pyscf_config.config['verbose']
        mol.build()
        return mol
    
    @staticmethod
    def validate_geometry(atoms: list) -> tuple[bool, str]:
        """
        Validate molecular geometry
        
        Returns:
            (is_valid, error_message)
        """
        if not atoms:
            return False, "No atoms provided"
        
        if len(atoms) < 1:
            return False, "At least one atom required"
        
        # Check atom format
        for i, atom in enumerate(atoms):
            if len(atom) != 4:
                return False, f"Atom {i}: Expected (symbol, x, y, z), got {atom}"
            
            symbol, x, y, z = atom
            if not isinstance(symbol, str):
                return False, f"Atom {i}: Symbol must be string, got {type(symbol)}"
            
            try:
                float(x), float(y), float(z)
            except (ValueError, TypeError):
                return False, f"Atom {i}: Coordinates must be numeric"
        
        # Check for duplicate atoms (too close)
        coords = np.array([[float(atom[1]), float(atom[2]), float(atom[3])] 
                          for atom in atoms])
        
        if len(coords) > 1:
            distances = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
            np.fill_diagonal(distances, np.inf)  # Ignore self-distances
            
            min_distance = np.min(distances)
            if min_distance < 0.1:  # Less than 0.1 Angstrom
                return False, f"Atoms too close: minimum distance {min_distance:.3f} Å"
        
        return True, "Geometry is valid"


class CalculationEngine:
    """Main calculation engine for quantum chemistry methods"""
    
    def __init__(self):
        self.config = pyscf_config
    
    def run_hf_calculation(self, mol: gto.Mole, **kwargs) -> dict:
        """Run Hartree-Fock calculation"""
        try:
            # Create HF object
            if mol.spin == 0:
                mf = scf.RHF(mol)
            else:
                mf = scf.UHF(mol)
            
            # Set convergence parameters
            mf.conv_tol = kwargs.get('conv_tol', self.config.config['convergence_threshold'])
            mf.max_cycle = kwargs.get('max_cycle', self.config.config['max_cycles'])
            
            # Run calculation
            logger.info(f"Starting HF calculation for {mol.formula}")
            energy = mf.kernel()
            
            if not mf.converged:
                logger.warning("HF calculation did not converge")
            
            # Extract results
            results = {
                'converged': mf.converged,
                'total_energy': float(energy),
                'nuclear_repulsion': float(mol.energy_nuc()),
                'electronic_energy': float(energy - mol.energy_nuc()),
                'homo_energy': float(mf.mo_energy[mf.mo_occ > 0][-1]) if hasattr(mf, 'mo_energy') else None,
                'lumo_energy': float(mf.mo_energy[mf.mo_occ == 0][0]) if hasattr(mf, 'mo_energy') and np.any(mf.mo_occ == 0) else None,
                'orbital_energies': mf.mo_energy.tolist() if hasattr(mf, 'mo_energy') else [],
                'occupations': mf.mo_occ.tolist() if hasattr(mf, 'mo_occ') else [],
                'cycles': mf.cycle if hasattr(mf, 'cycle') else 0
            }
            
            logger.info(f"HF calculation completed: E = {energy:.6f} Eh")
            return results
            
        except Exception as e:
            logger.error(f"HF calculation failed: {e}")
            raise
    
    def run_dft_calculation(self, mol: gto.Mole, functional: str = 'b3lyp', **kwargs) -> dict:
        """Run DFT calculation"""
        try:
            # Create DFT object
            if mol.spin == 0:
                mf = dft.RKS(mol)
            else:
                mf = dft.UKS(mol)
            
            # Set functional
            mf.xc = functional
            
            # Set convergence parameters
            mf.conv_tol = kwargs.get('conv_tol', self.config.config['convergence_threshold'])
            mf.max_cycle = kwargs.get('max_cycle', self.config.config['max_cycles'])
            
            # Run calculation
            logger.info(f"Starting DFT({functional}) calculation for {mol.formula}")
            energy = mf.kernel()
            
            if not mf.converged:
                logger.warning("DFT calculation did not converge")
            
            # Extract results
            results = {
                'converged': mf.converged,
                'total_energy': float(energy),
                'nuclear_repulsion': float(mol.energy_nuc()),
                'electronic_energy': float(energy - mol.energy_nuc()),
                'functional': functional,
                'homo_energy': float(mf.mo_energy[mf.mo_occ > 0][-1]) if hasattr(mf, 'mo_energy') else None,
                'lumo_energy': float(mf.mo_energy[mf.mo_occ == 0][0]) if hasattr(mf, 'mo_energy') and np.any(mf.mo_occ == 0) else None,
                'orbital_energies': mf.mo_energy.tolist() if hasattr(mf, 'mo_energy') else [],
                'occupations': mf.mo_occ.tolist() if hasattr(mf, 'mo_occ') else [],
                'cycles': mf.cycle if hasattr(mf, 'cycle') else 0
            }
            
            logger.info(f"DFT calculation completed: E = {energy:.6f} Eh")
            return results
            
        except Exception as e:
            logger.error(f"DFT calculation failed: {e}")
            raise
    
    def run_mp2_calculation(self, mol: gto.Mole, **kwargs) -> dict:
        """Run MP2 calculation"""
        try:
            # First run HF
            if mol.spin == 0:
                mf = scf.RHF(mol)
            else:
                mf = scf.UHF(mol)
            
            mf.conv_tol = kwargs.get('conv_tol', self.config.config['convergence_threshold'])
            mf.kernel()
            
            if not mf.converged:
                raise RuntimeError("HF reference calculation did not converge")
            
            # Run MP2
            logger.info(f"Starting MP2 calculation for {mol.formula}")
            mp2_solver = mp.MP2(mf)
            e_corr, t2 = mp2_solver.kernel()
            
            total_energy = mf.e_tot + e_corr
            
            results = {
                'converged': True,
                'total_energy': float(total_energy),
                'hf_energy': float(mf.e_tot),
                'correlation_energy': float(e_corr),
                'nuclear_repulsion': float(mol.energy_nuc()),
                'electronic_energy': float(total_energy - mol.energy_nuc())
            }
            
            logger.info(f"MP2 calculation completed: E = {total_energy:.6f} Eh")
            return results
            
        except Exception as e:
            logger.error(f"MP2 calculation failed: {e}")
            raise


# Global calculation engine instance
calculation_engine = CalculationEngine()