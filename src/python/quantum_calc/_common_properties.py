from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

class CommonPropertiesMixin:
    """Mixin for orbital, Mulliken, and shared property extraction helpers."""
    def _analyze_orbitals(self) -> Tuple[int, int]:
        """Analyze molecular orbitals to find HOMO and LUMO indices."""
        from .exceptions import CalculationError
        
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mo_occ is None:
            raise CalculationError("Orbital occupations not available")
        
        # Handle both RKS/RHF (1D array) and UKS/UHF (2D array) cases
        mo_occ = self._as_numpy_array(self.mf.mo_occ)
        if hasattr(mo_occ, 'ndim') and mo_occ.ndim == 2:
            # UKS/UHF case: use alpha orbitals
            mo_occ = mo_occ[0]
        
        # Find HOMO (highest occupied molecular orbital)
        occupied_indices = np.where(mo_occ > 0)[0]
        if len(occupied_indices) == 0:
            raise CalculationError("No occupied orbitals found")
        homo_idx = occupied_indices[-1]
        
        # Find LUMO (lowest unoccupied molecular orbital)
        unoccupied_indices = np.where(mo_occ == 0)[0]
        if len(unoccupied_indices) == 0:
            raise CalculationError("No unoccupied orbitals found")
        lumo_idx = unoccupied_indices[0]
        
        return int(homo_idx), int(lumo_idx)

    def _count_occupied_orbitals(self) -> int:
        """Count the number of occupied orbitals."""
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mo_occ is None:
            return 0
        
        mo_occ = self._as_numpy_array(self.mf.mo_occ)
        return int(np.sum(mo_occ > 0))

    def _count_virtual_orbitals(self) -> int:
        """Count the number of virtual orbitals."""
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mo_occ is None:
            return 0
        
        mo_occ = self._as_numpy_array(self.mf.mo_occ)
        return int(np.sum(mo_occ == 0))

    def _calculate_mulliken_charges(self) -> Optional[List[Dict[str, Any]]]:
        """Calculate Mulliken population analysis charges for each atom."""
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mol is None:
            return None

        # Use self.mf.mol to ensure consistency with the converged calculation
        mol = self.mf.mol

        try:
            # Perform Mulliken population analysis
            # This returns (pop, charges) where pop are populations and charges are atomic charges
            pop, charges = self.mf.mulliken_pop()
            charges = self._as_numpy_array(charges)

            # Extract charges for each atom
            mulliken_charges = []
            for i in range(mol.natm):
                atom_symbol = mol.atom_symbol(i)
                # Convert numpy float to Python float for JSON serialization
                charge = float(charges[i])

                mulliken_charges.append({
                    'atom_index': i,
                    'element': atom_symbol,
                    'charge': charge
                })

            # Verify total charge conservation (should equal molecular charge)
            total_charge = sum(item['charge'] for item in mulliken_charges)
            expected_charge = float(mol.charge)
            logger.info(f"Mulliken analysis: calculated total charge = {total_charge:.6f}, expected = {expected_charge:.6f}")
            
            if abs(total_charge - expected_charge) > 0.01:
                logger.warning(f"Mulliken total charge ({total_charge:.6f}) differs from molecular charge ({expected_charge:.6f}) by more than 0.01")
            
            return mulliken_charges
            
        except Exception as e:
            logger.warning(f"Mulliken population analysis failed: {str(e)}")
            return None

    def _perform_orbital_analysis(self) -> None:
        """Perform orbital analysis and store results."""
        logger.info("Performing orbital analysis...")
        homo_idx, lumo_idx = self._analyze_orbitals()
        
        self.results.update({
            'homo_index': homo_idx,
            'lumo_index': lumo_idx,
            'num_occupied_orbitals': int(self._count_occupied_orbitals()),
            'num_virtual_orbitals': int(self._count_virtual_orbitals())
        })
        logger.info("Orbital analysis completed")

    def _perform_mulliken_analysis(self) -> None:
        """Perform Mulliken population analysis."""
        logger.info("Performing Mulliken population analysis...")
        mulliken_charges = self._calculate_mulliken_charges()
        
        self.results['mulliken_charges'] = mulliken_charges
        
        if mulliken_charges is not None:
            logger.info(f"Calculated Mulliken charges for {len(mulliken_charges)} atoms")
        else:
            logger.warning("Mulliken population analysis failed or was skipped")

    def _extract_common_additional_properties(self) -> Dict[str, Any]:
        """
        Extract common additional properties available for all calculation methods.

        This method provides comprehensive electronic structure information including:
        - Dipole moment (x, y, z components and total, in Debye and a.u.)
        - HOMO-LUMO energies and gap (in hartree and eV)
        - Energy components (nuclear repulsion, electronic energy)
        - Basis set information (number of basis functions, primitive Gaussians)

        Returns:
            Dictionary containing all available common properties
        """
        properties = {}

        if self.mf is None or self.mf.mol is None:
            logger.warning("Mean field or molecular object not available for additional properties extraction")
            return properties

        # Use self.mf.mol to ensure consistency with the converged calculation
        mol = self.mf.mol

        try:
            # 1. Dipole Moment (双極子モーメント)
            try:
                logger.info("Calculating dipole moment...")
                # Get dipole moment in both Debye and atomic units
                dip_debye = np.asarray(self._to_numpy(self.mf.dip_moment(unit='Debye')))
                dip_au = np.asarray(self._to_numpy(self.mf.dip_moment(unit='A.U.')))

                properties['dipole_moment_x_debye'] = float(dip_debye[0])
                properties['dipole_moment_y_debye'] = float(dip_debye[1])
                properties['dipole_moment_z_debye'] = float(dip_debye[2])
                properties['dipole_moment_total_debye'] = float(np.linalg.norm(dip_debye))
                properties['dipole_moment_x_au'] = float(dip_au[0])
                properties['dipole_moment_y_au'] = float(dip_au[1])
                properties['dipole_moment_z_au'] = float(dip_au[2])
                properties['dipole_moment_total_au'] = float(np.linalg.norm(dip_au))

                logger.info(f"Dipole moment: {properties['dipole_moment_total_debye']:.4f} Debye")
            except Exception as e:
                logger.warning(f"Failed to calculate dipole moment: {e}")

            # 2. HOMO-LUMO Gap and Individual Energies (HOMO-LUMOギャップと個別エネルギー)
            try:
                logger.info("Calculating HOMO-LUMO energies and gap...")
                mo_energy = self._as_numpy_array(self.mf.mo_energy)

                # Handle both RKS/RHF (1D array) and UKS/UHF (2D array) cases
                if hasattr(mo_energy, 'ndim') and mo_energy.ndim == 2:
                    # UKS/UHF case: use alpha orbitals
                    mo_energy = mo_energy[0]

                homo_idx = self.results.get('homo_index')
                lumo_idx = self.results.get('lumo_index')

                if homo_idx is not None and lumo_idx is not None:
                    homo_energy_hartree = float(mo_energy[homo_idx])
                    lumo_energy_hartree = float(mo_energy[lumo_idx])
                    gap_hartree = lumo_energy_hartree - homo_energy_hartree

                    # Convert to eV (1 hartree = 27.2114 eV)
                    HARTREE_TO_EV = 27.2114
                    properties['homo_energy_hartree'] = homo_energy_hartree
                    properties['homo_energy_ev'] = homo_energy_hartree * HARTREE_TO_EV
                    properties['lumo_energy_hartree'] = lumo_energy_hartree
                    properties['lumo_energy_ev'] = lumo_energy_hartree * HARTREE_TO_EV
                    properties['homo_lumo_gap_hartree'] = gap_hartree
                    properties['homo_lumo_gap_ev'] = gap_hartree * HARTREE_TO_EV

                    logger.info(f"HOMO energy: {properties['homo_energy_ev']:.4f} eV")
                    logger.info(f"LUMO energy: {properties['lumo_energy_ev']:.4f} eV")
                    logger.info(f"HOMO-LUMO gap: {properties['homo_lumo_gap_ev']:.4f} eV")
            except Exception as e:
                logger.warning(f"Failed to calculate HOMO-LUMO energies: {e}")

            # 3. Energy Components (エネルギー成分)
            try:
                logger.info("Calculating energy components...")
                properties['nuclear_repulsion_energy'] = float(self.mf.energy_nuc())
                e_tot = self._to_numpy(self.mf.e_tot)
                e_nuc = properties['nuclear_repulsion_energy']
                properties['electronic_energy'] = float(e_tot - e_nuc)

                logger.info(f"Nuclear repulsion energy: {properties['nuclear_repulsion_energy']:.6f} hartree")
                logger.info(f"Electronic energy: {properties['electronic_energy']:.6f} hartree")
            except Exception as e:
                logger.warning(f"Failed to calculate energy components: {e}")

            # 4. Basis Set Information (基底関数情報)
            try:
                logger.info("Extracting basis set information...")
                properties['num_basis_functions'] = int(mol.nao)
                properties['num_primitive_gaussians'] = int(mol.npgto_nr())
                properties['total_electrons'] = int(mol.nelectron)

                logger.info(f"Number of basis functions: {properties['num_basis_functions']}")
                logger.info(f"Number of primitive Gaussians: {properties['num_primitive_gaussians']}")
                logger.info(f"Total electrons: {properties['total_electrons']}")
            except Exception as e:
                logger.warning(f"Failed to extract basis set information: {e}")

            logger.info("Common additional properties extraction completed")

        except Exception as e:
            logger.error(f"Failed to extract common additional properties: {e}")

        return properties
