from __future__ import annotations

import logging
from typing import List

import numpy as np

logger = logging.getLogger(__name__)

class GeometryOptimizationMixin:
    """Mixin for geometry optimization and coordinate alignment workflows."""
    def _atoms_to_string(self, atoms: List[List]) -> str:
        """Convert atoms list to PySCF atom string format."""
        from .exceptions import GeometryError
        
        atom_lines = []
        for atom_data in atoms:
            symbol = atom_data[0]
            coords = atom_data[1]
            if len(coords) != 3:
                raise GeometryError(f"Invalid coordinates for atom {symbol}: expected 3 coordinates, got {len(coords)}")
            atom_lines.append(f"{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        return "\n".join(atom_lines)

    def _geometry_to_xyz_string(self) -> str:
        """Convert optimized geometry to XYZ format string."""
        if not hasattr(self, 'optimized_geometry') or self.optimized_geometry is None:
            return ""
        if not hasattr(self, 'mf') or self.mf is None or self.mf.mol is None:
            return ""

        # Use self.mf.mol to ensure consistency with optimized geometry
        mol = self.mf.mol
        atom_symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
        lines = [str(mol.natm)]
        lines.append("Optimized geometry from PySCF calculation")
        
        for i, (symbol, coords) in enumerate(zip(atom_symbols, self.optimized_geometry)):
            lines.append(f"{symbol:2s} {coords[0]:12.6f} {coords[1]:12.6f} {coords[2]:12.6f}")
        
        return "\n".join(lines)

    def _perform_geometry_optimization(self) -> None:
        """Perform geometry optimization using geometric_solver."""
        from pyscf.geomopt import geometric_solver

        logger.info("Starting geometry optimization...")

        # Integrate callback for geometry optimization
        # Note: geometric_solver.optimize accepts callback parameter
        try:
            optimized_mol = geometric_solver.optimize(self.mf, callback=self._geometry_optimization_callback)
        except TypeError:
            # Fallback if callback is not supported by this version of PySCF
            logger.warning("Geometry optimization callback not supported, using standard optimization")
            optimized_mol = geometric_solver.optimize(self.mf)

        self.optimized_geometry = optimized_mol.atom_coords(unit="ANG")
        logger.info("Geometry optimization completed")

        # Apply coordinate alignment after optimization
        self._align_optimized_geometry()

    def _align_optimized_geometry(self) -> None:
        """
        Align optimized geometry using ASE with correct rotation matrix approach:
        1. Move center of mass to origin
        2. Calculate principal axes using inertia tensor
        3. Create rotation matrix to align principal axes with coordinate axes
        4. Apply rotation matrix to all atomic coordinates
        """
        if not hasattr(self, 'optimized_geometry') or self.optimized_geometry is None:
            logger.warning("No optimized geometry available for alignment")
            return
        
        if not hasattr(self, 'mol') or self.mol is None:
            logger.warning("Molecular object not available for alignment")
            return
        
        try:
            # Import required libraries
            from ase import Atoms
            
            # Get atom symbols from the molecule
            atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
            positions = self.optimized_geometry.copy()
            
            logger.info("Starting molecular alignment with ASE...")
            logger.info(f"Original coordinates (first few atoms): {positions[:min(3, len(positions))]}")
            
            # Create ASE Atoms object
            atoms = Atoms(symbols=atom_symbols, positions=positions)
            
            # Move center of mass to origin
            atoms.center()
            centered_positions = atoms.get_positions()
            logger.info("Moved center of mass to origin")
            
            # Calculate moments of inertia and get principal axes
            inertia_moments, principal_axes = atoms.get_moments_of_inertia(vectors=True)
            logger.info(f"Original inertia moments: {inertia_moments}")
            logger.info(f"Principal axes shape: {principal_axes.shape}")
            
            # Check for special cases
            if len(atoms) == 1:
                logger.info("Single atom molecule - only centering applied")
                aligned_positions = centered_positions
                rotation_applied = False
            elif np.any(inertia_moments < 1e-10):
                logger.info("Linear molecule detected - limited alignment applied")
                # For linear molecules, align the molecular axis with Z-axis
                aligned_positions, rotation_applied = self._align_linear_molecule(centered_positions, principal_axes)
            else:
                # Apply rotation to align principal axes with coordinate axes
                logger.info("Applying rotation to align principal axes with coordinate axes")
                aligned_positions, rotation_applied = self._apply_rotation_matrix(centered_positions, principal_axes)
            
            # Update the optimized geometry with aligned coordinates
            self.optimized_geometry = aligned_positions
            
            # Verify alignment by checking the new inertia tensor
            atoms_aligned = Atoms(symbols=atom_symbols, positions=aligned_positions)
            new_inertia_moments = atoms_aligned.get_moments_of_inertia()
            
            if rotation_applied:
                # Calculate rotation angle for logging
                rotation_angle = self._calculate_rotation_angle(inertia_moments, new_inertia_moments)
                logger.info(f"Molecular alignment completed successfully (rotation angle: {rotation_angle:.2f}°)")
                logger.info(f"Aligned coordinates (first few atoms): {aligned_positions[:min(3, len(aligned_positions))]}")
                logger.info(f"Original inertia moments: {inertia_moments}")
                logger.info(f"Aligned inertia moments: {new_inertia_moments}")
                
                # Verify that rotation was meaningful
                if not self._validate_rotation(inertia_moments, new_inertia_moments):
                    logger.warning("Rotation validation failed - inertia moments unchanged")
            else:
                logger.info("No rotation applied - molecule already aligned or special case")
                logger.info(f"Final coordinates (first few atoms): {aligned_positions[:min(3, len(aligned_positions))]}")
            
        except ImportError as e:
            logger.error(f"ASE library not available for molecular alignment: {str(e)}")
            logger.info("Proceeding without molecular alignment")
        except Exception as e:
            logger.error(f"Molecular alignment failed: {str(e)}")
            logger.info("Proceeding with unaligned geometry")

    def _apply_rotation_matrix(self, positions: np.ndarray, principal_axes: np.ndarray) -> tuple[np.ndarray, bool]:
        """
        Apply rotation matrix to align principal axes with coordinate axes.
        
        Args:
            positions: Atomic positions (N x 3 array)
            principal_axes: Principal axes from inertia tensor (3 x 3 matrix)
            
        Returns:
            Tuple of (rotated_positions, rotation_applied_flag)
        """
        try:
            # The principal_axes matrix contains eigenvectors as columns
            # We want to rotate so that these axes align with [1,0,0], [0,1,0], [0,0,1]
            # This means we want: principal_axes @ rotation_matrix = identity
            # So: rotation_matrix = principal_axes^T
            rotation_matrix = principal_axes.T
            
            # Apply rotation to all atomic positions
            rotated_positions = positions @ rotation_matrix
            
            logger.info(f"Applied rotation matrix with determinant: {np.linalg.det(rotation_matrix):.6f}")
            
            return rotated_positions, True
            
        except Exception as e:
            logger.error(f"Failed to apply rotation matrix: {str(e)}")
            return positions, False

    def _align_linear_molecule(self, positions: np.ndarray, principal_axes: np.ndarray) -> tuple[np.ndarray, bool]:
        """
        Align linear molecule with Z-axis.
        
        Args:
            positions: Atomic positions (N x 3 array)  
            principal_axes: Principal axes from inertia tensor (3 x 3 matrix)
            
        Returns:
            Tuple of (aligned_positions, rotation_applied_flag)
        """
        try:
            # For linear molecules, find the axis with the largest moment of inertia
            # and align it with the Z-axis
            molecular_axis = principal_axes[:, -1]  # Last eigenvector (largest moment)
            z_axis = np.array([0, 0, 1])
            
            # Calculate rotation axis and angle
            rotation_axis = np.cross(molecular_axis, z_axis)
            rotation_angle = np.arccos(np.clip(np.dot(molecular_axis, z_axis), -1, 1))
            
            if rotation_angle < 1e-6:  # Already aligned
                logger.info("Linear molecule already aligned with Z-axis")
                return positions, False
            
            # Create rotation matrix using Rodrigues' formula
            rotation_axis_normalized = rotation_axis / np.linalg.norm(rotation_axis)
            cos_angle = np.cos(rotation_angle)
            sin_angle = np.sin(rotation_angle)
            
            # Rodrigues' rotation matrix formula
            K = np.array([[0, -rotation_axis_normalized[2], rotation_axis_normalized[1]],
                         [rotation_axis_normalized[2], 0, -rotation_axis_normalized[0]],
                         [-rotation_axis_normalized[1], rotation_axis_normalized[0], 0]])
            
            rotation_matrix = np.eye(3) + sin_angle * K + (1 - cos_angle) * K @ K
            
            # Apply rotation
            aligned_positions = positions @ rotation_matrix.T
            
            logger.info(f"Linear molecule aligned with Z-axis (rotation angle: {np.degrees(rotation_angle):.2f}°)")
            
            return aligned_positions, True
            
        except Exception as e:
            logger.error(f"Failed to align linear molecule: {str(e)}")
            return positions, False

    def _calculate_rotation_angle(self, original_moments: np.ndarray, aligned_moments: np.ndarray) -> float:
        """Calculate approximate rotation angle from inertia moment changes."""
        try:
            # Use the change in the smallest moment as an indicator
            original_sorted = np.sort(original_moments)
            aligned_sorted = np.sort(aligned_moments)
            
            # Calculate relative change
            relative_change = np.abs(original_sorted - aligned_sorted) / (original_sorted + 1e-10)
            max_change = np.max(relative_change)
            
            # Rough estimate: larger changes suggest larger rotations
            estimated_angle = np.degrees(np.arctan(max_change * 10))  # Heuristic scaling
            
            return min(estimated_angle, 180.0)  # Cap at 180 degrees
            
        except Exception:
            return 0.0

    def _validate_rotation(self, original_moments: np.ndarray, aligned_moments: np.ndarray) -> bool:
        """Validate that rotation was meaningful by checking inertia moment changes."""
        try:
            # Check if any moment changed significantly (more than 0.1%)
            relative_changes = np.abs(original_moments - aligned_moments) / (original_moments + 1e-10)
            max_change = np.max(relative_changes)
            
            # If max change is less than 0.001 (0.1%), rotation was likely ineffective
            if max_change < 0.001:
                return False
                
            # Check if moments are properly ordered (sorted)
            aligned_sorted = np.sort(aligned_moments)
            is_properly_ordered = np.allclose(aligned_moments, aligned_sorted, rtol=1e-6)
            
            return is_properly_ordered
            
        except Exception:
            return False

    def _create_optimized_molecule(self):
        """Create molecular object with optimized geometry."""
        from pyscf import gto
        
        # Get atom symbols from original molecule
        atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
        
        # Create atom string with optimized coordinates
        atom_lines = []
        for symbol, coords in zip(atom_symbols, self.optimized_geometry):
            atom_lines.append(f"{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        atom_string = "\n".join(atom_lines)
        
        # Create new molecular object
        optimized_mol = gto.M(
            atom=atom_string,
            basis=self.mol.basis,
            charge=self.mol.charge,
            spin=self.mol.spin,
            verbose=0
        )
        optimized_mol.max_memory = self.mol.max_memory
        
        return optimized_mol
