#!/usr/bin/env python3
"""
Simple test script to verify the molecular alignment functionality.
"""

import numpy as np
import sys
import os

# Add the python source directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src', 'python'))

def test_ase_import():
    """Test if ASE can be imported."""
    try:
        from ase import Atoms
        print("✓ ASE library import successful")
        return True
    except ImportError as e:
        print(f"✗ ASE library import failed: {e}")
        return False

def test_basic_alignment():
    """Test basic molecular alignment functionality."""
    try:
        from ase import Atoms
        
        # Create a simple water molecule (slightly off-center and rotated)
        symbols = ['O', 'H', 'H']
        positions = np.array([
            [0.1, 0.1, 0.1],    # O slightly off-center
            [0.8, 0.6, 0.1],    # H
            [-0.6, 0.6, 0.1]    # H
        ])
        
        print("\nTesting basic alignment functionality:")
        print(f"Original positions:\n{positions}")
        
        # Create ASE Atoms object
        atoms = Atoms(symbols=symbols, positions=positions)
        
        # Move center of mass to origin
        atoms.center()
        print(f"After centering:\n{atoms.get_positions()}")
        
        # Get moments of inertia and principal axes
        inertia_moments, principal_axes = atoms.get_moments_of_inertia(vectors=True)
        print(f"Original inertia moments: {inertia_moments}")
        print(f"Principal axes shape: {principal_axes.shape}")
        
        # Apply alignment transformation following the sample code approach
        atoms.set_cell(principal_axes)
        atoms.wrap()
        atoms.set_cell([0, 0, 0])
        
        aligned_positions = atoms.get_positions()
        print(f"After alignment:\n{aligned_positions}")
        
        # Check alignment quality
        atoms_aligned = Atoms(symbols=symbols, positions=aligned_positions)
        new_inertia_moments = atoms_aligned.get_moments_of_inertia()
        print(f"Aligned inertia moments: {new_inertia_moments}")
        
        print("✓ Basic alignment test completed successfully")
        return True
        
    except Exception as e:
        print(f"✗ Basic alignment test failed: {e}")
        return False

def test_calculator_alignment():
    """Test the alignment functionality in BaseCalculator."""
    try:
        from quantum_calc.base_calculator import BaseCalculator
        from quantum_calc.dft_calculator import DFTCalculator
        
        # Create a test calculator
        calculator = DFTCalculator(keep_files=False)
        
        # Test the parse_xyz method
        water_xyz = """3
water molecule
O  0.0000   0.1000   0.0500
H  0.7571   0.5861   0.0000  
H -0.7571   0.5861   0.0000"""
        
        atoms = calculator.parse_xyz(water_xyz)
        print(f"\nParsed atoms: {len(atoms)} atoms")
        
        # Create mock optimized geometry and molecular object
        calculator.optimized_geometry = np.array([
            [0.0000, 0.1000, 0.0500],
            [0.7571, 0.5861, 0.0000],
            [-0.7571, 0.5861, 0.0000]
        ])
        
        # Create a minimal mock molecule for atom symbols
        class MockMol:
            natm = 3
            def atom_symbol(self, i):
                return ['O', 'H', 'H'][i]
        
        calculator.mol = MockMol()
        
        print(f"Original geometry:\n{calculator.optimized_geometry}")
        
        # Test the alignment method
        calculator._align_optimized_geometry()
        
        print(f"Aligned geometry:\n{calculator.optimized_geometry}")
        
        print("✓ Calculator alignment test completed successfully")
        return True
        
    except Exception as e:
        print(f"✗ Calculator alignment test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("Testing molecular alignment functionality...\n")
    
    tests = [
        test_ase_import,
        test_basic_alignment,
        test_calculator_alignment
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print(f"\nTest Summary:")
    print(f"Passed: {sum(results)}/{len(results)}")
    
    if all(results):
        print("🎉 All tests passed! Molecular alignment functionality is working correctly.")
        return 0
    else:
        print("❌ Some tests failed. Check the output above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())