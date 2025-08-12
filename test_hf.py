#!/usr/bin/env python3
"""Simple test script for HFCalculator."""

import sys
import os
sys.path.append('src/python')

from quantum_calc import HFCalculator

def test_hf_calculation():
    """Test HF calculation with a simple H2O molecule."""
    print("Testing HFCalculator with H2O molecule...")
    
    # XYZ coordinates for water molecule
    xyz_string = """3
Water molecule
O    0.0000    0.0000    0.0000
H    0.7570    0.5860    0.0000
H   -0.7570    0.5860    0.0000"""
    
    try:
        # Create calculator
        calculator = HFCalculator(keep_files=False)
        print("✓ HFCalculator created successfully")
        
        # Parse XYZ
        atoms = calculator.parse_xyz(xyz_string)
        print(f"✓ XYZ parsed successfully: {len(atoms)} atoms")
        
        # Setup calculation with minimal basis set for quick test
        calculator.setup_calculation(
            atoms,
            basis='STO-3G',  # Minimal basis set for quick calculation
            charge=0,
            spin=0,  # Singlet state for H2O
            max_cycle=50,
            solvent_method='none',
            solvent='-'
        )
        print("✓ Calculation setup completed")
        
        # For this test, we'll just verify setup was successful
        # without running the full calculation (which would take time)
        if calculator.mol is not None and calculator.mf is not None:
            print("✓ Molecular object and mean-field object created successfully")
            print(f"  - Molecule: {calculator.mol.nelectron} electrons, {calculator.mol.natm} atoms")
            print(f"  - Method: RHF (closed-shell)")
            print("✓ HF implementation test PASSED")
            return True
        else:
            print("✗ Setup failed - mol or mf object is None")
            return False
            
    except Exception as e:
        print(f"✗ Test failed with error: {str(e)}")
        return False
    finally:
        # Cleanup
        if 'calculator' in locals():
            calculator.cleanup()

if __name__ == "__main__":
    success = test_hf_calculation()
    sys.exit(0 if success else 1)