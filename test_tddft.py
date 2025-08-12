#!/usr/bin/env python3
"""
Test script for TDDFT implementation with NTO analysis
"""

import sys
import os
import logging

# Add the Python source directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src', 'python'))

from quantum_calc.tddft_calculator import TDDFTCalculator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Simple water molecule XYZ coordinates
water_xyz = """3
Water molecule
O    0.00000000    0.00000000    0.11779000
H    0.00000000    0.75545000   -0.47116000
H    0.00000000   -0.75545000   -0.47116000"""

def parse_xyz(xyz_string):
    """Parse XYZ string into atoms list"""
    lines = xyz_string.strip().split('\n')
    num_atoms = int(lines[0])
    atoms = []
    
    for i in range(2, 2 + num_atoms):  # Skip first two lines
        parts = lines[i].split()
        symbol = parts[0]
        coords = [float(x) for x in parts[1:4]]
        atoms.append([symbol, coords])
    
    return atoms

def test_tddft():
    """Test TDDFT calculation with NTO analysis"""
    print("Testing TDDFT implementation...")
    
    try:
        # Create calculator (let it use default temporary directory)
        calc = TDDFTCalculator(keep_files=True, molecule_name="water_test")
        
        # Parse water molecule
        atoms = parse_xyz(water_xyz)
        print(f"Parsed molecule with {len(atoms)} atoms")
        
        # Setup calculation parameters
        params = {
            'basis': '6-31G',  # Smaller basis for testing
            'xc': 'B3LYP',
            'charge': 0,
            'spin': 0,
            'nstates': 5,  # Request 5 excited states
            'tddft_method': 'TDDFT',
            'analyze_nto': True  # Enable NTO analysis
        }
        
        print("Setting up TDDFT calculation...")
        calc.setup_calculation(atoms, **params)
        
        print("Running TDDFT calculation...")
        results = calc.run_calculation()
        
        print("\n=== TDDFT Results ===")
        print(f"SCF Energy: {results.get('scf_energy', 'N/A')} hartree")
        print(f"Convergence: {results.get('converged', 'N/A')}")
        print(f"Number of excited states: {len(results.get('excitation_energies', []))}")
        
        if results.get('excitation_energies'):
            print("\n=== Excitation Energies ===")
            for i, energy in enumerate(results['excitation_energies'][:5]):
                wavelength = results['excitation_wavelengths'][i] if i < len(results.get('excitation_wavelengths', [])) else 'N/A'
                osc_strength = results['oscillator_strengths'][i] if i < len(results.get('oscillator_strengths', [])) else 'N/A'
                print(f"S{i+1}: {energy:.4f} eV, {wavelength:.1f} nm, f={osc_strength}")
        
        if results.get('nto_analysis'):
            print("\n=== NTO Analysis ===")
            for state_data in results['nto_analysis'][:3]:  # Show first 3 states
                print(f"State S{state_data['state']}: {len(state_data.get('nto_pairs', []))} NTO pairs")
                for pair in state_data.get('nto_pairs', [])[:3]:  # Show first 3 pairs
                    print(f"  {pair['hole_orbital']} -> {pair['particle_orbital']}: {pair['weight']:.6f} ({pair['contribution']:.1f}%)")
        
        print("\n=== Test completed successfully! ===")
        return True
        
    except Exception as e:
        print(f"\n=== Test failed with error ===")
        print(f"Error: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_tddft()
    sys.exit(0 if success else 1)