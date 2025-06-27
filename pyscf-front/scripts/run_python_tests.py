#!/usr/bin/env python3
"""
Python test runner script
Runs all Python unit tests and integration tests
"""

import sys
import os
import unittest
import subprocess
from pathlib import Path

def main():
    """Run all Python tests"""
    # Get project root
    project_root = Path(__file__).parent.parent
    
    # Add src/python to Python path
    python_src = project_root / "src" / "python"
    sys.path.insert(0, str(python_src))
    
    # Discover and run tests
    test_dir = project_root / "tests"
    
    print("=" * 60)
    print("PySCF Front - Python Test Suite")
    print("=" * 60)
    
    # Unit tests
    print("\n🧪 Running Unit Tests...")
    print("-" * 40)
    
    unit_test_dir = test_dir / "unit"
    loader = unittest.TestLoader()
    
    # Discover unit tests
    unit_suite = unittest.TestSuite()
    
    for test_file in unit_test_dir.glob("test_*.py"):
        module_name = test_file.stem
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location(module_name, test_file)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            suite = loader.loadTestsFromModule(module)
            unit_suite.addTest(suite)
            print(f"✓ Loaded {module_name}")
        except Exception as e:
            print(f"✗ Failed to load {module_name}: {e}")
    
    # Run unit tests
    runner = unittest.TextTestRunner(verbosity=2)
    unit_result = runner.run(unit_suite)
    
    # Integration tests
    print("\n🔧 Running Integration Tests...")
    print("-" * 40)
    
    integration_test_dir = test_dir / "integration"
    integration_suite = unittest.TestSuite()
    
    for test_file in integration_test_dir.glob("test_*.py"):
        module_name = test_file.stem
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location(module_name, test_file)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            suite = loader.loadTestsFromModule(module)
            integration_suite.addTest(suite)
            print(f"✓ Loaded {module_name}")
        except Exception as e:
            print(f"✗ Failed to load {module_name}: {e}")
    
    # Run integration tests
    integration_result = runner.run(integration_suite)
    
    # Summary
    print("\n📊 Test Summary")
    print("-" * 40)
    
    total_tests = unit_result.testsRun + integration_result.testsRun
    total_failures = len(unit_result.failures) + len(integration_result.failures)
    total_errors = len(unit_result.errors) + len(integration_result.errors)
    
    print(f"Total Tests: {total_tests}")
    print(f"Passed: {total_tests - total_failures - total_errors}")
    print(f"Failed: {total_failures}")
    print(f"Errors: {total_errors}")
    
    if total_failures == 0 and total_errors == 0:
        print("\n🎉 All Python tests passed!")
        return 0
    else:
        print("\n❌ Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())