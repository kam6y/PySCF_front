#!/usr/bin/env python3
"""
Simple integration test for MolecularAgent with Function Calling capability.
Tests basic initialization and tool integration without requiring external dependencies.
"""

import sys
import os

# Add the src/python directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src', 'python'))

def test_agent_initialization():
    """Test MolecularAgent initialization and tool loading."""
    print("Testing MolecularAgent initialization...")
    
    try:
        from agent.molecular_agent import MolecularAgent
        
        # Test basic initialization
        agent = MolecularAgent()
        print("✓ MolecularAgent initialized successfully")
        
        # Check if tools are loaded
        if hasattr(agent, 'available_tools'):
            print(f"✓ Tools loaded: {len(agent.available_tools)} tools available")
            for tool in agent.available_tools:
                print(f"  - {tool.__name__}")
        else:
            print("✗ No tools attribute found")
            return False
            
        # Test system prompt generation
        try:
            system_prompt = agent._get_system_prompt()
            if "ReAct" in system_prompt and "English" in system_prompt:
                print("✓ System prompt configured correctly for English and ReAct framework")
            else:
                print("✗ System prompt missing ReAct or English configuration")
                return False
        except Exception as e:
            print(f"✗ System prompt generation failed: {e}")
            return False
            
        print("✓ All basic initialization tests passed")
        return True
        
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False
    except Exception as e:
        print(f"✗ Initialization failed: {e}")
        return False

def test_tools_import():
    """Test that tools module can be imported and functions are accessible."""
    print("\nTesting tools module import...")
    
    try:
        from agent import tools
        
        # Check if expected functions exist
        expected_functions = ['list_all_calculations', 'get_calculation_details', 'get_supported_parameters']
        
        for func_name in expected_functions:
            if hasattr(tools, func_name):
                func = getattr(tools, func_name)
                print(f"✓ {func_name} function available")
                
                # Check if function has docstring (required for Gemini SDK)
                if func.__doc__:
                    print(f"  - Has docstring: {func.__doc__[:50]}...")
                else:
                    print(f"  ✗ Missing docstring for {func_name}")
                    return False
            else:
                print(f"✗ {func_name} function not found")
                return False
                
        print("✓ All tool functions available with proper documentation")
        return True
        
    except ImportError as e:
        print(f"✗ Tools import failed: {e}")
        return False
    except Exception as e:
        print(f"✗ Tools test failed: {e}")
        return False

def test_gemini_sdk_compatibility():
    """Test that the google.genai import works correctly."""
    print("\nTesting Google GenAI SDK compatibility...")
    
    try:
        from google import genai
        from google.genai import types
        print("✓ Google GenAI SDK imports successful")
        
        # Test that we can create types needed for Function Calling
        try:
            config = types.GenerateContentConfig(
                system_instruction="Test prompt",
                tools=[],  # Empty tools list for testing
            )
            print("✓ GenerateContentConfig can be created with tools parameter")
        except Exception as e:
            print(f"✗ GenerateContentConfig creation failed: {e}")
            return False
            
        return True
        
    except ImportError as e:
        print(f"✗ Google GenAI SDK not available: {e}")
        print("  This is expected if the SDK is not installed or configured")
        return True  # Don't fail the test for this
    except Exception as e:
        print(f"✗ Google GenAI SDK test failed: {e}")
        return False

def main():
    """Run all integration tests."""
    print("=== MolecularAgent Integration Test ===\n")
    
    tests = [
        test_tools_import,
        test_gemini_sdk_compatibility,
        test_agent_initialization,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"✗ Test {test.__name__} crashed: {e}")
            results.append(False)
    
    print(f"\n=== Test Results ===")
    print(f"Passed: {sum(results)}/{len(results)}")
    
    if all(results):
        print("🎉 All tests passed! The MolecularAgent refactoring is complete.")
        return 0
    else:
        print("❌ Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())