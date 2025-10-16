"""
Computational Chemist Agent Module

This module contains the Computational Chemist implementation
for quantum chemistry calculations.

Components:
- computational_chemist_agent.py: Computational Chemist implementation
- tools.py: Tool wrapper functions for quantum chemistry operations
"""

from .computational_chemist_agent import create_computational_chemist

__all__ = ["create_computational_chemist", "tools"]
