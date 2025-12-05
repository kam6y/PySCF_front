"""
AI Agent module for PySCF_frontlication.
Provides multi-agent supervisor system for molecular analysis and research assistance.
"""

# The agent system uses a Supervisor pattern with specialized workers:
# - Quantum Calculation Worker: Handles quantum chemistry and molecular analysis
# - Research Agent: Handles academic literature search

# Main entry point
from .graph import get_compiled_graph

__all__ = ['get_compiled_graph']
