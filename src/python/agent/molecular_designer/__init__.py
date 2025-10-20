"""
Molecular Designer Agent Module

This module contains the Molecular Designer implementation
for generating novel molecular structures using cheminformatics tools.

Components:
- molecular_designer_agent.py: Molecular Designer implementation
- tools.py: RDKit-based tool wrapper functions for molecular generation
"""

from .molecular_designer_agent import create_molecular_designer

__all__ = ["create_molecular_designer", "tools"]
