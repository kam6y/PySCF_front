"""
Deep Researcher Agent Module

This module contains the Deep Researcher implementation for academic paper search
and summarization. It specializes in querying arXiv and other scientific databases
to retrieve and summarize relevant research papers.

Components:
- tools.py: LangChain tools for academic search
- deep_researcher_agent.py: Deep Researcher implementation
"""

from .deep_researcher_agent import create_deep_researcher

__all__ = ["create_deep_researcher", "tools"]
