"""
Literature Surveyor Agent Module

This module contains the Literature Surveyor implementation for academic paper search
and summarization. It specializes in querying arXiv and other scientific databases
to retrieve and summarize relevant research papers.

Components:
- tools.py: LangChain tools for academic search
- literature_surveyor_agent.py: Literature Surveyor implementation
"""

from .literature_surveyor_agent import create_literature_surveyor

__all__ = ["create_literature_surveyor", "tools"]
