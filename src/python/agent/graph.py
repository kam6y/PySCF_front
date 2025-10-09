"""
LangGraph Multi-Agent Supervisor

This module implements the multi-agent architecture using LangGraph's Supervisor pattern.
It provides a stateful graph that coordinates specialized worker agents:
- Quantum Calculation Worker: Handles quantum chemistry and molecular analysis
- Research Agent: Handles academic paper search and literature reviews

The Supervisor intelligently delegates tasks to the appropriate worker based on
the user's intent and coordinates complex workflows.

Architecture:
- Supervisor: Central coordinator that analyzes requests and delegates to workers
- Workers: Specialized agents with domain-specific tools and expertise
- Hierarchical Control: Supervisor can call workers sequentially for complex tasks

This replaces the previous Router-based architecture with a more flexible and
powerful Supervisor pattern that enables complex multi-step workflows.
"""

import logging
from agent.supervisor import get_compiled_supervisor

# Set up logging
logger = logging.getLogger(__name__)


def get_compiled_graph():
    """
    Get a compiled instance of the multi-agent supervisor graph.

    This is the main entry point for the API to obtain the graph executor.
    The function delegates to the Supervisor module which handles all
    agent coordination and workflow management.

    Returns:
        Compiled LangGraph Supervisor application ready for streaming execution

    Raises:
        ValueError: If API key is not configured
        Exception: If agent initialization fails

    Example:
        >>> graph_app = get_compiled_graph()
        >>> result = graph_app.invoke({"messages": [{"role": "user", "content": "Run a DFT calculation"}]})
    """
    logger.info("Getting compiled Supervisor graph")
    return get_compiled_supervisor()
