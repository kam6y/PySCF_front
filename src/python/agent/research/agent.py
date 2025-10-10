"""
Research Agent Module

This module implements the Research Agent, which specializes in searching and
summarizing academic papers from sources like arXiv.

The Research Agent uses LangGraph's create_react_agent pattern with Google Gemini
for intelligent literature search and analysis.
"""

import logging
from pathlib import Path
from flask import current_app
from langchain_google_genai import ChatGoogleGenerativeAI

from .tools import search_arxiv
from agent.quantum_calc.quantum_calc_worker import get_gemini_api_key

# Set up logging
logger = logging.getLogger(__name__)


def _load_system_prompt() -> str:
    """Load the system prompt for the Research Expert."""
    try:
        # Try to load system prompt from file
        prompt_path = Path(__file__).parent / "prompts" / "system_prompt.txt"

        if prompt_path.exists():
            with open(prompt_path, 'r', encoding='utf-8') as f:
                prompt = f.read()
            logger.debug(f"Loaded system prompt from {prompt_path}")
            return prompt
        else:
            logger.warning(f"System prompt file not found at {prompt_path}, using fallback")
            return _get_fallback_system_prompt()

    except Exception as e:
        logger.error(f"Failed to load system prompt from file: {e}, using fallback")
        return _get_fallback_system_prompt()


def _get_fallback_system_prompt() -> str:
    """Fallback system prompt in case file loading fails."""
    return (
        "You are a Research Expert specializing in academic literature search for "
        "quantum chemistry and computational chemistry. Use the search_arxiv tool to find "
        "relevant papers and provide detailed explanations of their findings. Always explain "
        "papers in detail, not just listing titles."
    )


def create_research_agent_runnable():
    """
    Create a Research Expert agent using LangGraph's create_react_agent.

    This worker specializes in:
    - Academic literature search (arXiv)
    - Paper analysis and summarization
    - Identifying trends and developments in quantum chemistry research
    - Providing formatted citations with PDF links

    The agent uses LangGraph's create_react_agent pattern to automatically handle
    tool execution loops, allowing the LLM to call search_arxiv and process results.

    Returns:
        A compiled LangGraph ReAct agent ready for use in supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    from langgraph.prebuilt import create_react_agent

    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Research Expert.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    # Initialize LLM with configured model
    model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
    llm = ChatGoogleGenerativeAI(
        model=model_name,
        api_key=api_key
    )
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Research Expert")

    # Load system prompt
    system_prompt = _load_system_prompt()

    # Create ReAct agent with arXiv search tool
    # The agent will automatically:
    # 1. Decide when to use tools
    # 2. Execute tool calls
    # 3. Process tool outputs
    # 4. Generate final response
    agent = create_react_agent(
        model=llm,
        tools=[search_arxiv],
        name="research_expert",
        prompt=system_prompt
    )

    logger.info("Research Expert created successfully")
    return agent
