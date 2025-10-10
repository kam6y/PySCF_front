"""
Report Writer Worker

This module implements the Report Writer agent, which specializes in creating
comprehensive scientific reports by synthesizing information from other specialized agents.

The worker uses LangGraph's create_react_agent pattern with Google Gemini for
intelligent report generation without requiring external tools.
"""

import logging
from pathlib import Path
from typing import Optional
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.prebuilt import create_react_agent

# Set up logging
logger = logging.getLogger(__name__)


def get_gemini_api_key() -> Optional[str]:
    """
    Get Gemini API key from environment variable or settings file.

    Priority:
        1. Environment variable: GEMINI_API_KEY
        2. Settings file: settings.gemini_api_key
        3. None (fallback)

    Returns:
        Optional[str]: API key if found, None otherwise
    """
    import os
    from quantum_calc.settings_manager import get_current_settings

    # Priority 1: Environment variable
    api_key = os.environ.get("GEMINI_API_KEY")
    if api_key:
        logger.debug("Using Gemini API key from environment variable")
        return api_key

    # Priority 2: Settings file
    try:
        settings = get_current_settings()
        if settings.gemini_api_key:
            logger.debug("Using Gemini API key from settings file")
            return settings.gemini_api_key
    except Exception as e:
        logger.warning(f"Failed to load settings for API key: {e}")

    # Priority 3: None (fallback)
    return None


def _load_system_prompt() -> str:
    """Load the system prompt for the Report Writer."""
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
        "You are a Report Writer, an AI assistant specialized in creating comprehensive "
        "scientific reports for molecular science and computational chemistry. "
        "You synthesize information from other agents and produce well-structured, "
        "scientifically accurate reports in Japanese or English based on user preference."
    )


def create_report_writer():
    """
    Create a Report Writer agent using LangGraph's create_react_agent.

    This worker specializes in:
    - Creating comprehensive scientific reports
    - Synthesizing information from multiple sources
    - Generating calculation reports and literature reviews
    - Formatting professional documentation in Japanese or English

    The Report Writer does not use tools, relying entirely on the LLM's
    language generation capabilities to produce high-quality reports.

    Returns:
        A compiled LangGraph ReAct agent ready for use in supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Report Writer.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    # Initialize LLM with Gemini 2.5 Flash
    # Note: Consider using gemini-2.5-pro for higher quality reports if needed
    llm = ChatGoogleGenerativeAI(
        model="gemini-2.5-flash",
        api_key=api_key
    )
    logger.info("Initialized ChatGoogleGenerativeAI with gemini-2.5-flash for Report Writer")

    # Load system prompt
    system_prompt = _load_system_prompt()

    # Report Writer does not use tools - it purely generates text based on LLM capabilities
    tools = []
    logger.debug("Report Writer initialized with no tools (pure text generation)")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="report_writer",
        prompt=system_prompt
    )

    logger.info("Report Writer created successfully")
    return agent
