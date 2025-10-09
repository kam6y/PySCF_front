"""
Quantum Calculation Worker

This module implements the Quantum Calculation Worker agent, which specializes in
quantum chemistry calculations, molecular analysis, and computational chemistry tasks.

The worker uses LangGraph's create_react_agent pattern with Google Gemini for
intelligent task execution and tool calling.
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
    """Load the system prompt for the Quantum Calculation Worker."""
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
        "You are a Quantum Calculation Worker, an AI assistant specialized in molecular analysis "
        "and quantum chemistry. You have access to tools for calculation management, molecular "
        "structure analysis, orbital visualization, and more. Use them proactively to help users "
        "with their computational chemistry tasks."
    )


def _initialize_tools():
    """Initialize and return the list of available tools for the Quantum Calculation Worker."""
    from agent import tools

    return [
        # Calculation management
        tools.list_all_calculations,
        tools.get_calculation_details,
        tools.start_quantum_calculation,
        tools.delete_calculation,  # HIL-enabled destructive tool

        # Molecular structure tools
        tools.search_pubchem_by_name,
        tools.convert_smiles_to_xyz,
        tools.validate_xyz_format,

        # Molecular orbital analysis
        tools.get_molecular_orbitals,
        tools.generate_orbital_cube,
        tools.list_cube_files,
        tools.delete_cube_files,

        # Spectroscopy
        tools.generate_ir_spectrum,

        # System and settings
        tools.get_supported_parameters,
        tools.get_app_settings,
        tools.update_app_settings,
        tools.get_system_resources,
    ]


def create_quantum_calculation_worker():
    """
    Create a Quantum Calculation Worker agent using LangGraph's create_react_agent.

    This worker specializes in:
    - Quantum chemistry calculations (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF)
    - Molecular structure analysis and visualization
    - Geometry optimization and frequency analysis
    - Molecular orbital and NTO analysis
    - IR spectrum generation

    Returns:
        A compiled LangGraph ReAct agent ready for use in supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Quantum Calculation Worker.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    # Initialize LLM with Gemini 2.5 Flash
    llm = ChatGoogleGenerativeAI(
        model="gemini-2.5-flash",
        api_key=api_key
    )
    logger.info("Initialized ChatGoogleGenerativeAI with gemini-2.5-flash")

    # Load system prompt
    system_prompt = _load_system_prompt()

    # Initialize tools
    tools = _initialize_tools()
    logger.debug(f"Initialized {len(tools)} tools for Quantum Calculation Worker")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="quantum_calculation_worker",
        prompt=system_prompt
    )

    logger.info("Quantum Calculation Worker created successfully")
    return agent
