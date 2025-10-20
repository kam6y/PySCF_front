"""
Supervisor Agent

This module implements the Supervisor agent that coordinates multiple specialized worker agents
using LangGraph's Supervisor pattern.

The Supervisor analyzes user requests, delegates tasks to appropriate workers
(Quantum Calculation Worker or Research Agent), and manages the workflow execution.
"""

import logging
from pathlib import Path
from flask import current_app
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph_supervisor import create_supervisor
from langgraph_supervisor.handoff import create_forward_message_tool
from agent.utils import get_gemini_api_key, load_prompt
from agent.quantum_calculator import create_quantum_calculator
from agent.literature_surveyor import create_literature_surveyor
from agent.science_analyst import create_science_analyst
from agent.molecular_designer import create_molecular_designer

# Set up logging
logger = logging.getLogger(__name__)


def create_supervisor_agent():
    """
    Create a Supervisor agent using LangGraph's create_supervisor.

    The Supervisor coordinates between:
    - Quantum Calculator (quantum chemistry and molecular analysis)
    - Literature Surveyor (academic literature search)
    - Science Analyst (scientific report generation)
    - Molecular Designer (molecular structure generation and design)

    Returns:
        A compiled LangGraph Supervisor workflow ready for execution

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Supervisor.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    logger.info("Initializing Supervisor agent...")

    # Create worker agents
    try:
        quantum_calculator = create_quantum_calculator()
        logger.info("Quantum Calculator initialized")
    except Exception as e:
        logger.error(f"Failed to create Quantum Calculator: {e}", exc_info=True)
        raise

    try:
        literature_surveyor = create_literature_surveyor()
        logger.info("Literature Surveyor initialized")
    except Exception as e:
        logger.error(f"Failed to create Literature Surveyor: {e}", exc_info=True)
        raise

    try:
        science_analyst = create_science_analyst()
        logger.info("Science Analyst initialized")
    except Exception as e:
        logger.error(f"Failed to create Science Analyst: {e}", exc_info=True)
        raise

    try:
        molecular_designer = create_molecular_designer()
        logger.info("Molecular Designer initialized")
    except Exception as e:
        logger.error(f"Failed to create Molecular Designer: {e}", exc_info=True)
        raise

    # Initialize Supervisor LLM with configured model
    model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
    supervisor_llm = ChatGoogleGenerativeAI(
        model=model_name,
        api_key=api_key
    )
    logger.info(f"Supervisor LLM initialized with {model_name}")

    # Create forward message tool to prevent supervisor from paraphrasing worker responses
    forward_tool = create_forward_message_tool("supervisor")
    logger.info("Created forward_message tool to bypass supervisor paraphrasing")

    # Load system prompt using unified utility function
    prompt_dir = Path(__file__).parent / "prompts"
    system_prompt = load_prompt(prompt_dir, "system_prompt.txt")
    if not system_prompt:
        # Fallback prompt if file not found
        system_prompt = (
            "You are a Supervisor coordinating specialized AI agents for molecular science "
            "and quantum computational chemistry research. Analyze user requests, delegate to the "
            "appropriate specialist (quantum_calculator, literature_surveyor, or science_analyst), "
            "and orchestrate multi-step workflows. Use the forward_message tool to relay worker "
            "responses without modification."
        )
        logger.warning("Using fallback system prompt for Supervisor")

    # Create supervisor workflow using create_supervisor
    workflow = create_supervisor(
        [quantum_calculator, literature_surveyor, science_analyst, molecular_designer],
        model=supervisor_llm,
        prompt=system_prompt,
        tools=[forward_tool],  # Add forward_message tool
    )

    logger.info("Supervisor workflow created successfully with forward_message tool")
    return workflow


def get_compiled_supervisor():
    """
    Get a compiled Supervisor agent ready for execution.

    This is the main entry point for the API to obtain the supervisor graph executor.

    Returns:
        Compiled LangGraph Supervisor application
    """
    workflow = create_supervisor_agent()
    app = workflow.compile()
    logger.info("Supervisor compiled and ready for execution")
    return app
