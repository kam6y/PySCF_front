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
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph_supervisor import create_supervisor
from langgraph_supervisor.handoff import create_handoff_tool
from agent.utils import get_gemini_api_key, load_prompt, detect_language_cached
from agent.quantum_calculator import create_quantum_calculator
from agent.literature_surveyor import create_literature_surveyor
from agent.science_analyst import create_science_analyst
from agent.molecular_designer import create_molecular_designer

# Set up logging
logger = logging.getLogger(__name__)


def _initialize_tools():
    """
    Initialize and return the list of available tools for the Supervisor.

    The Supervisor uses these tools to search for and identify calculations,
    allowing it to provide specific calculation IDs to worker agents.

    Returns:
        list: List of tool functions for the Supervisor
    """
    from agent.common_tools.analysis_tools import (
        list_all_calculations,
        find_calculations,
    )

    return [
        list_all_calculations,
        find_calculations,
    ]


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

    # Note: We no longer use forward_message tool
    # The supervisor now adds value by summarizing findings and providing context
    logger.info("Supervisor configured to add contextual value to worker responses")

    # Create custom handoff tools with explicit descriptions
    handoff_tools = [
        create_handoff_tool(
            agent_name="quantum_calculator",
            name="transfer_to_quantum_calculator",
            description=(
                "Transfer/delegate the task to quantum_calculator agent. "
                "Use this when the user requests quantum chemistry calculations, "
                "molecular structure preparation, or system management. "
                "IMPORTANT: Calling this tool executes the transfer - "
                "do not just describe the transfer in text."
            )
        ),
        create_handoff_tool(
            agent_name="literature_surveyor",
            name="transfer_to_literature_surveyor",
            description=(
                "Transfer/delegate the task to literature_surveyor agent. "
                "Use this when the user requests academic literature research. "
                "IMPORTANT: Calling this tool executes the transfer - "
                "do not just describe the transfer in text."
            )
        ),
        create_handoff_tool(
            agent_name="science_analyst",
            name="transfer_to_science_analyst",
            description=(
                "Transfer/delegate the task to science_analyst agent. "
                "Use this when the user requests result interpretation, "
                "report generation, or data analysis. "
                "IMPORTANT: Calling this tool executes the transfer - "
                "do not just describe the transfer in text."
            )
        ),
        create_handoff_tool(
            agent_name="molecular_designer",
            name="transfer_to_molecular_designer",
            description=(
                "Transfer/delegate the task to molecular_designer agent. "
                "Use this when the user requests molecular design or "
                "structure generation tasks. "
                "IMPORTANT: Calling this tool executes the transfer - "
                "do not just describe the transfer in text."
            )
        ),
    ]
    logger.info(f"Created {len(handoff_tools)} custom handoff tools")

    # Initialize Supervisor's tools for calculation search and identification
    supervisor_tools = _initialize_tools()
    supervisor_tools.extend(handoff_tools)  # Add custom handoff tools
    logger.info(f"Initialized {len(supervisor_tools)} tools for Supervisor (including custom handoffs)")

    # Load system prompt using unified utility function
    prompt_dir = Path(__file__).parent / "prompts"
    system_prompt = load_prompt(prompt_dir, "system_prompt.txt")
    if not system_prompt:
        # Fallback prompt if file not found
        system_prompt = (
            "You are a Supervisor coordinating specialized AI agents for molecular science "
            "and quantum computational chemistry research. Analyze user requests, delegate to the "
            "appropriate specialist (quantum_calculator, literature_surveyor, science_analyst, "
            "or molecular_designer), and orchestrate multi-step workflows. When workers complete "
            "their tasks, add context by summarizing key findings, explaining current status, "
            "suggesting next steps, and including the full worker response between *** separators."
        )
        logger.warning("Using fallback system prompt for Supervisor")

    # Create supervisor workflow using create_supervisor
    workflow = create_supervisor(
        [quantum_calculator, literature_surveyor, science_analyst, molecular_designer],
        model=supervisor_llm,
        prompt=system_prompt,
        tools=supervisor_tools,  # Add all Supervisor tools
    )

    logger.info("Supervisor workflow created successfully with intelligent coordination")
    return workflow


class LanguageDetectionWrapper:
    """
    Wrapper that detects user language once and shares it with all workers.

    This wrapper:
    1. Detects the language from the first HumanMessage
    2. Inserts a SystemMessage with "DETECTED_LANGUAGE:xx" marker
    3. Invokes the Supervisor workflow with the enhanced messages

    This eliminates redundant language detection across multiple workers.
    """
    def __init__(self, supervisor_app):
        """Initialize with the compiled Supervisor application."""
        self.supervisor_app = supervisor_app
        logger.info("Language detection wrapper initialized")

    def invoke(self, state, config=None, **kwargs):
        """
        Invoke the Supervisor with language detection.

        Args:
            state: MessagesState or dict with "messages" key
            config: Optional configuration dict
            **kwargs: Additional arguments to pass to the Supervisor

        Returns:
            The result from the Supervisor workflow
        """
        try:
            # Extract messages from state
            if isinstance(state, dict):
                messages = state.get("messages", [])
            else:
                messages = state.messages if hasattr(state, "messages") else []

            # Find the first HumanMessage to detect language
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]

            if human_messages:
                user_query = human_messages[0].content
                logger.info(f"Detecting language for query: {user_query[:100]}...")

                # Detect language using LLM with caching
                api_key = get_gemini_api_key()
                model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
                llm = ChatGoogleGenerativeAI(model=model_name, api_key=api_key)
                user_language = detect_language_cached(user_query, llm)

                logger.info(f"Detected user language: {user_language}")

                # Insert language marker as SystemMessage at the beginning
                language_marker = SystemMessage(
                    content=f"DETECTED_LANGUAGE:{user_language}",
                    name="supervisor_system"
                )

                # Create enhanced state with language marker
                enhanced_messages = [language_marker] + messages
                enhanced_state = {"messages": enhanced_messages} if isinstance(state, dict) else state.__class__(messages=enhanced_messages)
            else:
                logger.warning("No HumanMessage found for language detection")
                enhanced_state = state

            # Invoke the Supervisor workflow
            result = self.supervisor_app.invoke(enhanced_state, config=config, **kwargs)
            return result

        except Exception as e:
            logger.error(f"Error in language detection wrapper: {str(e)}", exc_info=True)
            # On error, continue without language detection
            return self.supervisor_app.invoke(state, config=config, **kwargs)

    def stream(self, state, config=None, **kwargs):
        """
        Stream the Supervisor execution with language detection.

        Args:
            state: MessagesState or dict with "messages" key
            config: Optional configuration dict
            **kwargs: Additional arguments to pass to the Supervisor (e.g., stream_mode)

        Yields:
            Stream chunks from the Supervisor workflow
        """
        try:
            # Extract messages from state
            if isinstance(state, dict):
                messages = state.get("messages", [])
            else:
                messages = state.messages if hasattr(state, "messages") else []

            # Find the first HumanMessage to detect language
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]

            if human_messages:
                user_query = human_messages[0].content
                logger.info(f"Detecting language for query (streaming): {user_query[:100]}...")

                # Detect language using LLM with caching
                api_key = get_gemini_api_key()
                model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
                llm = ChatGoogleGenerativeAI(model=model_name, api_key=api_key)
                user_language = detect_language_cached(user_query, llm)

                logger.info(f"Detected user language: {user_language}")

                # Insert language marker as SystemMessage at the beginning
                language_marker = SystemMessage(
                    content=f"DETECTED_LANGUAGE:{user_language}",
                    name="supervisor_system"
                )

                # Create enhanced state with language marker
                enhanced_messages = [language_marker] + messages
                enhanced_state = {"messages": enhanced_messages} if isinstance(state, dict) else state.__class__(messages=enhanced_messages)
            else:
                logger.warning("No HumanMessage found for language detection")
                enhanced_state = state

            # Stream the Supervisor workflow
            for chunk in self.supervisor_app.stream(enhanced_state, config=config, **kwargs):
                yield chunk

        except Exception as e:
            logger.error(f"Error in language detection wrapper (streaming): {str(e)}", exc_info=True)
            # On error, continue without language detection
            for chunk in self.supervisor_app.stream(state, config=config, **kwargs):
                yield chunk


def get_compiled_supervisor():
    """
    Get a compiled Supervisor agent ready for execution with language detection.

    This is the main entry point for the API to obtain the supervisor graph executor.
    The returned supervisor includes automatic language detection that runs once
    and shares the detected language with all workers.

    Returns:
        LanguageDetectionWrapper: Wrapper around the compiled Supervisor with language detection
    """
    # Create base Supervisor workflow
    base_workflow = create_supervisor_agent()
    base_app = base_workflow.compile()

    # Wrap with language detection functionality
    wrapped_app = LanguageDetectionWrapper(base_app)

    logger.info("Supervisor compiled with language detection wrapper and ready for execution")
    return wrapped_app
