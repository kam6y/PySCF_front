"""
Quantum Calculator Agent

This module implements the Quantum Calculator agent, which specializes in
quantum chemistry calculations, molecular analysis, and quantum computational chemistry tasks.

The agent uses LangGraph's create_react_agent pattern with Google Gemini for
intelligent task execution and tool calling.
"""

import logging
from pathlib import Path
from flask import current_app
from langchain_core.messages import HumanMessage, AIMessage
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.prebuilt import create_react_agent
from langgraph.graph import StateGraph, START, END, MessagesState

from agent.utils import get_gemini_api_key, extract_language_from_messages, get_language_name, load_prompt

# Set up logging
logger = logging.getLogger(__name__)


def _initialize_tools():
    """Initialize and return the list of available tools for the Quantum Calculator."""
    from agent.common_tools.execution_tools import (
        start_quantum_calculation,
        delete_calculation,
        convert_smiles_to_xyz,
        validate_xyz_format,
        search_pubchem_by_name,
        update_app_settings,
    )
    from agent.common_tools.analysis_tools import (
        list_all_calculations,
        get_calculation_details,
        get_supported_parameters,
        get_app_settings,
        get_system_resources,
    )

    return [
        # Calculation execution and management (execution tools)
        start_quantum_calculation,
        delete_calculation,  # HIL-enabled destructive tool

        # Molecular structure tools (execution tools)
        search_pubchem_by_name,
        convert_smiles_to_xyz,
        validate_xyz_format,

        # Calculation data retrieval (analysis tools)
        list_all_calculations,
        get_calculation_details,

        # System and settings (analysis + execution tools)
        get_supported_parameters,
        get_app_settings,
        update_app_settings,
        get_system_resources,
    ]


def _create_supervisor_wrapper(core_agent):
    """
    Create a wrapper that makes the Quantum Calculator compatible with Supervisor's message-based interface.

    This wrapper:
    1. Receives messages from Supervisor
    2. Detects user's language from the latest message
    3. Enhances the system prompt with language requirements
    4. Invokes the core ReAct agent
    5. Returns results as a message
    """
    def wrapper_node(state: MessagesState) -> dict:
        """Wrapper node that detects language and enhances prompts."""
        try:
            # Get messages from Supervisor
            messages = state.get("messages", [])
            if not messages:
                logger.error("No messages received from Supervisor")
                return {"messages": [AIMessage(content="Error: No query provided", name="quantum_calculator")]}

            # Filter for HumanMessages only (exclude ToolMessage forwarding)
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]
            if not human_messages:
                logger.error("No HumanMessage found in message history")
                return {"messages": [AIMessage(content="Error: 有効なクエリが見つかりませんでした", name="quantum_calculator")]}

            user_query = human_messages[-1].content
            logger.info(f"Quantum Calculator received query: {user_query}")

            # Extract user's language from Supervisor's marker message
            user_language = extract_language_from_messages(messages)
            language_name = get_language_name(user_language)
            logger.info(f"Using language: {user_language} ({language_name})")

            # Enhance messages with language instruction
            # Create a system message with language requirement
            language_instruction = f"""
## USER LANGUAGE REQUIREMENT (CRITICAL)
**You MUST write your ENTIRE response in: {language_name}**
**Language Code: {user_language}**
**This is the language the user used in their query.**

When executing calculations, providing data, or explaining results:
- Write ALL content in {language_name}
- Use appropriate technical terminology for {language_name}
- Maintain the same language throughout the entire response
"""

            # Pass to the core agent with enhanced context
            enhanced_messages = [
                AIMessage(content=language_instruction, name="system"),
                *messages
            ]

            result = core_agent.invoke({"messages": enhanced_messages})

            # Extract the final response
            response_messages = result.get("messages", [])
            if response_messages:
                # Get the last AI message
                final_message = response_messages[-1]
                logger.info(f"Quantum Calculator completed. Response length: {len(final_message.content)} characters")
                return {"messages": [AIMessage(content=final_message.content, name="quantum_calculator")]}
            else:
                logger.warning("No response from core agent")
                return {"messages": [AIMessage(content="No response generated", name="quantum_calculator")]}

        except Exception as e:
            logger.error(f"Error in Quantum Calculator wrapper: {str(e)}", exc_info=True)
            error_message = f"Error during calculation execution: {str(e)}"
            return {"messages": [AIMessage(content=error_message, name="quantum_calculator")]}

    # Build a simple graph with the wrapper node
    wrapper_builder = StateGraph(MessagesState)
    wrapper_builder.add_node("quantum_calculator", wrapper_node)
    wrapper_builder.add_edge(START, "quantum_calculator")
    wrapper_builder.add_edge("quantum_calculator", END)

    return wrapper_builder.compile(name="quantum_calculator")


def create_quantum_calculator():
    """
    Create a Quantum Calculator agent compatible with LangGraph Supervisor.

    This function creates the core Quantum Calculator agent and wraps it in a
    message-based interface that the Supervisor can interact with, with
    automatic language detection and adaptation.

    The agent specializes in:
    - Quantum chemistry calculations (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF)
    - Molecular structure analysis and visualization in user's language
    - Geometry optimization and frequency analysis
    - Molecular orbital and NTO analysis
    - IR spectrum generation

    Returns:
        A compiled LangGraph application compatible with Supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Create the core Quantum Calculator agent
    core_agent = _create_core_quantum_calculator()

    # Wrap it for Supervisor compatibility with language detection
    supervisor_compatible_agent = _create_supervisor_wrapper(core_agent)

    logger.info("Quantum Calculator with Supervisor wrapper and language detection created successfully")
    return supervisor_compatible_agent


def _create_core_quantum_calculator():
    """
    Create the core Quantum Calculator agent (internal implementation).

    This is the original create_quantum_calculator function, now renamed
    to distinguish it from the public API.
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Quantum Calculator.")
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
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Quantum Calculator")

    # Load system prompt using unified utility function
    prompt_dir = Path(__file__).parent / "prompts"
    system_prompt = load_prompt(prompt_dir, "system_prompt.txt")
    if not system_prompt:
        # Fallback prompt if file not found
        system_prompt = (
            "You are a Quantum Calculator, an AI assistant specialized in molecular analysis "
            "and quantum chemistry. You have access to tools for calculation management, molecular "
            "structure analysis, orbital visualization, and more. Use them proactively to help users "
            "with their quantum computational chemistry tasks."
        )
        logger.warning("Using fallback system prompt for Quantum Calculator")

    # Initialize tools
    tools = _initialize_tools()
    logger.debug(f"Initialized {len(tools)} tools for Quantum Calculator")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="quantum_calculator",
        prompt=system_prompt
    )

    logger.info("Core Quantum Calculator agent created successfully")
    return agent
