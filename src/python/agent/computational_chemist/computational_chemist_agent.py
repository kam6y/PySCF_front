"""
Computational Chemist Agent

This module implements the Computational Chemist agent, which specializes in
quantum chemistry calculations, molecular analysis, and computational chemistry tasks.

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

from agent.utils import get_gemini_api_key, detect_language, get_language_name

# Set up logging
logger = logging.getLogger(__name__)


def _load_system_prompt() -> str:
    """Load the system prompt for the Computational Chemist."""
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
        "You are a Computational Chemist, an AI assistant specialized in molecular analysis "
        "and quantum chemistry. You have access to tools for calculation management, molecular "
        "structure analysis, orbital visualization, and more. Use them proactively to help users "
        "with their computational chemistry tasks."
    )


def _initialize_tools():
    """Initialize and return the list of available tools for the Computational Chemist."""
    from . import tools

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

        # System and settings
        tools.get_supported_parameters,
        tools.get_app_settings,
        tools.update_app_settings,
        tools.get_system_resources,
    ]


def _create_supervisor_wrapper(core_agent):
    """
    Create a wrapper that makes the Computational Chemist compatible with Supervisor's message-based interface.

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
                return {"messages": [AIMessage(content="Error: No query provided", name="computational_chemist")]}

            # Filter for HumanMessages only (exclude ToolMessage forwarding)
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]
            if not human_messages:
                logger.error("No HumanMessage found in message history")
                return {"messages": [AIMessage(content="Error: 有効なクエリが見つかりませんでした", name="computational_chemist")]}

            user_query = human_messages[-1].content
            logger.info(f"Computational Chemist received query: {user_query}")

            # Detect user's language using LLM
            api_key = get_gemini_api_key()
            model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
            llm = ChatGoogleGenerativeAI(model=model_name, api_key=api_key)
            user_language = detect_language(user_query, llm)
            language_name = get_language_name(user_language)
            logger.info(f"Detected user language: {user_language} ({language_name})")

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
                logger.info(f"Computational Chemist completed. Response length: {len(final_message.content)} characters")
                return {"messages": [AIMessage(content=final_message.content, name="computational_chemist")]}
            else:
                logger.warning("No response from core agent")
                return {"messages": [AIMessage(content="No response generated", name="computational_chemist")]}

        except Exception as e:
            logger.error(f"Error in Computational Chemist wrapper: {str(e)}", exc_info=True)
            error_message = f"Error during calculation execution: {str(e)}"
            return {"messages": [AIMessage(content=error_message, name="computational_chemist")]}

    # Build a simple graph with the wrapper node
    wrapper_builder = StateGraph(MessagesState)
    wrapper_builder.add_node("computational_chemist", wrapper_node)
    wrapper_builder.add_edge(START, "computational_chemist")
    wrapper_builder.add_edge("computational_chemist", END)

    return wrapper_builder.compile(name="computational_chemist")


def create_computational_chemist():
    """
    Create a Computational Chemist agent compatible with LangGraph Supervisor.

    This function creates the core Computational Chemist agent and wraps it in a
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
    # Create the core Computational Chemist agent
    core_agent = _create_core_computational_chemist()

    # Wrap it for Supervisor compatibility with language detection
    supervisor_compatible_agent = _create_supervisor_wrapper(core_agent)

    logger.info("Computational Chemist with Supervisor wrapper and language detection created successfully")
    return supervisor_compatible_agent


def _create_core_computational_chemist():
    """
    Create the core Computational Chemist agent (internal implementation).

    This is the original create_computational_chemist function, now renamed
    to distinguish it from the public API.
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Computational Chemist.")
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
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Computational Chemist")

    # Load system prompt
    system_prompt = _load_system_prompt()

    # Initialize tools
    tools = _initialize_tools()
    logger.debug(f"Initialized {len(tools)} tools for Computational Chemist")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="computational_chemist",
        prompt=system_prompt
    )

    logger.info("Core Computational Chemist agent created successfully")
    return agent
