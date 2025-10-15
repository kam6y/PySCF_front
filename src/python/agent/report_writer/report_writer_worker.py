"""
Report Writer Worker

This module implements the Report Writer agent, which specializes in creating
comprehensive scientific reports by synthesizing information from other specialized agents.

The worker uses LangGraph's create_react_agent pattern with Google Gemini for
intelligent report generation without requiring external tools.
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


def _initialize_tools():
    """Initialize and return the list of available tools for the Report Writer."""
    from . import tools

    return [
        # Calculation data retrieval
        tools.get_calculation_details,

        # Molecular orbital analysis
        tools.get_molecular_orbitals,
        tools.generate_orbital_cube,
        tools.list_cube_files,
        tools.delete_cube_files,

        # Spectroscopy
        tools.generate_ir_spectrum,
    ]


def _create_supervisor_wrapper(core_agent):
    """
    Create a wrapper that makes the Report Writer compatible with Supervisor's message-based interface.

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
                return {"messages": [AIMessage(content="Error: No query provided", name="report_writer")]}

            # Filter for HumanMessages only (exclude ToolMessage forwarding)
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]
            if not human_messages:
                logger.error("No HumanMessage found in message history")
                return {"messages": [AIMessage(content="Error: 有効なクエリが見つかりませんでした", name="report_writer")]}

            user_query = human_messages[-1].content
            logger.info(f"Report Writer received query: {user_query}")

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

When creating reports, analyzing data, or explaining results:
- Write ALL content in {language_name}
- Use appropriate scientific terminology for {language_name}
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
                logger.info(f"Report Writer completed. Response length: {len(final_message.content)} characters")
                return {"messages": [AIMessage(content=final_message.content, name="report_writer")]}
            else:
                logger.warning("No response from core agent")
                return {"messages": [AIMessage(content="No response generated", name="report_writer")]}

        except Exception as e:
            logger.error(f"Error in Report Writer wrapper: {str(e)}", exc_info=True)
            error_message = f"Error during report generation: {str(e)}"
            return {"messages": [AIMessage(content=error_message, name="report_writer")]}

    # Build a simple graph with the wrapper node
    wrapper_builder = StateGraph(MessagesState)
    wrapper_builder.add_node("report", wrapper_node)
    wrapper_builder.add_edge(START, "report")
    wrapper_builder.add_edge("report", END)

    return wrapper_builder.compile(name="report_writer")


def create_report_writer():
    """
    Create a Report Writer agent compatible with LangGraph Supervisor.

    This function creates the core Report Writer agent and wraps it in a
    message-based interface that the Supervisor can interact with, with
    automatic language detection and adaptation.

    The agent specializes in:
    - Retrieving and interpreting quantum chemistry calculation results
    - Analyzing molecular orbitals and spectroscopy data
    - Creating comprehensive scientific reports in user's language
    - Synthesizing information from multiple sources
    - Generating calculation reports and literature reviews
    - Formatting professional documentation

    Returns:
        A compiled LangGraph application compatible with Supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Create the core Report Writer agent
    core_agent = _create_core_report_writer()

    # Wrap it for Supervisor compatibility with language detection
    supervisor_compatible_agent = _create_supervisor_wrapper(core_agent)

    logger.info("Report Writer with Supervisor wrapper and language detection created successfully")
    return supervisor_compatible_agent


def _create_core_report_writer():
    """
    Create the core Report Writer agent (internal implementation).

    This is the original create_report_writer function, now renamed
    to distinguish it from the public API.
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Report Writer.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    # Initialize LLM with configured model
    # Note: Consider using gemini-2.5-pro for higher quality reports if needed
    model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
    llm = ChatGoogleGenerativeAI(
        model=model_name,
        api_key=api_key
    )
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Report Writer")

    # Load system prompt
    system_prompt = _load_system_prompt()

    # Initialize tools for accessing and analyzing calculation results
    tools = _initialize_tools()
    logger.debug(f"Report Writer initialized with {len(tools)} tools for result analysis")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="report_writer",
        prompt=system_prompt
    )

    logger.info("Core Report Writer agent created successfully")
    return agent
