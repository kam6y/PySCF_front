"""
Science Analyst Agent

This module implements the Science Analyst agent, which specializes in creating
comprehensive scientific reports by synthesizing information from other specialized agents.

The agent uses LangGraph's create_react_agent pattern with Google Gemini for
intelligent report generation and analysis.
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
    """Initialize and return the list of available tools for the Science Analyst."""
    from agent.common_tools.analysis_tools import (
        get_calculation_details,
        get_molecular_orbitals,
        generate_orbital_cube,
        list_cube_files,
        generate_ir_spectrum,
    )

    return [
        # Calculation data retrieval (analysis tools only)
        get_calculation_details,

        # Molecular orbital analysis (analysis tools only)
        get_molecular_orbitals,
        generate_orbital_cube,
        list_cube_files,
        # Note: delete_cube_files removed - Science Analyst focuses on analysis, not deletion

        # Spectroscopy (analysis tools only)
        generate_ir_spectrum,
    ]


def _create_supervisor_wrapper(core_agent):
    """
    Create a wrapper that makes the Science Analyst compatible with Supervisor's message-based interface.

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
                return {"messages": [AIMessage(content="Error: No query provided", name="science_analyst")]}

            # Filter for HumanMessages only (exclude ToolMessage forwarding)
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]
            if not human_messages:
                logger.error("No HumanMessage found in message history")
                return {"messages": [AIMessage(content="Error: 有効なクエリが見つかりませんでした", name="science_analyst")]}

            user_query = human_messages[-1].content
            logger.info(f"Science Analyst received query: {user_query}")

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
                logger.info(f"Science Analyst completed. Response length: {len(final_message.content)} characters")
                return {"messages": [AIMessage(content=final_message.content, name="science_analyst")]}
            else:
                logger.warning("No response from core agent")
                return {"messages": [AIMessage(content="No response generated", name="science_analyst")]}

        except Exception as e:
            logger.error(f"Error in Science Analyst wrapper: {str(e)}", exc_info=True)
            error_message = f"Error during report generation: {str(e)}"
            return {"messages": [AIMessage(content=error_message, name="science_analyst")]}

    # Build a simple graph with the wrapper node
    wrapper_builder = StateGraph(MessagesState)
    wrapper_builder.add_node("science_analyst", wrapper_node)
    wrapper_builder.add_edge(START, "science_analyst")
    wrapper_builder.add_edge("science_analyst", END)

    return wrapper_builder.compile(name="science_analyst")


def create_science_analyst():
    """
    Create a Science Analyst agent compatible with LangGraph Supervisor.

    This function creates the core Science Analyst agent and wraps it in a
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
    # Create the core Science Analyst agent
    core_agent = _create_core_science_analyst()

    # Wrap it for Supervisor compatibility with language detection
    supervisor_compatible_agent = _create_supervisor_wrapper(core_agent)

    logger.info("Science Analyst with Supervisor wrapper and language detection created successfully")
    return supervisor_compatible_agent


def _create_core_science_analyst():
    """
    Create the core Science Analyst agent (internal implementation).

    This is the original create_science_analyst function, now renamed
    to distinguish it from the public API.
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Science Analyst.")
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
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Science Analyst")

    # Load system prompt using unified utility function
    prompt_dir = Path(__file__).parent / "prompts"
    system_prompt = load_prompt(prompt_dir, "system_prompt.txt")
    if not system_prompt:
        # Fallback prompt if file not found
        system_prompt = (
            "You are a Science Analyst, an AI assistant specialized in creating comprehensive "
            "scientific reports for molecular science and quantum computational chemistry. "
            "You synthesize information from other agents and produce well-structured, "
            "scientifically accurate reports in Japanese or English based on user preference."
        )
        logger.warning("Using fallback system prompt for Science Analyst")

    # Initialize tools for accessing and analyzing calculation results
    tools = _initialize_tools()
    logger.debug(f"Science Analyst initialized with {len(tools)} tools for result analysis")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="science_analyst",
        prompt=system_prompt
    )

    logger.info("Core Science Analyst agent created successfully")
    return agent
