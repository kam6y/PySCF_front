"""
Molecular Designer Agent

This module implements the Molecular Designer agent, which specializes in
generating novel molecular structures using cheminformatics tools (RDKit).

The agent uses LangGraph's create_react_agent pattern with Google Gemini for
intelligent molecular design and structure-property reasoning.
"""

import logging
from pathlib import Path
from flask import current_app
from langchain_core.messages import HumanMessage, AIMessage
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.prebuilt import create_react_agent
from langgraph.graph import StateGraph, START, END, MessagesState

from agent.utils import get_gemini_api_key, detect_language_cached, get_language_name, load_prompt

# Set up logging
logger = logging.getLogger(__name__)


def _initialize_tools():
    """Initialize and return the list of available tools for the Molecular Designer."""
    from . import tools

    return [
        # PubChem SMILES search (CRITICAL - must use this for compound names)
        tools.search_compound_smiles_pubchem,
        # Molecular generation
        tools.generate_analogs_rdkit,
        tools.predict_simple_properties_rdkit,
        # Fragment-based design
        tools.brics_decompose_rdkit,
        tools.recap_decompose_rdkit,
        # Systematic substitution
        tools.substitute_side_chains_rdkit,
    ]


def _create_supervisor_wrapper(core_agent):
    """
    Create a wrapper that makes the Molecular Designer compatible with Supervisor's message-based interface.

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
                return {"messages": [AIMessage(content="Error: No query provided", name="molecular_designer")]}

            # Filter for HumanMessages only (exclude ToolMessage forwarding)
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]
            if not human_messages:
                logger.error("No HumanMessage found in message history")
                return {"messages": [AIMessage(content="Error: 有効なクエリが見つかりませんでした", name="molecular_designer")]}

            user_query = human_messages[-1].content
            logger.info(f"Molecular Designer received query: {user_query}")

            # Detect user's language using LLM with caching
            api_key = get_gemini_api_key()
            model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
            llm = ChatGoogleGenerativeAI(model=model_name, api_key=api_key)
            user_language = detect_language_cached(user_query, llm)
            language_name = get_language_name(user_language)
            logger.info(f"Detected user language: {user_language} ({language_name})")

            # Enhance messages with language instruction
            # Create a system message with language requirement
            language_instruction = f"""
## USER LANGUAGE REQUIREMENT (CRITICAL)
**You MUST write your ENTIRE response in: {language_name}**
**Language Code: {user_language}**
**This is the language the user used in their query.**

When designing molecules, providing rationale, or explaining design choices:
- Write ALL content in {language_name}
- Use appropriate chemical terminology for {language_name}
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
                logger.info(f"Molecular Designer completed. Response length: {len(final_message.content)} characters")
                return {"messages": [AIMessage(content=final_message.content, name="molecular_designer")]}
            else:
                logger.warning("No response from core agent")
                return {"messages": [AIMessage(content="No response generated", name="molecular_designer")]}

        except Exception as e:
            logger.error(f"Error in Molecular Designer wrapper: {str(e)}", exc_info=True)
            error_message = f"Error during molecular design: {str(e)}"
            return {"messages": [AIMessage(content=error_message, name="molecular_designer")]}

    # Build a simple graph with the wrapper node
    wrapper_builder = StateGraph(MessagesState)
    wrapper_builder.add_node("molecular_designer", wrapper_node)
    wrapper_builder.add_edge(START, "molecular_designer")
    wrapper_builder.add_edge("molecular_designer", END)

    return wrapper_builder.compile(name="molecular_designer")


def create_molecular_designer():
    """
    Create a Molecular Designer agent compatible with LangGraph Supervisor.

    This function creates the core Molecular Designer agent and wraps it in a
    message-based interface that the Supervisor can interact with, with
    automatic language detection and adaptation.

    The agent specializes in:
    - Generating novel molecular structures using RDKit
    - Proposing analogs based on design strategies (conjugation, push-pull, heteroatom)
    - Predicting simple molecular properties for filtering
    - Providing design rationale based on structure-property relationships

    Returns:
        A compiled LangGraph application compatible with Supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Create the core Molecular Designer agent
    core_agent = _create_core_molecular_designer()

    # Wrap it for Supervisor compatibility with language detection
    supervisor_compatible_agent = _create_supervisor_wrapper(core_agent)

    logger.info("Molecular Designer with Supervisor wrapper and language detection created successfully")
    return supervisor_compatible_agent


def _create_core_molecular_designer():
    """
    Create the core Molecular Designer agent (internal implementation).

    This is the original create_molecular_designer function, now renamed
    to distinguish it from the public API.
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Molecular Designer.")
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
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Molecular Designer")

    # Load system prompt using unified utility function
    prompt_dir = Path(__file__).parent / "prompts"
    system_prompt = load_prompt(prompt_dir, "system_prompt.txt")
    if not system_prompt:
        # Fallback prompt if file not found
        system_prompt = (
            "You are a Molecular Designer, an AI assistant specialized in generating novel molecular "
            "structures using cheminformatics tools. You have access to RDKit-based tools for molecular "
            "generation and property prediction. Use them proactively to propose candidate molecules "
            "based on user design goals. Always provide SMILES format output with design rationale."
        )
        logger.warning("Using fallback system prompt for Molecular Designer")

    # Initialize tools
    tools = _initialize_tools()
    logger.debug(f"Initialized {len(tools)} tools for Molecular Designer")

    # Create ReAct agent using LangGraph prebuilt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        name="molecular_designer",
        prompt=system_prompt
    )

    logger.info("Core Molecular Designer agent created successfully")
    return agent
