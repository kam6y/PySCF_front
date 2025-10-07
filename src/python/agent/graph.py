"""
LangGraph Multi-Agent Dispatcher

This module implements the core multi-agent architecture using LangGraph.
It provides a stateful graph that routes user queries to specialized agents
(Molecular Agent or Research Agent) based on intent classification.

Architecture:
- Router Node: Classifies user intent and determines which agent to invoke
- Molecular Agent Node: Handles quantum chemistry calculations and molecular analysis
- Research Agent Node: Handles academic paper search and literature reviews
- Stateful Graph: Maintains conversation history and context across agent invocations

The graph uses LangGraph's conditional edges to implement dynamic routing,
allowing the system to intelligently dispatch tasks to the appropriate specialist.
"""

from typing import TypedDict, List, Literal, Dict, Any, Iterator, Annotated
import logging
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from langchain_core.pydantic_v1 import BaseModel, Field
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.graph import StateGraph, END
from langgraph.graph.message import add_messages
from agent.molecular_agent import get_gemini_api_key

# Set up logging
logger = logging.getLogger(__name__)


# ============================================================================
# State Definition
# ============================================================================

class AgentState(TypedDict):
    """
    Shared state object for the multi-agent graph.

    This state is passed between all nodes in the graph and maintains
    the conversation context.

    Attributes:
        messages: List of conversation messages (user inputs, AI responses, tool outputs)
                  Annotated with add_messages to enable proper LangGraph message accumulation
        next_node: The name of the next node to execute ("molecular", "research", or "end")
    """
    messages: Annotated[List[BaseMessage], add_messages]
    next_node: str


# ============================================================================
# Router Definition
# ============================================================================

class RouteQuery(BaseModel):
    """
    Structured output for the router's decision.

    This is used with Gemini's structured output feature to ensure
    the router returns a valid agent selection.
    """
    destination: Literal["molecular", "research"] = Field(
        description=(
            "The specialist agent to route the query to:\n"
            "- 'molecular': For quantum chemistry calculations, molecular properties, "
            "orbital analysis, geometry optimization, etc.\n"
            "- 'research': For finding academic papers, literature reviews, "
            "and research summaries."
        )
    )


def router_node(state: AgentState) -> AgentState:
    """
    Router node that classifies user intent and determines the next agent.

    This node uses Gemini's structured output to classify the user's query
    and decide whether to route to the Molecular Agent or Research Agent.

    Args:
        state: Current agent state containing conversation messages

    Returns:
        Updated state with next_node field set to the chosen agent
    """
    # Create an LLM with structured output for reliable routing
    api_key = get_gemini_api_key()
    llm_router = ChatGoogleGenerativeAI(
        model="gemini-2.5-flash",
        api_key=api_key
    ).with_structured_output(RouteQuery)

    # Get the latest user message
    user_message = state["messages"][-1].content if state["messages"] else ""

    # Invoke the router LLM with classification instructions
    result = llm_router.invoke(
        f"Given the user query: '{user_message}', which agent should handle it?\n\n"
        "Use 'molecular' for:\n"
        "- Quantum chemistry calculations (DFT, HF, MP2, CASSCF, etc.)\n"
        "- Molecular properties and orbital analysis\n"
        "- Geometry optimization\n"
        "- Vibrational frequency analysis\n"
        "- Any chemistry computation or molecular structure tasks\n\n"
        "Use 'research' for:\n"
        "- Finding academic papers\n"
        "- Literature searches\n"
        "- Research summaries and reviews\n"
        "- Questions about scientific publications"
    )

    # Update state with routing decision
    state["next_node"] = result.destination
    return state


# ============================================================================
# Message Conversion Utilities
# ============================================================================

def _convert_langchain_to_dict_format(
    messages: List[BaseMessage]
) -> List[Dict[str, Any]]:
    """
    Convert LangChain BaseMessage list to MolecularAgent dictionary format.

    LangChain format:
        [HumanMessage(content="..."), AIMessage(content="...")]

    MolecularAgent format:
        [
            {"role": "user", "parts": [{"text": "..."}]},
            {"role": "model", "parts": [{"text": "..."}]}
        ]

    Args:
        messages: List of LangChain BaseMessage objects

    Returns:
        List of dictionaries in MolecularAgent format
    """
    from typing import Dict, Any
    
    converted = []

    for msg in messages:
        if isinstance(msg, HumanMessage):
            role = "user"
        elif isinstance(msg, AIMessage):
            role = "model"
        else:
            # Skip SystemMessage or other message types
            logger.debug(f"Skipping message type: {type(msg).__name__}")
            continue

        converted.append({
            "role": role,
            "parts": [{"text": msg.content}]
        })

    return converted


# ============================================================================
# Agent Nodes
# ============================================================================

def molecular_agent_node(state: AgentState) -> AgentState:
    """
    Molecular Agent node - handles quantum chemistry calculations and molecular analysis.

    This node integrates the existing MolecularAgent into the LangGraph workflow.
    It converts LangChain messages to the format expected by MolecularAgent,
    invokes the agent with streaming, and collects the complete response.

    Args:
        state: Current agent state containing conversation messages

    Returns:
        Updated state with the molecular agent's response appended to messages

    Raises:
        Logs errors but does not raise exceptions; adds error messages to state instead
    """
    from agent.molecular_agent import MolecularAgent

    logger.info("Executing molecular_agent_node")

    # Extract messages and current user query
    messages = state["messages"]
    if not messages:
        logger.warning("No messages in state, returning empty response")
        state["messages"].append(AIMessage(content="メッセージが見つかりませんでした。"))
        return state

    user_message = messages[-1].content

    # Convert LangChain message history to MolecularAgent format (dict list)
    history = _convert_langchain_to_dict_format(messages[:-1])
    logger.debug(f"Converted {len(messages)-1} messages to MolecularAgent format")

    try:
        # Initialize and invoke MolecularAgent
        agent = MolecularAgent()

        # Collect streaming response chunks
        response_chunks = []
        for chunk in agent.chat(user_message, history):
            response_chunks.append(chunk)

        # Combine all chunks into complete response
        full_response = "".join(response_chunks)

        if not full_response.strip():
            logger.warning("MolecularAgent returned empty response")
            full_response = "応答を生成できませんでした。もう一度お試しください。"

        # Append complete response to state
        state["messages"].append(AIMessage(content=full_response))
        logger.info(f"MolecularAgent response added to state (length: {len(full_response)})")

    except Exception as e:
        logger.error(f"Error in molecular_agent_node: {e}", exc_info=True)
        error_message = (
            f"分子エージェントでエラーが発生しました: {str(e)}\n\n"
            "APIキーが正しく設定されているか確認してください。"
        )
        state["messages"].append(AIMessage(content=error_message))

    return state


def research_agent_node(state: AgentState) -> AgentState:
    """
    Research Agent node - handles academic paper search and literature reviews.

    This node integrates the ResearchAgent runnable into the LangGraph workflow.
    It invokes the agent with the user's query and LangChain automatically
    handles tool calling (search_arxiv).

    Args:
        state: Current agent state containing conversation messages

    Returns:
        Updated state with the research agent's response appended to messages

    Raises:
        Logs errors but does not raise exceptions; adds error messages to state instead
    """
    from research.agent import create_research_agent_runnable

    logger.info("Executing research_agent_node")

    # Extract current user query
    messages = state["messages"]
    if not messages:
        logger.warning("No messages in state, returning empty response")
        state["messages"].append(AIMessage(content="メッセージが見つかりませんでした。"))
        return state

    user_message = messages[-1].content

    try:
        # Create and invoke research agent runnable
        agent_runnable = create_research_agent_runnable()
        logger.debug(f"Invoking research agent with query: {user_message[:100]}...")

        # Invoke agent (ReAct agent automatically executes tool calls)
        # The agent returns a dict with 'messages' key containing the conversation
        response = agent_runnable.invoke({"messages": [("user", user_message)]})

        # Extract the final AI message from the response
        # create_react_agent returns: {"messages": [..., AIMessage(content="final response")]}
        if isinstance(response, dict) and "messages" in response:
            final_message = response["messages"][-1]
            if hasattr(final_message, 'content'):
                response_text = final_message.content
            else:
                response_text = str(final_message)
        elif hasattr(response, 'content'):
            response_text = response.content
        elif isinstance(response, str):
            response_text = response
        else:
            # Fallback: convert to string
            response_text = str(response)
            logger.warning(f"Unexpected response type: {type(response).__name__}")

        if not response_text.strip():
            logger.warning("ResearchAgent returned empty response")
            response_text = "論文検索結果が見つかりませんでした。別のキーワードでお試しください。"

        # Append response to state
        state["messages"].append(AIMessage(content=response_text))
        logger.info(f"ResearchAgent response added to state (length: {len(response_text)})")

    except Exception as e:
        logger.error(f"Error in research_agent_node: {e}", exc_info=True)
        error_message = (
            f"リサーチエージェントでエラーが発生しました: {str(e)}\n\n"
            "APIキーが正しく設定されているか、またはarXivサービスが利用可能か確認してください。"
        )
        state["messages"].append(AIMessage(content=error_message))

    return state


# ============================================================================
# Graph Construction
# ============================================================================

def create_dispatcher_graph():
    """
    Create and compile the LangGraph multi-agent dispatcher.

    Returns:
        Compiled LangGraph application ready for execution
    """
    # Initialize the state graph
    workflow = StateGraph(AgentState)

    # Add nodes
    workflow.add_node("router", router_node)
    workflow.add_node("molecular", molecular_agent_node)
    workflow.add_node("research", research_agent_node)

    # Set entry point
    workflow.set_entry_point("router")

    # Add conditional edges from router to specialist agents
    workflow.add_conditional_edges(
        "router",
        lambda state: state["next_node"],
        {
            "molecular": "molecular",
            "research": "research",
        }
    )

    # Both specialist agents return to END after execution
    workflow.add_edge("molecular", END)
    workflow.add_edge("research", END)

    # Compile and return the graph
    return workflow.compile()


# ============================================================================
# Utility Functions
# ============================================================================

def get_compiled_graph():
    """
    Get a compiled instance of the dispatcher graph.

    This is the main entry point for the API to obtain the graph executor.

    Returns:
        Compiled LangGraph application
    """
    return create_dispatcher_graph()
