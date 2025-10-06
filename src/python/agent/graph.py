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

from typing import TypedDict, List, Literal
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from langchain_core.pydantic_v1 import BaseModel, Field
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.graph import StateGraph, END


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
        next_node: The name of the next node to execute ("molecular", "research", or "end")
    """
    messages: List[BaseMessage]
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
    llm_router = ChatGoogleGenerativeAI(
        model="gemini-1.5-flash"
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
# Agent Nodes (Placeholders for Phase 1)
# ============================================================================

def molecular_agent_node(state: AgentState) -> AgentState:
    """
    Molecular Agent node (placeholder).

    This will invoke the existing MolecularAgent logic.
    To be implemented in Phase 3.

    Args:
        state: Current agent state

    Returns:
        Updated state with molecular agent's response
    """
    # TODO: Implement in Phase 3 - Integration with existing MolecularAgent
    # For now, return a placeholder message
    state["messages"].append(
        AIMessage(content="[Molecular Agent - To be implemented in Phase 3]")
    )
    state["next_node"] = "end"
    return state


def research_agent_node(state: AgentState) -> AgentState:
    """
    Research Agent node (placeholder).

    This will invoke the ResearchAgent logic defined in research/agent.py.
    To be implemented in Phase 2.

    Args:
        state: Current agent state

    Returns:
        Updated state with research agent's response
    """
    # TODO: Implement in Phase 2 - Connect to ResearchAgent
    # For now, return a placeholder message
    state["messages"].append(
        AIMessage(content="[Research Agent - To be implemented in Phase 2]")
    )
    state["next_node"] = "end"
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
