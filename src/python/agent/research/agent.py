"""
Research Agent Module

This module implements the Research Agent, which specializes in searching and
summarizing academic papers from sources like arXiv.

The Research Agent uses LangChain and Google's Gemini model to provide intelligent
responses to research-related queries.
"""

from langchain_google_genai import ChatGoogleGenerativeAI

from .tools import search_arxiv
from agent.quantum_calc.quantum_calc_worker import get_gemini_api_key


# System prompt for the Research Agent
RESEARCH_AGENT_PROMPT = """You are a specialized research assistant for academic literature search.

Your primary role is to help users find and understand relevant academic papers,
particularly in the fields of:
- Quantum Chemistry
- Computational Chemistry
- Molecular Physics
- Density Functional Theory
- Electronic Structure Theory

When responding to queries:
1. Use the search_arxiv tool to find relevant papers
2. Summarize the key findings in a clear, accessible way
3. Provide clickable links to PDF versions of papers
4. If multiple papers are found, organize them by relevance
5. Format your responses using Markdown for readability

Always cite your sources and provide direct links to papers when available.
"""


def create_research_agent_runnable():
    """
    Create a LangChain ReAct agent for the research agent.

    This agent uses LangGraph's create_react_agent to automatically handle
    tool execution loops, allowing the LLM to call search_arxiv and process results.

    Returns:
        A compiled LangGraph agent that processes research queries with tool execution
    """
    from langgraph.prebuilt import create_react_agent
    
    # Use Gemini 2.5 Flash for fast, cost-effective research queries
    api_key = get_gemini_api_key()
    llm = ChatGoogleGenerativeAI(model="gemini-2.5-flash", api_key=api_key)

    # Create a ReAct agent with the arXiv search tool
    # The agent will automatically:
    # 1. Decide when to use tools
    # 2. Execute tool calls
    # 3. Process tool outputs
    # 4. Generate final response
    agent = create_react_agent(
        llm,
        tools=[search_arxiv],
        name="research_expert",
        prompt=RESEARCH_AGENT_PROMPT
    )

    return agent
