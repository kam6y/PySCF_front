"""
Research Agent Module

This module implements the Research Agent, which specializes in searching and
summarizing academic papers from sources like arXiv.

The Research Agent uses LangChain and Google's Gemini model to provide intelligent
responses to research-related queries.
"""

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.prompts import ChatPromptTemplate
from .tools import search_arxiv


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
    Create a LangChain runnable for the research agent.

    Returns:
        A LangChain runnable chain that processes research queries
    """
    prompt = ChatPromptTemplate.from_messages([
        ("system", RESEARCH_AGENT_PROMPT),
        ("human", "{input}")
    ])

    # Use Gemini 1.5 Flash for fast, cost-effective research queries
    llm = ChatGoogleGenerativeAI(model="gemini-1.5-flash")

    # Bind the arXiv search tool to the LLM
    llm_with_tools = llm.bind_tools([search_arxiv])

    return prompt | llm_with_tools
