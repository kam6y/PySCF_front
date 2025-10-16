"""
Deep Research Agent Module

This module implements a Deep Research agent that performs iterative, in-depth
investigation of topics using multiple search tools (arXiv, Tavily, PubMed).

The agent uses LangGraph's StateGraph for custom iteration control, allowing it to:
1. Analyze current findings
2. Decide what to search next
3. Perform searches across multiple sources
4. Synthesize all findings into a comprehensive report

This replaces the simple ReAct agent with a more sophisticated deep research workflow.
"""

import logging
import json
import re
from pathlib import Path
from typing import Optional, Literal
from typing_extensions import TypedDict
from flask import current_app
from langchain_core.messages import HumanMessage, AIMessage
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.graph import StateGraph, START, END, MessagesState

from .tools import (
    search_arxiv,
    search_tavily,
    search_pubmed,
    search_chemrxiv,
    search_openalex,
    search_semantic_scholar
)
from agent.utils import get_gemini_api_key, detect_language, get_language_name

# Set up logging
logger = logging.getLogger(__name__)


# State schema for Deep Research workflow
class DeepResearchState(TypedDict):
    """State for the deep research iteration loop."""
    # Input
    topic: str                          # Original research topic
    depth: int                          # Maximum iterations
    user_language: str                  # Detected language code (e.g., 'ja', 'en', 'es')

    # Iteration management
    current_iteration: int              # Current iteration number
    searched_topics: list[str]          # Topics already searched
    findings: list[dict]                # Accumulated search results

    # Control flow
    next_search_topic: Optional[str]    # Next topic to search (None if done)
    should_continue: bool               # Whether to continue iterating
    semantic_scholar_rate_limited: bool # Track if Semantic Scholar rate limit was hit

    # Output
    final_report: str                   # Synthesized final report


def _load_prompt(filename: str) -> str:
    """Load a prompt from the prompts directory."""
    try:
        prompt_path = Path(__file__).parent / "prompts" / filename
        if prompt_path.exists():
            with open(prompt_path, 'r', encoding='utf-8') as f:
                return f.read()
        else:
            logger.warning(f"Prompt file not found: {prompt_path}")
            return ""
    except Exception as e:
        logger.error(f"Failed to load prompt {filename}: {e}")
        return ""


def _create_analyze_node(llm: ChatGoogleGenerativeAI):
    """
    Create the analysis node that decides what to search next.

    This node:
    1. Reviews current findings
    2. Identifies knowledge gaps
    3. Decides next search topic or concludes the research
    """
    analyze_prompt = _load_prompt("analyze_prompt.txt")

    def analyze_findings(state: DeepResearchState) -> dict:
        """Analyze current findings and decide next steps."""
        try:
            logger.info(f"Analyzing findings (iteration {state['current_iteration']}/{state['depth']})")

            # Check if we've reached max depth
            if state['current_iteration'] >= state['depth']:
                logger.info("Reached maximum depth, concluding research")
                return {
                    "should_continue": False,
                    "next_search_topic": None
                }

            # Prepare context for analysis
            findings_summary = "\n\n".join([
                f"Search {i+1} (Topic: {f['topic']}):\n{f['results']}"
                for i, f in enumerate(state['findings'])
            ])

            searched_topics_str = ", ".join(state['searched_topics']) if state['searched_topics'] else "None"

            # Create analysis prompt
            prompt = f"""{analyze_prompt}

## Current Research Context

**Original Topic:** {state['topic']}

**Already Searched Topics:** {searched_topics_str}

**Current Findings:**
{findings_summary if findings_summary else "No findings yet."}

## Your Task

Analyze the current state and decide:
1. Is more information needed?
2. If yes, what specific aspect should be researched next (in English)?
3. If no, we will proceed to synthesis.

Respond in JSON format:
{{
    "nextSearchTopic": "specific topic to search (in English)" or null,
    "shouldContinue": true or false,
    "reasoning": "brief explanation of your decision"
}}

IMPORTANT:
- Do NOT output a topic that is exactly the same as any topic in "Already Searched Topics"
- Make nextSearchTopic more specific or explore a different angle
- If sufficient information has been gathered, set shouldContinue to false
"""

            # Get LLM decision with JSON response
            response = llm.invoke(prompt)
            response_text = response.content.strip()

            # Parse JSON response
            # Handle potential markdown code blocks
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0].strip()
            elif "```" in response_text:
                response_text = response_text.split("```")[1].split("```")[0].strip()

            decision = json.loads(response_text)

            next_topic = decision.get("nextSearchTopic")
            should_continue = decision.get("shouldContinue", False)
            reasoning = decision.get("reasoning", "")

            logger.info(f"Analysis decision: continue={should_continue}, next_topic={next_topic}")
            logger.debug(f"Reasoning: {reasoning}")

            return {
                "next_search_topic": next_topic,
                "should_continue": should_continue,
                "current_iteration": state['current_iteration'] + 1
            }

        except Exception as e:
            logger.error(f"Error in analyze_findings: {str(e)}", exc_info=True)
            # On error, stop iteration
            return {
                "should_continue": False,
                "next_search_topic": None
            }

    return analyze_findings


def _create_search_node():
    """
    Create the search node that executes searches across multiple sources.

    This node:
    1. Takes the next_search_topic from state
    2. Searches across 6 sources: Tavily (web), arXiv, PubMed, ChemRxiv, OpenAlex, Semantic Scholar
    3. Accumulates results in findings
    """
    def search_knowledge(state: DeepResearchState) -> dict:
        """Execute searches across all available sources."""
        try:
            next_topic = state.get('next_search_topic')
            if not next_topic:
                logger.warning("No next_search_topic provided, skipping search")
                return {}

            logger.info(f"Searching for: '{next_topic}'")

            all_results = []
            rate_limited = state.get('semantic_scholar_rate_limited', False)

            # Search Tavily (web search)
            try:
                tavily_results = search_tavily.invoke({"query": next_topic, "max_results": 3})
                all_results.append(f"## Web Search (Tavily)\n\n{tavily_results}")
                logger.debug("Tavily search completed")
            except Exception as e:
                logger.warning(f"Tavily search failed: {str(e)}")

            # Search arXiv
            try:
                arxiv_results = search_arxiv.invoke({"query": next_topic, "max_results": 3})
                all_results.append(f"## Academic Papers (arXiv)\n\n{arxiv_results}")
                logger.debug("arXiv search completed")
            except Exception as e:
                logger.warning(f"arXiv search failed: {str(e)}")

            # Search PubMed
            try:
                pubmed_results = search_pubmed.invoke({"query": next_topic, "max_results": 3})
                all_results.append(f"## Biomedical Literature (PubMed)\n\n{pubmed_results}")
                logger.debug("PubMed search completed")
            except Exception as e:
                logger.warning(f"PubMed search failed: {str(e)}")

            # Search ChemRxiv (chemistry preprints)
            try:
                chemrxiv_results = search_chemrxiv.invoke({"query": next_topic, "max_results": 3})
                all_results.append(f"## Chemistry Papers (ChemRxiv)\n\n{chemrxiv_results}")
                logger.debug("ChemRxiv search completed")
            except Exception as e:
                logger.warning(f"ChemRxiv search failed: {str(e)}")

            # Search OpenAlex (comprehensive academic database)
            try:
                openalex_results = search_openalex.invoke({"query": next_topic, "max_results": 3})
                all_results.append(f"## Academic Database (OpenAlex)\n\n{openalex_results}")
                logger.debug("OpenAlex search completed")
            except Exception as e:
                logger.warning(f"OpenAlex search failed: {str(e)}")

            # Search Semantic Scholar (AI/ML/science papers with citations)
            # Skip if already rate limited in a previous iteration
            if not rate_limited:
                try:
                    semantic_results = search_semantic_scholar.invoke({"query": next_topic, "max_results": 3})
                    all_results.append(f"## Research Papers (Semantic Scholar)\n\n{semantic_results}")
                    logger.debug("Semantic Scholar search completed")

                    # Check if rate limit was hit
                    if "rate limit" in semantic_results.lower():
                        logger.info("Semantic Scholar rate limit detected, will skip in future iterations")
                        rate_limited = True
                except Exception as e:
                    logger.warning(f"Semantic Scholar search failed: {str(e)}")
            else:
                logger.info("Skipping Semantic Scholar search (rate limited in previous iteration)")

            # Combine all results
            combined_results = "\n\n---\n\n".join(all_results)

            # Create finding entry
            finding = {
                "topic": next_topic,
                "results": combined_results
            }

            # Update state
            new_findings = state['findings'] + [finding]
            new_searched_topics = state['searched_topics'] + [next_topic]

            logger.info(f"Search completed. Total findings: {len(new_findings)}")

            return {
                "findings": new_findings,
                "searched_topics": new_searched_topics,
                "semantic_scholar_rate_limited": rate_limited
            }

        except Exception as e:
            logger.error(f"Error in search_knowledge: {str(e)}", exc_info=True)
            return {}

    return search_knowledge


def _create_synthesize_node(llm: ChatGoogleGenerativeAI):
    """
    Create the synthesis node that generates the final comprehensive report.

    This node:
    1. Reviews all accumulated findings
    2. Synthesizes them into a coherent report
    3. Includes citations and key insights
    """
    synthesize_prompt = _load_prompt("synthesize_prompt.txt")

    def synthesize_report(state: DeepResearchState) -> dict:
        """Synthesize all findings into a comprehensive report."""
        try:
            logger.info("Synthesizing final report from all findings")

            # Prepare all findings for synthesis
            findings_text = "\n\n".join([
                f"### Search {i+1}: {f['topic']}\n\n{f['results']}"
                for i, f in enumerate(state['findings'])
            ])

            # Get language name for clarity
            language_name = get_language_name(state['user_language'])

            # Create synthesis prompt
            prompt = f"""{synthesize_prompt}

## USER LANGUAGE REQUIREMENT (CRITICAL)
**You MUST write your ENTIRE report in: {language_name}**
**Language Code: {state['user_language']}**
**This is the language the user used in their original query.**

## Original Research Topic
{state['topic']}

## All Research Findings
{findings_text}

## Your Task
Create a comprehensive analysis based on the findings above. Your report should:
1. Provide key insights and conclusions
2. Cite sources appropriately
3. Identify remaining uncertainties or gaps
4. Be detailed and thorough (long-form content expected)
5. Use markdown formatting

**REMEMBER: Write the ENTIRE report in {language_name} ({state['user_language']})**

Generate the report now:
"""

            # Generate report
            response = llm.invoke(prompt)
            final_report = response.content.strip()

            logger.info(f"Final report generated ({len(final_report)} characters) in language: {state['user_language']}")

            return {
                "final_report": final_report
            }

        except Exception as e:
            logger.error(f"Error in synthesize_report: {str(e)}", exc_info=True)
            return {
                "final_report": f"Error generating report: {str(e)}"
            }

    return synthesize_report


def _should_continue_research(state: DeepResearchState) -> Literal["search", "synthesize"]:
    """
    Conditional edge function to determine next node.

    Returns:
        "search" if should continue iteration
        "synthesize" if research is complete
    """
    if state.get('should_continue', False) and state.get('next_search_topic'):
        return "search"
    else:
        return "synthesize"


def _create_supervisor_wrapper(research_graph):
    """
    Create a wrapper that makes the Deep Research graph compatible with Supervisor's message-based interface.

    This wrapper:
    1. Receives messages from Supervisor
    2. Extracts topic and depth from the latest message
    3. Invokes the Deep Research graph
    4. Returns the final_report as a message
    """
    def wrapper_node(state: MessagesState) -> dict:
        """Wrapper node that converts messages to DeepResearchState and back."""
        try:
            # Get the latest message from Supervisor
            messages = state.get("messages", [])
            if not messages:
                logger.error("No messages received from Supervisor")
                return {"messages": [AIMessage(content="Error: No research topic provided", name="deep_researcher")]}

            # LangGraphベストプラクティス: メッセージタイプでフィルタリング
            # 転送メッセージ（ToolMessage）を除外し、実際のユーザークエリ（HumanMessage）のみを取得
            human_messages = [msg for msg in messages if isinstance(msg, HumanMessage)]
            if not human_messages:
                logger.error("No HumanMessage found in message history")
                return {"messages": [AIMessage(content="Error: 有効なリサーチトピックが見つかりませんでした", name="deep_researcher")]}

            user_query = human_messages[-1].content
            logger.info(f"Deep Research received query: {user_query}")

            # Detect user's language using LLM
            api_key = get_gemini_api_key()
            model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
            llm = ChatGoogleGenerativeAI(model=model_name, api_key=api_key)
            user_language = detect_language(user_query, llm)
            logger.info(f"Detected user language: {user_language}")

            # Extract depth from query (e.g., "research topic, depth 5")
            depth_match = re.search(r'depth\s*[=:]?\s*(\d+)', user_query, re.IGNORECASE)
            depth = int(depth_match.group(1)) if depth_match else 3  # Default depth: 3

            # Remove depth specification from topic
            topic = re.sub(r',?\s*depth\s*[=:]?\s*\d+', '', user_query, flags=re.IGNORECASE).strip()

            logger.info(f"Extracted topic: '{topic}', depth: {depth}, language: {user_language}")

            # Initialize Deep Research state
            initial_state: DeepResearchState = {
                "topic": topic,
                "depth": depth,
                "user_language": user_language,
                "current_iteration": 0,
                "searched_topics": [],
                "findings": [],
                "next_search_topic": None,
                "should_continue": True,
                "semantic_scholar_rate_limited": False,
                "final_report": ""
            }

            # Invoke the Deep Research graph
            logger.info(f"Starting Deep Research with depth={depth}, language={user_language}")
            result = research_graph.invoke(initial_state, config={"recursion_limit": depth * 3 + 10})

            # Extract final report
            final_report = result.get("final_report", "")

            if not final_report:
                final_report = "Research completed but no report was generated."

            logger.info(f"Deep Researcher completed. Report length: {len(final_report)} characters")

            # Return as AIMessage
            return {"messages": [AIMessage(content=final_report, name="deep_researcher")]}

        except Exception as e:
            logger.error(f"Error in Deep Researcher wrapper: {str(e)}", exc_info=True)
            error_message = f"Error during deep research: {str(e)}"
            return {"messages": [AIMessage(content=error_message, name="deep_researcher")]}

    # Build a simple graph with the wrapper node
    wrapper_builder = StateGraph(MessagesState)
    wrapper_builder.add_node("deep_researcher", wrapper_node)
    wrapper_builder.add_edge(START, "deep_researcher")
    wrapper_builder.add_edge("deep_researcher", END)

    return wrapper_builder.compile(name="deep_researcher")


def create_deep_researcher():
    """
    Create a Deep Researcher agent compatible with LangGraph Supervisor.

    This function creates the core Deep Researcher graph and wraps it in a
    message-based interface that the Supervisor can interact with.

    The agent accepts queries in the format:
    - "Research [topic]" → Uses default depth (3)
    - "Research [topic], depth 5" → Uses specified depth

    Returns:
        A compiled LangGraph application compatible with Supervisor workflows

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Create the core Deep Researcher graph
    core_graph = _create_core_deep_researcher_graph()

    # Wrap it for Supervisor compatibility
    supervisor_compatible_graph = _create_supervisor_wrapper(core_graph)

    logger.info("Deep Researcher with Supervisor wrapper created successfully")
    return supervisor_compatible_graph


def _create_core_deep_researcher_graph():
    """
    Create the core Deep Researcher graph (internal implementation).

    This is the original create_deep_researcher function, now renamed
    to distinguish it from the public API.
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Deep Researcher.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    # Initialize LLM
    model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
    llm = ChatGoogleGenerativeAI(
        model=model_name,
        api_key=api_key,
        temperature=0.7  # Slightly higher temperature for creative research
    )
    logger.info(f"Initialized ChatGoogleGenerativeAI with {model_name} for Deep Researcher")

    # Create nodes
    analyze_node = _create_analyze_node(llm)
    search_node = _create_search_node()
    synthesize_node = _create_synthesize_node(llm)

    # Build the graph
    builder = StateGraph(DeepResearchState)

    # Add nodes
    builder.add_node("analyze", analyze_node)
    builder.add_node("search", search_node)
    builder.add_node("synthesize", synthesize_node)

    # Add edges
    builder.add_edge(START, "analyze")
    builder.add_conditional_edges(
        "analyze",
        _should_continue_research,
        {
            "search": "search",
            "synthesize": "synthesize"
        }
    )
    builder.add_edge("search", "analyze")  # Loop back to analyze after search
    builder.add_edge("synthesize", END)

    # Compile the graph
    graph = builder.compile()

    logger.info("Core Deep Researcher graph created successfully")
    return graph
