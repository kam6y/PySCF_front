"""
Supervisor Agent

This module implements the Supervisor agent that coordinates multiple specialized worker agents
using LangGraph's Supervisor pattern.

The Supervisor analyzes user requests, delegates tasks to appropriate workers
(Quantum Calculation Worker or Research Agent), and manages the workflow execution.
"""

import logging
from typing import Optional
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph_supervisor import create_supervisor
from agent.quantum_calc.quantum_calc_worker import create_quantum_calculation_worker, get_gemini_api_key
from agent.research.agent import create_research_agent_runnable
from agent.report_writer.report_writer_worker import create_report_writer

# Set up logging
logger = logging.getLogger(__name__)


# Supervisor system prompt
SUPERVISOR_PROMPT = """You are an intelligent Supervisor managing a team of specialized AI agents for molecular science and research.

Your role is to:
1. Analyze user requests and understand their intent
2. Delegate tasks to the most appropriate specialist worker
3. Coordinate complex workflows that may require multiple workers
4. Ensure tasks are completed effectively

**Available Workers:**

**quantum_calculation_worker** - Quantum Chemistry and Molecular Analysis Expert
- Handles: Quantum chemistry calculations (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF)
- Capabilities:
  * Running and managing quantum chemistry calculations
  * Molecular structure analysis and validation
  * PubChem searches and SMILES conversions
  * Molecular orbital analysis and visualization
  * Geometry optimization and frequency analysis
  * IR spectrum generation
  * System settings and resource management
- Use for: Any computational chemistry task, molecular calculations, orbital analysis, structure optimization

**research_expert** - Academic Research and Literature Search Expert
- Handles: Academic paper searches and literature reviews
- Capabilities:
  * Searching arXiv for scientific papers
  * Finding papers on quantum chemistry, computational chemistry, DFT, etc.
  * Summarizing research findings
  * Providing formatted citations with PDF links
- Use for: Literature searches, finding papers, research summaries, academic references

**report_writer** - Scientific Report Generation Expert
- Handles: Creating comprehensive scientific reports and documentation
- Capabilities:
  * Generating structured calculation reports from quantum chemistry results
  * Creating literature review documents from research findings
  * Synthesizing information from multiple sources into coherent reports
  * Formatting scientific documents with proper structure, tables, and citations
  * Writing in Japanese or English based on user preference
  * Producing executive summaries and technical documentation
- Use for: Report generation, documentation creation, summary writing, formal scientific reports

**Decision Guidelines:**

For QUANTUM CALCULATIONS and MOLECULAR TASKS → Use **quantum_calculation_worker**:
- "Run a DFT calculation on water"
- "Show me the HOMO-LUMO gap of benzene"
- "Optimize the geometry of methane"
- "What calculations are available?"
- "Generate an IR spectrum"
- "Search PubChem for aspirin"
- "Convert this SMILES to XYZ"

For LITERATURE SEARCH and ACADEMIC PAPERS → Use **research_expert**:
- "Find papers about time-dependent DFT"
- "What are recent advances in CASSCF?"
- "Show me research on excited states"
- "Literature review on quantum chemistry methods"

For REPORT GENERATION and DOCUMENTATION → Use **report_writer**:
- "Create a report on this calculation"
- "Write a summary of these findings"
- "Generate a comprehensive report"
- "Prepare documentation for this analysis"
- "Create a literature review report"
- "Write a technical summary"

For COMPLEX WORKFLOWS → Delegate sequentially:
1. **Calculation + Report**: Use **quantum_calculation_worker** first, then **report_writer** to document results
2. **Research + Report**: Use **research_expert** first, then **report_writer** to create literature review
3. **Calculation + Research + Report**: Use **quantum_calculation_worker**, then **research_expert** for context, finally **report_writer** to integrate findings
4. **Research + Calculation**: Use **research_expert** to gather literature context, then **quantum_calculation_worker** to perform calculations based on findings

**Important:**
- Always choose the most specialized worker for the task
- Be decisive - analyze the request and delegate immediately
- Provide clear task descriptions when delegating
- If a task requires multiple workers, coordinate them sequentially
- Trust your workers - they are experts in their domains

**Critical - Response Handling:**
- **NEVER modify, summarize, or rephrase worker responses**
- **ALWAYS pass through worker responses exactly as received**
- Workers provide complete, detailed answers - do not alter them
- Your only job after delegation is to relay the worker's response unchanged
- If the worker provides a detailed explanation with multiple papers, return ALL of it
- Do not add your own commentary or interpretation to worker responses
"""


def create_supervisor_agent():
    """
    Create a Supervisor agent using LangGraph's create_supervisor.

    The Supervisor coordinates between:
    - Quantum Calculation Worker (quantum chemistry and molecular analysis)
    - Research Agent (academic literature search)
    - Report Writer (scientific report generation)

    Returns:
        A compiled LangGraph Supervisor workflow ready for execution

    Raises:
        ValueError: If Gemini API key is not configured
    """
    # Get API key
    api_key = get_gemini_api_key()
    if not api_key:
        logger.error("Gemini API key not found. Cannot initialize Supervisor.")
        raise ValueError(
            "Gemini API key not configured. Please set GEMINI_API_KEY environment variable "
            "or configure it in application settings."
        )

    logger.info("Initializing Supervisor agent...")

    # Create worker agents
    try:
        quantum_worker = create_quantum_calculation_worker()
        logger.info("Quantum Calculation Worker initialized")
    except Exception as e:
        logger.error(f"Failed to create Quantum Calculation Worker: {e}", exc_info=True)
        raise

    try:
        research_worker = create_research_agent_runnable()
        logger.info("Research Agent initialized")
    except Exception as e:
        logger.error(f"Failed to create Research Agent: {e}", exc_info=True)
        raise

    try:
        report_writer_worker = create_report_writer()
        logger.info("Report Writer initialized")
    except Exception as e:
        logger.error(f"Failed to create Report Writer: {e}", exc_info=True)
        raise

    # Initialize Supervisor LLM (using more capable model for coordination)
    supervisor_llm = ChatGoogleGenerativeAI(
        model="gemini-2.5-flash",
        api_key=api_key
    )
    logger.info("Supervisor LLM initialized with gemini-2.5-flash")

    # Create supervisor workflow using create_supervisor
    workflow = create_supervisor(
        [quantum_worker, research_worker, report_writer_worker],
        model=supervisor_llm,
        prompt=SUPERVISOR_PROMPT,
    )

    logger.info("Supervisor workflow created successfully")
    return workflow


def get_compiled_supervisor():
    """
    Get a compiled Supervisor agent ready for execution.

    This is the main entry point for the API to obtain the supervisor graph executor.

    Returns:
        Compiled LangGraph Supervisor application
    """
    workflow = create_supervisor_agent()
    app = workflow.compile()
    logger.info("Supervisor compiled and ready for execution")
    return app
