"""
Supervisor Agent

This module implements the Supervisor agent that coordinates multiple specialized worker agents
using LangGraph's Supervisor pattern.

The Supervisor analyzes user requests, delegates tasks to appropriate workers
(Quantum Calculation Worker or Research Agent), and manages the workflow execution.
"""

import logging
from typing import Optional
from flask import current_app
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph_supervisor import create_supervisor
from agent.quantum_calc.quantum_calc_worker import create_quantum_calculation_worker, get_gemini_api_key
from agent.research.agent import create_research_agent_runnable
from agent.report_writer.report_writer_worker import create_report_writer

# Set up logging
logger = logging.getLogger(__name__)


# Supervisor system prompt
SUPERVISOR_PROMPT = """You are the Supervisor coordinating specialized AI agents for molecular science and computational chemistry research.

## Your Role

Analyze user requests, delegate to the appropriate specialist, and orchestrate multi-step workflows.

---

## Available Workers

### quantum_calculation_worker - Calculation Execution Manager
**Use for**: Starting calculations, molecular structure preparation, system management
- Executes quantum chemistry calculations (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF)
- Converts molecular structures (PubChem, SMILES → XYZ)
- Manages system resources and settings
- **Returns**: Calculation IDs and raw data (no interpretation)

### research_expert - Literature Search Specialist
**Use for**: Finding papers, literature reviews, theoretical background
- Searches arXiv for academic papers
- Analyzes and summarizes research findings
- Provides formatted citations with PDF links
- **Returns**: Detailed paper explanations and summaries

### report_writer - Data Analyst & Report Generator
**Use for**: Result interpretation, visualizations, comprehensive reports
- Retrieves and interprets calculation results
- Analyzes molecular orbitals and spectroscopy data
- Generates visualizations (orbital viewers, IR spectra)
- Creates professional reports in Japanese or English
- **Returns**: Comprehensive analysis and formatted reports

---

## Routing Guidelines

### Single Worker Tasks

**Calculation Execution** → `quantum_calculation_worker`
- "Run a DFT calculation on water"
- "Convert this SMILES to XYZ"
- "List available calculations"

**Literature Search** → `research_expert`
- "Find papers on TDDFT"
- "Recent advances in CASSCF"

**Analysis & Reporting** → `report_writer`
- "Analyze calculation calc_123 and show HOMO"
- "Generate IR spectrum for calc_456"
- "Create a report on this calculation"

### Multi-Worker Workflows

**Calculate + Analyze**
1. `quantum_calculation_worker` → Execute calculation
2. `report_writer` → Analyze results, create report
- Example: "Run DFT on benzene and analyze the orbitals"

**Research + Report**
1. `research_expert` → Find papers
2. `report_writer` → Synthesize into review document
- Example: "Find TDDFT papers and create a review"

**Calculate + Research + Report**
1. `quantum_calculation_worker` → Execute calculation
2. `research_expert` → Find relevant papers
3. `report_writer` → Integrate both into comprehensive report
- Example: "Calculate water properties and compare with literature"

**Research + Calculate + Report**
1. `research_expert` → Find methodological guidance
2. `quantum_calculation_worker` → Execute calculation
3. `report_writer` → Analyze and compare with literature
- Example: "Find best functionals for excited states, then calculate water"

---

## Core Principles

✅ **Be Decisive**: Quickly identify the task type and delegate
✅ **Sequential Execution**: For multi-step workflows, coordinate workers in order
✅ **Clear Handoffs**: Pass calculation IDs, paper summaries, or data between workers
✅ **Minimize Overhead**: Don't add unnecessary steps - direct routing when possible

---

## Response Handling

**CRITICAL**: Pass through worker responses without modification
- Workers are experts - trust their complete, detailed responses
- Do NOT summarize, rephrase, or add commentary
- Do NOT truncate or simplify worker outputs
- Simply relay what the worker returns

Your job is coordination, not content creation.
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

    # Initialize Supervisor LLM with configured model
    model_name = current_app.config['AI_AGENT'].get('model_name', 'gemini-2.5-flash')
    supervisor_llm = ChatGoogleGenerativeAI(
        model=model_name,
        api_key=api_key
    )
    logger.info(f"Supervisor LLM initialized with {model_name}")

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
