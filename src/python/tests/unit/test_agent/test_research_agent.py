"""
Unit tests for Research Agent.

Tests the Research Agent implementation with mocked LangChain and Google Gemini
components to verify runnable creation, prompt configuration, and tool binding.
"""

import pytest
from unittest.mock import MagicMock, patch, PropertyMock

from research.agent import create_research_agent_runnable, RESEARCH_AGENT_PROMPT


# ============================================================================
# Research Agent Runnable Creation Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_runnable_returns_runnable(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN all components are properly initialized
    THEN it should return a compiled LangGraph agent
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    assert result is not None
    assert result == mock_agent
    mock_llm_class.assert_called_once()
    mock_create_react_agent.assert_called_once()


@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.get_gemini_api_key')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_uses_correct_model(
    mock_llm_class,
    mock_get_api_key,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN LLM is initialized
    THEN it should use gemini-2.5-flash model
    """
    # ARRANGE
    mock_get_api_key.return_value = "test-api-key"
    
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    create_research_agent_runnable()

    # ASSERT
    mock_llm_class.assert_called_once_with(model="gemini-2.5-flash", api_key="test-api-key")


# ============================================================================
# ReAct Agent Creation Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_calls_create_react_agent(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN creating the agent
    THEN it should call create_react_agent with correct parameters
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    create_research_agent_runnable()

    # ASSERT
    mock_create_react_agent.assert_called_once()
    call_args = mock_create_react_agent.call_args
    
    # Verify LLM is passed
    assert call_args[0][0] == mock_llm
    # Verify tools parameter
    assert "tools" in call_args[1]
    # Verify state_modifier parameter
    assert "state_modifier" in call_args[1]
    assert call_args[1]["state_modifier"] == RESEARCH_AGENT_PROMPT


@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_uses_research_system_prompt(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN agent is created with state_modifier
    THEN it should use RESEARCH_AGENT_PROMPT as system prompt
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    create_research_agent_runnable()

    # ASSERT
    call_args = mock_create_react_agent.call_args
    assert call_args[1]["state_modifier"] == RESEARCH_AGENT_PROMPT


# ============================================================================
# Tool Binding Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.search_arxiv')
def test_create_research_agent_binds_search_tool(
    mock_search_arxiv,
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN tools are passed to create_react_agent
    THEN it should include search_arxiv tool
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    create_research_agent_runnable()

    # ASSERT
    call_args = mock_create_react_agent.call_args
    tools_list = call_args[1]["tools"]
    assert len(tools_list) == 1
    assert tools_list[0] == mock_search_arxiv


# ============================================================================
# System Prompt Content Tests
# ============================================================================

def test_research_agent_prompt_mentions_arxiv():
    """
    GIVEN RESEARCH_AGENT_PROMPT constant
    WHEN checking its content
    THEN it should mention arXiv as a data source
    """
    # ASSERT
    assert "arxiv" in RESEARCH_AGENT_PROMPT.lower() or "arXiv" in RESEARCH_AGENT_PROMPT


def test_research_agent_prompt_mentions_research_role():
    """
    GIVEN RESEARCH_AGENT_PROMPT constant
    WHEN checking its content
    THEN it should describe research assistant role
    """
    # ASSERT
    prompt_lower = RESEARCH_AGENT_PROMPT.lower()
    assert "research" in prompt_lower or "literature" in prompt_lower


def test_research_agent_prompt_mentions_quantum_chemistry():
    """
    GIVEN RESEARCH_AGENT_PROMPT constant
    WHEN checking its content
    THEN it should mention quantum chemistry focus areas
    """
    # ASSERT
    prompt_lower = RESEARCH_AGENT_PROMPT.lower()
    assert "quantum chemistry" in prompt_lower or "chemistry" in prompt_lower


def test_research_agent_prompt_provides_usage_instructions():
    """
    GIVEN RESEARCH_AGENT_PROMPT constant
    WHEN checking its content
    THEN it should provide clear instructions for tool usage
    """
    # ASSERT
    prompt_lower = RESEARCH_AGENT_PROMPT.lower()
    # Should mention using the search tool
    assert "search" in prompt_lower or "search_arxiv" in prompt_lower


def test_research_agent_prompt_mentions_markdown_formatting():
    """
    GIVEN RESEARCH_AGENT_PROMPT constant
    WHEN checking its content
    THEN it should instruct agent to use Markdown formatting
    """
    # ASSERT
    assert "Markdown" in RESEARCH_AGENT_PROMPT or "markdown" in RESEARCH_AGENT_PROMPT


def test_research_agent_prompt_is_not_empty():
    """
    GIVEN RESEARCH_AGENT_PROMPT constant
    WHEN checking its length
    THEN it should contain substantial instructions
    """
    # ASSERT
    assert len(RESEARCH_AGENT_PROMPT) > 100
    assert len(RESEARCH_AGENT_PROMPT.strip()) > 0


# ============================================================================
# Agent Construction Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_returns_agent(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN agent is created
    THEN it should return the compiled agent
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    assert result == mock_agent


# ============================================================================
# Integration Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_full_pipeline(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN all components are initialized in sequence
    THEN complete pipeline should be set up correctly
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_agent = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    # Verify call sequence
    assert mock_llm_class.called
    assert mock_create_react_agent.called

    # Verify final result
    assert result is not None
    assert result == mock_agent


# ============================================================================
# Error Handling Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_handles_llm_initialization_error(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN ChatGoogleGenerativeAI raises an error during initialization
    WHEN create_research_agent_runnable is called
    THEN it should propagate the error
    """
    # ARRANGE
    mock_llm_class.side_effect = ValueError("Invalid model configuration")

    # ACT & ASSERT
    with pytest.raises(ValueError, match="Invalid model configuration"):
        create_research_agent_runnable()


@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_handles_agent_creation_error(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_react_agent raises an error during creation
    WHEN create_research_agent_runnable is called
    THEN it should propagate the error
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm
    
    mock_create_react_agent.side_effect = RuntimeError(
        "Invalid agent configuration"
    )

    # ACT & ASSERT
    with pytest.raises(RuntimeError, match="Invalid agent configuration"):
        create_research_agent_runnable()


# ============================================================================
# Agent Interface Tests
# ============================================================================

@patch('langgraph.prebuilt.create_react_agent')
@patch('research.agent.ChatGoogleGenerativeAI')
def test_create_research_agent_returns_invokable_object(
    mock_llm_class,
    mock_create_react_agent
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN the returned object is inspected
    THEN it should support LangGraph agent interface
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_llm_class.return_value = mock_llm

    # Create an agent mock with invoke method
    mock_agent = MagicMock()
    mock_agent.invoke = MagicMock()
    mock_create_react_agent.return_value = mock_agent

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    # Verify the result has agent interface methods
    assert hasattr(result, 'invoke')
