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

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_runnable_returns_runnable(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN all components are properly initialized
    THEN it should return a LangChain runnable chain
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm_with_tools = MagicMock()
    mock_llm.bind_tools.return_value = mock_llm_with_tools
    mock_llm_class.return_value = mock_llm

    # Mock the | operator for chaining
    mock_chain = MagicMock()
    mock_prompt.__or__ = MagicMock(return_value=mock_chain)

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    assert result is not None
    assert result == mock_chain
    mock_prompt_template_class.from_messages.assert_called_once()
    mock_llm_class.assert_called_once()
    mock_llm.bind_tools.assert_called_once()


@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_uses_correct_model(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN LLM is initialized
    THEN it should use gemini-1.5-flash model
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm.bind_tools.return_value = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_prompt.__or__ = MagicMock(return_value=MagicMock())

    # ACT
    create_research_agent_runnable()

    # ASSERT
    mock_llm_class.assert_called_once_with(model="gemini-1.5-flash")


# ============================================================================
# Prompt Template Tests
# ============================================================================

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_configures_prompt_template(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN prompt template is created
    THEN it should use system and human message templates
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm.bind_tools.return_value = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_prompt.__or__ = MagicMock(return_value=MagicMock())

    # ACT
    create_research_agent_runnable()

    # ASSERT
    mock_prompt_template_class.from_messages.assert_called_once()
    call_args = mock_prompt_template_class.from_messages.call_args[0][0]

    # Verify message structure
    assert len(call_args) == 2
    assert call_args[0][0] == "system"
    assert call_args[0][1] == RESEARCH_AGENT_PROMPT
    assert call_args[1][0] == "human"
    assert call_args[1][1] == "{input}"


@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_uses_research_system_prompt(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN system prompt is set
    THEN it should use RESEARCH_AGENT_PROMPT constant
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm.bind_tools.return_value = MagicMock()
    mock_llm_class.return_value = mock_llm

    mock_prompt.__or__ = MagicMock(return_value=MagicMock())

    # ACT
    create_research_agent_runnable()

    # ASSERT
    call_args = mock_prompt_template_class.from_messages.call_args[0][0]
    system_message = call_args[0]
    assert system_message[1] == RESEARCH_AGENT_PROMPT


# ============================================================================
# Tool Binding Tests
# ============================================================================

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
@patch('research.agent.search_arxiv')
def test_create_research_agent_binds_search_tool(
    mock_search_arxiv,
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN tools are bound to LLM
    THEN it should bind search_arxiv tool
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm_with_tools = MagicMock()
    mock_llm.bind_tools.return_value = mock_llm_with_tools
    mock_llm_class.return_value = mock_llm

    mock_prompt.__or__ = MagicMock(return_value=MagicMock())

    # ACT
    create_research_agent_runnable()

    # ASSERT
    mock_llm.bind_tools.assert_called_once()
    call_args = mock_llm.bind_tools.call_args[0][0]
    assert len(call_args) == 1
    assert call_args[0] == mock_search_arxiv


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
# Chain Construction Tests
# ============================================================================

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_chains_prompt_and_llm(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN runnable chain is created
    THEN it should chain prompt template with tool-bound LLM
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm_with_tools = MagicMock()
    mock_llm.bind_tools.return_value = mock_llm_with_tools
    mock_llm_class.return_value = mock_llm

    mock_chain = MagicMock()
    mock_prompt.__or__ = MagicMock(return_value=mock_chain)

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    # Verify that prompt | llm_with_tools was created
    mock_prompt.__or__.assert_called_once_with(mock_llm_with_tools)
    assert result == mock_chain


# ============================================================================
# Integration Tests
# ============================================================================

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_full_pipeline(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN all components are initialized in sequence
    THEN complete pipeline should be set up correctly
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm_with_tools = MagicMock()
    mock_llm.bind_tools.return_value = mock_llm_with_tools
    mock_llm_class.return_value = mock_llm

    mock_chain = MagicMock()
    mock_prompt.__or__ = MagicMock(return_value=mock_chain)

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    # Verify call sequence
    assert mock_prompt_template_class.from_messages.called
    assert mock_llm_class.called
    assert mock_llm.bind_tools.called
    assert mock_prompt.__or__.called

    # Verify final result
    assert result is not None
    assert result == mock_chain


# ============================================================================
# Error Handling Tests
# ============================================================================

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_handles_llm_initialization_error(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN ChatGoogleGenerativeAI raises an error during initialization
    WHEN create_research_agent_runnable is called
    THEN it should propagate the error
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm_class.side_effect = ValueError("Invalid model configuration")

    # ACT & ASSERT
    with pytest.raises(ValueError, match="Invalid model configuration"):
        create_research_agent_runnable()


@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_handles_prompt_creation_error(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN ChatPromptTemplate raises an error during creation
    WHEN create_research_agent_runnable is called
    THEN it should propagate the error
    """
    # ARRANGE
    mock_prompt_template_class.from_messages.side_effect = RuntimeError(
        "Invalid prompt template"
    )

    # ACT & ASSERT
    with pytest.raises(RuntimeError, match="Invalid prompt template"):
        create_research_agent_runnable()


# ============================================================================
# Runnable Interface Tests
# ============================================================================

@patch('research.agent.ChatGoogleGenerativeAI')
@patch('research.agent.ChatPromptTemplate')
def test_create_research_agent_returns_chainable_object(
    mock_prompt_template_class,
    mock_llm_class
):
    """
    GIVEN create_research_agent_runnable is called
    WHEN the returned object is inspected
    THEN it should support LangChain runnable interface
    """
    # ARRANGE
    mock_prompt = MagicMock()
    mock_prompt_template_class.from_messages.return_value = mock_prompt

    mock_llm = MagicMock()
    mock_llm_with_tools = MagicMock()
    mock_llm.bind_tools.return_value = mock_llm_with_tools
    mock_llm_class.return_value = mock_llm

    # Create a chain mock with runnable methods
    mock_chain = MagicMock()
    mock_chain.invoke = MagicMock()
    mock_chain.stream = MagicMock()
    mock_prompt.__or__ = MagicMock(return_value=mock_chain)

    # ACT
    result = create_research_agent_runnable()

    # ASSERT
    # Verify the result has runnable interface methods
    assert hasattr(result, 'invoke')
    assert hasattr(result, 'stream')
