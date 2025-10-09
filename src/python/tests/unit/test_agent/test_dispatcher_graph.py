"""
Unit tests for LangGraph Multi-Agent Dispatcher.

Tests the dispatcher graph implementation including router node, agent nodes,
message conversion utilities, and full graph execution.
"""

import pytest
from unittest.mock import MagicMock, patch, call
from langchain_core.messages import HumanMessage, AIMessage

from agent.graph import (
    AgentState,
    router_node,
    molecular_agent_node,
    research_agent_node,
    create_dispatcher_graph,
    _convert_langchain_to_dict_format,
)


# ============================================================================
# Message Conversion Tests
# ============================================================================

def test_convert_langchain_to_dict_format_single_user_message():
    """
    GIVEN a list with single HumanMessage
    WHEN _convert_langchain_to_dict_format is called
    THEN it should return correctly formatted dict
    """
    # ARRANGE
    messages = [HumanMessage(content="Test query")]

    # ACT
    result = _convert_langchain_to_dict_format(messages)

    # ASSERT
    assert len(result) == 1
    assert result[0]["role"] == "user"
    assert result[0]["parts"] == [{"text": "Test query"}]


def test_convert_langchain_to_dict_format_conversation():
    """
    GIVEN a conversation with user and AI messages
    WHEN _convert_langchain_to_dict_format is called
    THEN it should preserve conversation structure
    """
    # ARRANGE
    messages = [
        HumanMessage(content="Hello"),
        AIMessage(content="Hi there!"),
        HumanMessage(content="How are you?"),
        AIMessage(content="I'm doing well!")
    ]

    # ACT
    result = _convert_langchain_to_dict_format(messages)

    # ASSERT
    assert len(result) == 4
    assert result[0]["role"] == "user"
    assert result[0]["parts"][0]["text"] == "Hello"
    assert result[1]["role"] == "model"
    assert result[1]["parts"][0]["text"] == "Hi there!"
    assert result[2]["role"] == "user"
    assert result[3]["role"] == "model"


def test_convert_langchain_to_dict_format_empty_list():
    """
    GIVEN an empty message list
    WHEN _convert_langchain_to_dict_format is called
    THEN it should return empty list
    """
    # ARRANGE
    messages = []

    # ACT
    result = _convert_langchain_to_dict_format(messages)

    # ASSERT
    assert result == []


# ============================================================================
# Router Node Tests
# ============================================================================

@patch('agent.graph.ChatGoogleGenerativeAI')
def test_router_node_routes_to_molecular(mock_llm_class):
    """
    GIVEN a molecular chemistry query
    WHEN router_node is invoked
    THEN it should set next_node to 'molecular'
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_structured_llm = MagicMock()
    mock_llm.with_structured_output.return_value = mock_structured_llm

    # Mock the router decision
    mock_result = MagicMock()
    mock_result.destination = "molecular"
    mock_structured_llm.invoke.return_value = mock_result

    mock_llm_class.return_value = mock_llm

    state = {
        "messages": [HumanMessage(content="Calculate HOMO energy for water")],
        "next_node": ""
    }

    # ACT
    result = router_node(state)

    # ASSERT
    assert result["next_node"] == "molecular"
    mock_llm.with_structured_output.assert_called_once()
    mock_structured_llm.invoke.assert_called_once()


@patch('agent.graph.ChatGoogleGenerativeAI')
def test_router_node_routes_to_research(mock_llm_class):
    """
    GIVEN a research/paper search query
    WHEN router_node is invoked
    THEN it should set next_node to 'research'
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_structured_llm = MagicMock()
    mock_llm.with_structured_output.return_value = mock_structured_llm

    # Mock the router decision
    mock_result = MagicMock()
    mock_result.destination = "research"
    mock_structured_llm.invoke.return_value = mock_result

    mock_llm_class.return_value = mock_llm

    state = {
        "messages": [HumanMessage(content="Find papers about DFT")],
        "next_node": ""
    }

    # ACT
    result = router_node(state)

    # ASSERT
    assert result["next_node"] == "research"


@patch('agent.graph.ChatGoogleGenerativeAI')
def test_router_node_passes_correct_prompt(mock_llm_class):
    """
    GIVEN router_node is called
    WHEN LLM is invoked
    THEN it should receive classification instructions
    """
    # ARRANGE
    mock_llm = MagicMock()
    mock_structured_llm = MagicMock()
    mock_llm.with_structured_output.return_value = mock_structured_llm

    mock_result = MagicMock()
    mock_result.destination = "molecular"
    mock_structured_llm.invoke.return_value = mock_result

    mock_llm_class.return_value = mock_llm

    state = {
        "messages": [HumanMessage(content="Test query")],
        "next_node": ""
    }

    # ACT
    router_node(state)

    # ASSERT
    invoke_args = mock_structured_llm.invoke.call_args[0][0]
    assert "Test query" in invoke_args
    assert "molecular" in invoke_args
    assert "research" in invoke_args


# ============================================================================
# Molecular Agent Node Tests
# ============================================================================

@patch('agent.molecular_agent.MolecularAgent')
def test_molecular_agent_node_success(mock_agent_class):
    """
    GIVEN molecular_agent_node is called with valid state
    WHEN MolecularAgent returns response chunks
    THEN it should collect chunks and add AIMessage to state
    """
    # ARRANGE
    mock_agent = MagicMock()
    mock_agent.chat.return_value = iter(["Response ", "chunk ", "1"])
    mock_agent_class.return_value = mock_agent

    state = {
        "messages": [HumanMessage(content="Calculate HOMO for H2O")],
        "next_node": ""
    }

    # ACT
    result = molecular_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert isinstance(result["messages"][1], AIMessage)
    assert result["messages"][1].content == "Response chunk 1"
    mock_agent.chat.assert_called_once()


@patch('agent.molecular_agent.MolecularAgent')
def test_molecular_agent_node_converts_history_correctly(mock_agent_class):
    """
    GIVEN state with conversation history
    WHEN molecular_agent_node is called
    THEN it should convert history to dict format
    """
    # ARRANGE
    mock_agent = MagicMock()
    mock_agent.chat.return_value = iter(["Response"])
    mock_agent_class.return_value = mock_agent

    state = {
        "messages": [
            HumanMessage(content="Hello"),
            AIMessage(content="Hi"),
            HumanMessage(content="Calculate energy")
        ],
        "next_node": ""
    }

    # ACT
    molecular_agent_node(state)

    # ASSERT
    call_args = mock_agent.chat.call_args
    message_arg = call_args[0][0]
    history_arg = call_args[0][1]

    assert message_arg == "Calculate energy"
    assert len(history_arg) == 2
    assert history_arg[0]["role"] == "user"
    assert history_arg[0]["parts"][0]["text"] == "Hello"
    assert history_arg[1]["role"] == "model"
    assert history_arg[1]["parts"][0]["text"] == "Hi"


@patch('agent.molecular_agent.MolecularAgent')
def test_molecular_agent_node_handles_error(mock_agent_class):
    """
    GIVEN MolecularAgent raises an exception
    WHEN molecular_agent_node is called
    THEN it should add error message to state without crashing
    """
    # ARRANGE
    mock_agent = MagicMock()
    mock_agent.chat.side_effect = Exception("API error")
    mock_agent_class.return_value = mock_agent

    state = {
        "messages": [HumanMessage(content="Test")],
        "next_node": ""
    }

    # ACT
    result = molecular_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert isinstance(result["messages"][1], AIMessage)
    assert "error occurred" in result["messages"][1].content
    assert "API error" in result["messages"][1].content


@patch('agent.molecular_agent.MolecularAgent')
def test_molecular_agent_node_handles_empty_response(mock_agent_class):
    """
    GIVEN MolecularAgent returns empty string
    WHEN molecular_agent_node is called
    THEN it should add fallback message
    """
    # ARRANGE
    mock_agent = MagicMock()
    mock_agent.chat.return_value = iter(["", "  "])
    mock_agent_class.return_value = mock_agent

    state = {
        "messages": [HumanMessage(content="Test")],
        "next_node": ""
    }

    # ACT
    result = molecular_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert "Failed to generate a response" in result["messages"][1].content


# ============================================================================
# Research Agent Node Tests
# ============================================================================

@patch('agent.research.agent.create_research_agent_runnable')
def test_research_agent_node_success(mock_create_runnable):
    """
    GIVEN research_agent_node is called with valid state
    WHEN ResearchAgent returns response
    THEN it should add AIMessage to state
    """
    # ARRANGE
    mock_runnable = MagicMock()
    mock_response = MagicMock()
    mock_response.content = "Here are the top 3 papers..."
    mock_runnable.invoke.return_value = mock_response
    mock_create_runnable.return_value = mock_runnable

    state = {
        "messages": [HumanMessage(content="Find papers about CASSCF")],
        "next_node": ""
    }

    # ACT
    result = research_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert isinstance(result["messages"][1], AIMessage)
    assert result["messages"][1].content == "Here are the top 3 papers..."
    mock_runnable.invoke.assert_called_once_with({"input": "Find papers about CASSCF"})


@patch('agent.research.agent.create_research_agent_runnable')
def test_research_agent_node_handles_string_response(mock_create_runnable):
    """
    GIVEN ResearchAgent returns a plain string
    WHEN research_agent_node is called
    THEN it should handle the string response
    """
    # ARRANGE
    mock_runnable = MagicMock()
    mock_runnable.invoke.return_value = "String response"
    mock_create_runnable.return_value = mock_runnable

    state = {
        "messages": [HumanMessage(content="Test query")],
        "next_node": ""
    }

    # ACT
    result = research_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert result["messages"][1].content == "String response"


@patch('agent.research.agent.create_research_agent_runnable')
def test_research_agent_node_handles_error(mock_create_runnable):
    """
    GIVEN ResearchAgent raises an exception
    WHEN research_agent_node is called
    THEN it should add error message to state without crashing
    """
    # ARRANGE
    mock_runnable = MagicMock()
    mock_runnable.invoke.side_effect = Exception("arXiv connection error")
    mock_create_runnable.return_value = mock_runnable

    state = {
        "messages": [HumanMessage(content="Test")],
        "next_node": ""
    }

    # ACT
    result = research_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert isinstance(result["messages"][1], AIMessage)
    assert "error occurred" in result["messages"][1].content
    assert "arXiv connection error" in result["messages"][1].content


@patch('agent.research.agent.create_research_agent_runnable')
def test_research_agent_node_handles_empty_response(mock_create_runnable):
    """
    GIVEN ResearchAgent returns empty content
    WHEN research_agent_node is called
    THEN it should add fallback message
    """
    # ARRANGE
    mock_runnable = MagicMock()
    mock_response = MagicMock()
    mock_response.content = "   "
    mock_runnable.invoke.return_value = mock_response
    mock_create_runnable.return_value = mock_runnable

    state = {
        "messages": [HumanMessage(content="Test")],
        "next_node": ""
    }

    # ACT
    result = research_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 2
    assert "No papers found" in result["messages"][1].content


# ============================================================================
# Graph Construction Tests
# ============================================================================

def test_create_dispatcher_graph_returns_compiled_graph():
    """
    GIVEN create_dispatcher_graph is called
    WHEN graph is created
    THEN it should return a compiled LangGraph application
    """
    # ACT
    graph = create_dispatcher_graph()

    # ASSERT
    assert graph is not None
    # Compiled graph should have invoke and stream methods
    assert hasattr(graph, 'invoke')
    assert hasattr(graph, 'stream')


# ============================================================================
# End-to-End Graph Execution Tests
# ============================================================================

@patch('agent.molecular_agent.MolecularAgent')
@patch('agent.graph.ChatGoogleGenerativeAI')
def test_graph_execution_molecular_pathway(mock_llm_class, mock_agent_class):
    """
    GIVEN a molecular chemistry query
    WHEN graph is invoked
    THEN it should route to molecular agent and return response
    """
    # ARRANGE - Router mock
    mock_router_llm = MagicMock()
    mock_structured_llm = MagicMock()
    mock_router_llm.with_structured_output.return_value = mock_structured_llm

    mock_route_result = MagicMock()
    mock_route_result.destination = "molecular"
    mock_structured_llm.invoke.return_value = mock_route_result

    mock_llm_class.return_value = mock_router_llm

    # ARRANGE - MolecularAgent mock
    mock_agent = MagicMock()
    mock_agent.chat.return_value = iter(["Molecular response"])
    mock_agent_class.return_value = mock_agent

    # Create graph
    graph = create_dispatcher_graph()

    # ACT
    input_state = {
        "messages": [HumanMessage(content="Calculate HOMO for benzene")]
    }
    result = graph.invoke(input_state)

    # ASSERT
    assert "messages" in result
    assert len(result["messages"]) == 2
    assert isinstance(result["messages"][1], AIMessage)
    assert result["messages"][1].content == "Molecular response"


@patch('agent.research.agent.create_research_agent_runnable')
@patch('agent.graph.ChatGoogleGenerativeAI')
def test_graph_execution_research_pathway(mock_llm_class, mock_create_runnable):
    """
    GIVEN a research/paper query
    WHEN graph is invoked
    THEN it should route to research agent and return response
    """
    # ARRANGE - Router mock
    mock_router_llm = MagicMock()
    mock_structured_llm = MagicMock()
    mock_router_llm.with_structured_output.return_value = mock_structured_llm

    mock_route_result = MagicMock()
    mock_route_result.destination = "research"
    mock_structured_llm.invoke.return_value = mock_route_result

    mock_llm_class.return_value = mock_router_llm

    # ARRANGE - ResearchAgent mock
    mock_runnable = MagicMock()
    mock_response = MagicMock()
    mock_response.content = "Found 3 papers about DFT"
    mock_runnable.invoke.return_value = mock_response
    mock_create_runnable.return_value = mock_runnable

    # Create graph
    graph = create_dispatcher_graph()

    # ACT
    input_state = {
        "messages": [HumanMessage(content="Find papers about density functional theory")]
    }
    result = graph.invoke(input_state)

    # ASSERT
    assert "messages" in result
    assert len(result["messages"]) == 2
    assert isinstance(result["messages"][1], AIMessage)
    assert result["messages"][1].content == "Found 3 papers about DFT"


# ============================================================================
# State Management Tests
# ============================================================================

def test_agent_state_preserves_messages():
    """
    GIVEN AgentState with multiple messages
    WHEN new messages are appended
    THEN all messages should be preserved
    """
    # ARRANGE
    state: AgentState = {
        "messages": [
            HumanMessage(content="First message"),
            AIMessage(content="First response")
        ],
        "next_node": ""
    }

    # ACT
    state["messages"].append(HumanMessage(content="Second message"))

    # ASSERT
    assert len(state["messages"]) == 3
    assert state["messages"][0].content == "First message"
    assert state["messages"][1].content == "First response"
    assert state["messages"][2].content == "Second message"


# ============================================================================
# Edge Case Tests
# ============================================================================

@patch('agent.molecular_agent.MolecularAgent')
def test_molecular_agent_node_with_empty_state(mock_agent_class):
    """
    GIVEN state with no messages
    WHEN molecular_agent_node is called
    THEN it should handle gracefully
    """
    # ARRANGE
    state = {
        "messages": [],
        "next_node": ""
    }

    # ACT
    result = molecular_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 1
    assert "No messages found in conversation state" in result["messages"][0].content


@patch('agent.research.agent.create_research_agent_runnable')
def test_research_agent_node_with_empty_state(mock_create_runnable):
    """
    GIVEN state with no messages
    WHEN research_agent_node is called
    THEN it should handle gracefully
    """
    # ARRANGE
    state = {
        "messages": [],
        "next_node": ""
    }

    # ACT
    result = research_agent_node(state)

    # ASSERT
    assert len(result["messages"]) == 1
    assert "No messages found in conversation state" in result["messages"][0].content
