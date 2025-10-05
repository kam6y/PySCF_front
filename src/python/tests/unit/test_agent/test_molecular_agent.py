"""
Unit tests for Molecular Agent.

Tests the MolecularAgent class with mocked Google Gemini API to verify
chat handling, Function Calling workflows, and error handling.
"""

import pytest
from unittest.mock import MagicMock, patch, PropertyMock
from google.genai.errors import APIError, ServerError

from agent.molecular_agent import MolecularAgent
from generated_models import HistoryItem, Role, Part


# ============================================================================
# Initialization Tests
# ============================================================================

@patch('agent.molecular_agent.genai.Client')
def test_molecular_agent_initialization_with_api_key(mock_client_class):
    """
    GIVEN GEMINI_API_KEY environment variable is set
    WHEN MolecularAgent is initialized
    THEN it should create genai.Client with the API key
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    # ACT
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-api-key'}):
        agent = MolecularAgent()
    
    # ASSERT
    assert agent.client is not None
    assert agent.api_key == 'test-api-key'
    mock_client_class.assert_called_once_with(api_key='test-api-key')


@patch('agent.molecular_agent.get_current_settings')
@patch('agent.molecular_agent.genai.Client')
def test_molecular_agent_initialization_without_api_key(mock_client_class, mock_settings):
    """
    GIVEN GEMINI_API_KEY environment variable is not set
    WHEN MolecularAgent is initialized
    THEN it should not create a client and use fallback mode
    """
    # ARRANGE
    mock_settings.return_value.gemini_api_key = None

    # ACT
    with patch.dict('os.environ', {}, clear=True):
        agent = MolecularAgent()

    # ASSERT
    assert agent.client is None
    assert agent.api_key is None
    mock_client_class.assert_not_called()


@patch('agent.molecular_agent.genai.Client')
def test_molecular_agent_tools_initialization(mock_client_class):
    """
    GIVEN MolecularAgent is initialized with API key
    WHEN agent is created
    THEN it should initialize all available tools
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    # ACT
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ASSERT
    assert len(agent.available_tools) > 0
    # Verify some expected tools are present
    tool_names = [tool.__name__ for tool in agent.available_tools]
    assert 'list_all_calculations' in tool_names
    assert 'get_calculation_details' in tool_names
    assert 'start_quantum_calculation' in tool_names


# ============================================================================
# Chat - Simple Text Response Tests
# ============================================================================

@patch('agent.molecular_agent.genai.Client')
def test_chat_simple_text_response(mock_client_class):
    """
    GIVEN Gemini API returns a simple text response
    WHEN chat is called
    THEN it should stream the text response
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    # Create mock response chunks
    chunk1 = MagicMock()
    chunk1.text = "Hello! "
    chunk2 = MagicMock()
    chunk2.text = "I can help you with molecular calculations."
    
    mock_response_stream = [chunk1, chunk2]
    mock_client.models.generate_content_stream.return_value = mock_response_stream
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ACT
    response_chunks = list(agent.chat("Hello", []))
    
    # ASSERT
    assert len(response_chunks) == 2
    assert response_chunks[0] == "Hello! "
    assert response_chunks[1] == "I can help you with molecular calculations."


@patch('agent.molecular_agent.genai.Client')
def test_chat_with_history(mock_client_class):
    """
    GIVEN chat history is provided
    WHEN chat is called
    THEN it should convert history to Gemini format and include it
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    chunk = MagicMock()
    chunk.text = "Response"
    mock_client.models.generate_content_stream.return_value = [chunk]
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # Create history in frontend format
    history = [
        HistoryItem(
            role=Role.user,
            parts=[Part(text="Previous question")]
        ),
        HistoryItem(
            role=Role.model,
            parts=[Part(text="Previous answer")]
        )
    ]
    
    # ACT
    list(agent.chat("New question", history))
    
    # ASSERT
    # Verify generate_content_stream was called
    mock_client.models.generate_content_stream.assert_called_once()
    call_args = mock_client.models.generate_content_stream.call_args
    
    # Check that contents parameter includes history
    contents = call_args.kwargs['contents']
    assert len(contents) == 3  # 2 history + 1 new message
    assert contents[0]['role'] == 'user'
    assert contents[0]['parts'][0]['text'] == "Previous question"
    assert contents[1]['role'] == 'model'
    assert contents[2]['role'] == 'user'
    assert contents[2]['parts'][0]['text'] == "New question"


# ============================================================================
# Chat - Fallback Response Tests
# ============================================================================

@patch('agent.molecular_agent.get_current_settings')
@patch('agent.molecular_agent.genai.Client')
def test_chat_fallback_when_no_client(mock_client_class, mock_settings):
    """
    GIVEN MolecularAgent without API key (no client)
    WHEN chat is called
    THEN it should return fallback response
    """
    # ARRANGE
    mock_settings.return_value.gemini_api_key = None

    with patch.dict('os.environ', {}, clear=True):
        agent = MolecularAgent()

    # ACT
    response_chunks = list(agent.chat("Hello", []))

    # ASSERT
    assert len(response_chunks) == 1
    response = response_chunks[0]
    assert "not fully configured" in response.lower()
    assert "gemini_api_key" in response.lower()


# ============================================================================
# Chat - Error Handling Tests
# ============================================================================

@patch('agent.molecular_agent.genai.Client')
def test_chat_server_error_503_retry_exhausted(mock_client_class):
    """
    GIVEN Gemini API returns 503 errors on all retry attempts
    WHEN chat is called
    THEN it should return appropriate error message after retries
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client

    error = ServerError(code=503, response_json={'error': {'message': 'Service unavailable'}})
    mock_client.models.generate_content_stream.side_effect = error

    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()

    # ACT
    response_chunks = list(agent.chat("Hello", []))

    # ASSERT
    assert len(response_chunks) == 1
    response = response_chunks[0]
    assert "high demand" in response.lower() or "overloaded" in response.lower()

    # Verify retries were attempted (default is 3)
    assert mock_client.models.generate_content_stream.call_count == 3


@patch('agent.molecular_agent.genai.Client')
def test_chat_server_error_429_rate_limit(mock_client_class):
    """
    GIVEN Gemini API returns 429 rate limit error
    WHEN chat is called
    THEN it should retry and return appropriate error message
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client

    error = ServerError(code=429, response_json={'error': {'message': 'Too many requests'}})
    mock_client.models.generate_content_stream.side_effect = error

    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()

    # ACT
    response_chunks = list(agent.chat("Hello", []))

    # ASSERT
    response = response_chunks[0]
    assert "too many requests" in response.lower() or "wait" in response.lower()


@patch('agent.molecular_agent.genai.Client')
def test_chat_api_error_authentication(mock_client_class):
    """
    GIVEN Gemini API returns authentication error
    WHEN chat is called
    THEN it should return authentication error message
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client

    error = APIError(code=401, response_json={'error': {'message': 'Invalid API key'}})
    mock_client.models.generate_content_stream.side_effect = error

    with patch.dict('os.environ', {'GEMINI_API_KEY': 'invalid-key'}):
        agent = MolecularAgent()

    # ACT
    response_chunks = list(agent.chat("Hello", []))

    # ASSERT
    response = response_chunks[0]
    assert "authentication" in response.lower() or "api key" in response.lower()


@patch('agent.molecular_agent.genai.Client')
def test_chat_unexpected_error(mock_client_class):
    """
    GIVEN Gemini API raises unexpected exception
    WHEN chat is called
    THEN it should return generic error message
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    mock_client.models.generate_content_stream.side_effect = RuntimeError('Unexpected error')
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ACT
    response_chunks = list(agent.chat("Hello", []))
    
    # ASSERT
    response = response_chunks[0]
    assert "temporary issue" in response.lower() or "cannot process" in response.lower()


# ============================================================================
# API Key Reload Tests
# ============================================================================

@patch('agent.molecular_agent.get_current_settings')
@patch('agent.molecular_agent.genai.Client')
def test_reload_api_key_success(mock_client_class, mock_settings):
    """
    GIVEN MolecularAgent without API key initially
    WHEN reload_api_key is called with new API key in environment
    THEN it should reinitialize the client
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    mock_settings.return_value.gemini_api_key = None

    # Start without API key
    with patch.dict('os.environ', {}, clear=True):
        agent = MolecularAgent()

    assert agent.client is None

    # ACT - Add API key and reload
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'new-api-key'}):
        success = agent.reload_api_key()

    # ASSERT
    assert success is True
    assert agent.api_key == 'new-api-key'
    assert agent.client is not None


@patch('agent.molecular_agent.genai.Client')
def test_reload_api_key_no_change(mock_client_class):
    """
    GIVEN MolecularAgent with API key
    WHEN reload_api_key is called with same API key
    THEN it should not reinitialize (returns current status)
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    initial_client = agent.client
    initial_call_count = mock_client_class.call_count
    
    # ACT - Reload with same key
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        success = agent.reload_api_key()
    
    # ASSERT
    assert success is True
    assert agent.client is initial_client  # Same client instance
    assert mock_client_class.call_count == initial_call_count  # Not called again


@patch('agent.molecular_agent.get_current_settings')
@patch('agent.molecular_agent.genai.Client')
def test_reload_api_key_failure(mock_client_class, mock_settings):
    """
    GIVEN API key is removed from environment
    WHEN reload_api_key is called
    THEN it should return False and clear client
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    mock_settings.return_value.gemini_api_key = None

    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()

    assert agent.client is not None

    # ACT - Remove API key and reload
    with patch.dict('os.environ', {}, clear=True):
        success = agent.reload_api_key()

    # ASSERT
    assert success is False
    assert agent.client is None
    assert agent.api_key is None


# ============================================================================
# Configuration Tests
# ============================================================================

@patch('agent.molecular_agent.genai.Client')
def test_agent_loads_config(mock_client_class):
    """
    GIVEN server-config.json exists with ai_agent configuration
    WHEN MolecularAgent is initialized
    THEN it should load and use the configuration
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    # ACT
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ASSERT
    assert agent.config is not None
    assert 'model_name' in agent.config
    assert 'max_retries' in agent.config
    # Default values should be set
    assert agent.config['model_name'] in ['gemini-2.5-flash', 'gemini-2.0-flash-exp']
    assert agent.config['max_retries'] >= 1


@patch('agent.molecular_agent.genai.Client')
def test_agent_uses_system_prompt(mock_client_class):
    """
    GIVEN MolecularAgent is initialized
    WHEN chat is called
    THEN it should use system prompt in generate_content call
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    chunk = MagicMock()
    chunk.text = "Response"
    mock_client.models.generate_content_stream.return_value = [chunk]
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ACT
    list(agent.chat("Hello", []))
    
    # ASSERT
    call_args = mock_client.models.generate_content_stream.call_args
    assert 'config' in call_args.kwargs
    config = call_args.kwargs['config']
    assert hasattr(config, 'system_instruction') or 'system_instruction' in dir(config)


# ============================================================================
# Edge Cases and Validation Tests
# ============================================================================

@patch('agent.molecular_agent.genai.Client')
def test_chat_empty_message(mock_client_class):
    """
    GIVEN an empty message
    WHEN chat is called
    THEN it should still attempt to generate response
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    chunk = MagicMock()
    chunk.text = "Response"
    mock_client.models.generate_content_stream.return_value = [chunk]
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ACT & ASSERT - Should not raise exception
    response_chunks = list(agent.chat("", []))
    assert len(response_chunks) >= 1


@patch('agent.molecular_agent.genai.Client')
def test_chat_empty_history(mock_client_class):
    """
    GIVEN empty history list
    WHEN chat is called
    THEN it should handle it correctly
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    chunk = MagicMock()
    chunk.text = "Response"
    mock_client.models.generate_content_stream.return_value = [chunk]
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ACT
    response_chunks = list(agent.chat("Hello", []))
    
    # ASSERT
    assert len(response_chunks) >= 1
    # Verify only the current message was sent
    call_args = mock_client.models.generate_content_stream.call_args
    contents = call_args.kwargs['contents']
    assert len(contents) == 1  # Only current message


@patch('agent.molecular_agent.genai.Client')
def test_chat_handles_none_text_chunks(mock_client_class):
    """
    GIVEN response stream contains chunks with no text
    WHEN chat is called
    THEN it should skip chunks without text
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client_class.return_value = mock_client
    
    chunk1 = MagicMock()
    chunk1.text = "Hello"
    chunk2 = MagicMock()
    chunk2.text = None  # No text in this chunk
    chunk3 = MagicMock()
    chunk3.text = "World"
    
    mock_client.models.generate_content_stream.return_value = [chunk1, chunk2, chunk3]
    
    with patch.dict('os.environ', {'GEMINI_API_KEY': 'test-key'}):
        agent = MolecularAgent()
    
    # ACT
    response_chunks = list(agent.chat("Hello", []))
    
    # ASSERT
    # Should only yield chunks with text
    assert len(response_chunks) == 2
    assert response_chunks[0] == "Hello"
    assert response_chunks[1] == "World"
