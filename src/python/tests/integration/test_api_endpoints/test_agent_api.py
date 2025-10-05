"""
Integration tests for AI Agent API endpoints.

Tests the chat endpoint which uses Server-Sent Events (SSE) for streaming
responses from the molecular agent.
"""

import pytest
import json


class TestAgentChatAPI:
    """Integration tests for POST /api/agent/chat endpoint."""

    def test_chat_success_with_streaming(self, client, mocker):
        """
        GIVEN MolecularAgent is available and returns streaming response
        WHEN POST /api/agent/chat is called with valid message
        THEN SSE stream is returned with chunks and done event
        """
        # ARRANGE
        # Mock the agent's chat method to return an iterator
        def mock_chat_iterator(message, history):
            yield "This is "
            yield "a test "
            yield "response."
        
        mock_agent = mocker.MagicMock()
        mock_agent.chat.return_value = mock_chat_iterator("test", [])
        
        # Patch get_molecular_agent to return our mock
        mocker.patch('api.agent.get_molecular_agent', return_value=mock_agent)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'What is water?',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        assert response.content_type == 'text/event-stream'
        
        # Parse SSE stream
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Should have chunk events and a done event
        assert len(lines) >= 2  # At least some chunks + done
        
        # Verify last event is 'done'
        last_event = json.loads(lines[-1].replace('data: ', ''))
        assert last_event['type'] == 'done'

    def test_chat_with_history(self, client, mocker):
        """
        GIVEN chat history is provided
        WHEN POST /api/agent/chat is called
        THEN agent receives the history
        """
        # ARRANGE
        chat_history = [
            {'role': 'user', 'parts': [{'text': 'Hello'}]},
            {'role': 'model', 'parts': [{'text': 'Hi there!'}]}
        ]
        
        def mock_chat_iterator(message, history):
            yield f"History length: {len(history)}"
        
        mock_agent = mocker.MagicMock()
        mock_agent.chat.return_value = mock_chat_iterator("test", chat_history)
        mocker.patch('api.agent.get_molecular_agent', return_value=mock_agent)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Continue conversation',
            'history': chat_history
        })

        # ASSERT
        assert response.status_code == 200
        mock_agent.chat.assert_called_once()
        call_args = mock_agent.chat.call_args
        assert call_args[0][0] == 'Continue conversation'
        assert len(call_args[0][1]) == 2

    def test_chat_empty_message(self, client):
        """
        GIVEN empty message
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/agent/chat', json={
            'message': '',
            'history': []
        })

        # ASSERT
        assert response.status_code == 400

    def test_chat_whitespace_only_message(self, client):
        """
        GIVEN message with only whitespace
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/agent/chat', json={
            'message': '   ',
            'history': []
        })

        # ASSERT
        assert response.status_code == 400

    def test_chat_message_too_long(self, client):
        """
        GIVEN message exceeds maximum length
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        max_length = 100000
        very_long_message = 'a' * (max_length + 1)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': very_long_message,
            'history': []
        })

        # ASSERT
        assert response.status_code == 400

    def test_chat_missing_message_field(self, client):
        """
        GIVEN request is missing message field
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/agent/chat', json={
            'history': []
        })

        # ASSERT
        assert response.status_code == 400

    def test_chat_agent_not_available_fallback(self, client, mocker):
        """
        GIVEN MolecularAgent is not available
        WHEN POST /api/agent/chat is called
        THEN fallback SSE stream is returned
        """
        # ARRANGE
        # Patch get_molecular_agent to return None (agent not available)
        mocker.patch('api.agent.get_molecular_agent', return_value=None)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'What is water?',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        assert response.content_type == 'text/event-stream'
        
        # Parse SSE stream
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Should contain fallback message
        assert len(lines) >= 1
        first_event = json.loads(lines[0].replace('data: ', ''))
        assert 'not available' in first_event['payload']['text'].lower()

    def test_chat_agent_error_during_streaming(self, client, mocker):
        """
        GIVEN agent raises exception during streaming
        WHEN POST /api/agent/chat is called
        THEN error event is sent in SSE stream
        """
        # ARRANGE
        def mock_chat_error(message, history):
            yield "Starting..."
            raise Exception("Agent internal error")
        
        mock_agent = mocker.MagicMock()
        mock_agent.chat.return_value = mock_chat_error("test", [])
        mocker.patch('api.agent.get_molecular_agent', return_value=mock_agent)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Cause an error',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Should contain error event
        events = [json.loads(line.replace('data: ', '')) for line in lines]
        error_events = [e for e in events if e['type'] == 'error']
        assert len(error_events) > 0

    def test_chat_history_optional(self, client, mocker):
        """
        GIVEN history is empty list
        WHEN POST /api/agent/chat is called
        THEN agent receives empty history
        """
        # ARRANGE
        def mock_chat_iterator(message, history):
            yield "Response without history"
        
        mock_agent = mocker.MagicMock()
        mock_agent.chat.return_value = mock_chat_iterator("test", [])
        mocker.patch('api.agent.get_molecular_agent', return_value=mock_agent)

        # ACT - empty history list
        response = client.post('/api/agent/chat', json={
            'message': 'Test message',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        mock_agent.chat.assert_called_once()
        call_args = mock_agent.chat.call_args
        # History should be empty list
        assert call_args[0][1] == []

    def test_chat_invalid_json(self, client):
        """
        GIVEN invalid JSON payload
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post(
            '/api/agent/chat',
            data='invalid json',
            content_type='application/json'
        )

        # ASSERT
        assert response.status_code == 400

    def test_chat_multiple_chunks_ordering(self, client, mocker):
        """
        GIVEN agent returns multiple chunks
        WHEN POST /api/agent/chat is called
        THEN chunks are received in correct order
        """
        # ARRANGE
        expected_chunks = ["First ", "Second ", "Third"]
        
        def mock_chat_iterator(message, history):
            for chunk in expected_chunks:
                yield chunk
        
        mock_agent = mocker.MagicMock()
        mock_agent.chat.return_value = mock_chat_iterator("test", [])
        mocker.patch('api.agent.get_molecular_agent', return_value=mock_agent)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Test ordering',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Parse chunk events (exclude done event)
        chunk_events = [json.loads(line.replace('data: ', '')) for line in lines[:-1]]
        
        # Verify chunks are in correct order
        assert len(chunk_events) == 3
        for i, expected_chunk in enumerate(expected_chunks):
            assert chunk_events[i]['type'] == 'chunk'
            assert chunk_events[i]['payload']['text'] == expected_chunk

    def test_chat_response_format(self, client, mocker):
        """
        GIVEN agent returns response
        WHEN POST /api/agent/chat is called
        THEN SSE events have correct format
        """
        # ARRANGE
        def mock_chat_iterator(message, history):
            yield "Test response"
        
        mock_agent = mocker.MagicMock()
        mock_agent.chat.return_value = mock_chat_iterator("test", [])
        mocker.patch('api.agent.get_molecular_agent', return_value=mock_agent)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Test format',
            'history': []
        })

        # ASSERT
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Verify each line is valid JSON with expected structure
        for line in lines:
            event = json.loads(line.replace('data: ', ''))
            assert 'type' in event
            assert event['type'] in ['chunk', 'done', 'error']
            
            if event['type'] == 'chunk':
                assert 'payload' in event
                assert 'text' in event['payload']
