"""
Integration tests for AI Agent API endpoints.

Tests the chat endpoint which uses Server-Sent Events (SSE) for streaming
responses from the Gemini API.
"""

import json


class TestAgentChatAPI:
    """Integration tests for POST /api/agent/chat endpoint."""

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

    def test_chat_success_with_streaming_mock(self, client, mocker):
        """
        GIVEN Gemini API is mocked to return streaming response
        WHEN POST /api/agent/chat is called with valid message
        THEN SSE stream is returned with chunks and done event
        """
        # ARRANGE
        # Mock settings service to return API key
        mock_settings_service = mocker.patch('api.agent.SettingsService')
        mock_settings_service.return_value.get_settings.return_value = {
            'gemini_api_key': 'test-api-key'
        }

        # Mock Gemini model by injecting into sys.modules
        # This ensures the import statement inside _create_simple_chat_stream gets the mock
        mock_genai = mocker.MagicMock()
        mock_model = mocker.MagicMock()
        mock_chat = mocker.MagicMock()
        mock_response = mocker.MagicMock()

        # Configure mock response chunks
        mock_chunk1 = mocker.MagicMock()
        mock_chunk1.text = "Hello"
        mock_chunk2 = mocker.MagicMock()
        mock_chunk2.text = " World"
        mock_response.__iter__ = lambda self: iter([mock_chunk1, mock_chunk2])

        mock_chat.send_message.return_value = mock_response
        mock_model.start_chat.return_value = mock_chat
        mock_genai.GenerativeModel.return_value = mock_model
        mock_genai.configure = mocker.MagicMock()  # Mock the configure function

        # Inject mock into sys.modules so import statement gets the mock
        mocker.patch.dict('sys.modules', {'google.generativeai': mock_genai})

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

    def test_chat_no_api_key(self, client, mocker):
        """
        GIVEN Gemini API key is not configured
        WHEN POST /api/agent/chat is called
        THEN error event is returned in SSE stream
        """
        # ARRANGE
        mock_settings_service = mocker.patch('api.agent.SettingsService')
        mock_settings_service.return_value.get_settings.return_value = {
            'gemini_api_key': None
        }

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
        
        # Should contain error event about missing API key
        events = [json.loads(line.replace('data: ', '')) for line in lines]
        error_events = [e for e in events if e['type'] == 'error']
        assert len(error_events) > 0
        assert 'API key' in error_events[0]['payload']['message']
