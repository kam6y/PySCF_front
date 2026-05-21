"""
Integration tests for AI Agent API endpoints.

Tests the chat endpoint which uses Server-Sent Events (SSE) for streaming
responses from the Gemini API.
"""

import json

from generated_models import AgentChatRequest
from api.agent import stream_chat_response


def test_agent_chat_stream_preserves_sse_format(client, mocker):
    def fake_stream(*args, **kwargs):
        yield {'type': 'agent_status', 'payload': {'status': 'started'}}
        yield {'type': 'chunk', 'payload': {'text': 'hello'}}
        yield {'type': 'done', 'payload': {'message_id': 'msg-1'}}

    mocker.patch('api.agent.stream_chat_response', side_effect=fake_stream)

    with client.stream(
        'POST',
        '/api/agent/chat',
        json={'message': 'Hello', 'history': [], 'session_id': 'session-1'},
    ) as response:
        chunks = ''.join(response.iter_text())

    assert response.status_code == 200
    assert response.headers['content-type'].startswith('text/event-stream')
    assert 'data: {"type": "agent_status", "payload": {"status": "started"}}\n\n' in chunks
    assert 'data: {"type": "chunk", "payload": {"text": "hello"}}\n\n' in chunks
    assert 'data: {"type": "done", "payload": {"message_id": "msg-1"}}\n\n' in chunks


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
            content=b'invalid json',
            headers={'Content-Type': 'application/json'},
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
        # This ensures the import statement inside stream_chat_response gets the mock
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
        assert response.headers['content-type'].startswith('text/event-stream')
        
        # Parse SSE stream
        data_str = response.content.decode('utf-8')
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
        assert response.headers['content-type'].startswith('text/event-stream')
        
        # Parse SSE stream
        data_str = response.content.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Should contain error event about missing API key
        events = [json.loads(line.replace('data: ', '')) for line in lines]
        error_events = [e for e in events if e['type'] == 'error']
        assert len(error_events) > 0
        assert 'API key' in error_events[0]['payload']['message']

    def test_chat_stream_client_abort_does_not_save_partial_model_message(self, mocker):
        """
        GIVEN Gemini API starts streaming a model response for a saved session
        WHEN the SSE generator is closed after the first chunk
        THEN the partial model response is not saved to chat history
        """
        # ARRANGE
        session_id = "session-cancelled"
        mock_chat_service = mocker.MagicMock()
        mocker.patch(
            "api.agent.get_chat_history_service",
            return_value=mock_chat_service,
        )

        mock_settings_service = mocker.patch("api.agent.SettingsService")
        mock_settings_service.return_value.get_settings.return_value = {
            "gemini_api_key": "test-api-key"
        }

        mock_genai = mocker.MagicMock()
        mock_model = mocker.MagicMock()
        mock_chat = mocker.MagicMock()
        mock_chunk = mocker.MagicMock()
        mock_chunk.text = "Partial response"

        mock_chat.send_message.return_value = iter([mock_chunk])
        mock_model.start_chat.return_value = mock_chat
        mock_genai.GenerativeModel.return_value = mock_model
        mocker.patch.dict("sys.modules", {"google.generativeai": mock_genai})

        request = AgentChatRequest(
            message="What is water?",
            history=[],
            session_id=session_id,
        )
        generator = stream_chat_response(request)

        # ACT
        status_event = next(generator)
        chunk_event = next(generator)
        generator.close()

        # ASSERT
        assert status_event["type"] == "agent_status"
        assert chunk_event["type"] == "chunk"
        mock_chat_service.add_message.assert_any_call(
            session_id,
            "user",
            "What is water?",
        )
        model_saves = [
            call
            for call in mock_chat_service.add_message.call_args_list
            if call.args[1] == "model"
        ]
        assert model_saves == []
