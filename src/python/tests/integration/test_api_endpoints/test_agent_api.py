"""
Integration tests for AI Agent API endpoints.

Tests the chat endpoint which uses Server-Sent Events (SSE) for streaming
responses from the Gemini API.
"""

import json
import logging

from generated_models import AgentChatRequest
from api.agent import (
    stream_chat_response,
)


def test_agent_chat_stream_preserves_sse_format(client, mocker):
    def fake_stream(*args, **kwargs):
        yield {"type": "agent_status", "payload": {"status": "started"}}
        yield {"type": "chunk", "payload": {"text": "hello"}}
        yield {"type": "done", "payload": {"message_id": "msg-1"}}

    mocker.patch("api.agent.stream_chat_response", side_effect=fake_stream)

    with client.stream(
        "POST",
        "/api/agent/chat",
        json={"message": "Hello", "history": [], "session_id": "session-1"},
    ) as response:
        chunks = "".join(response.iter_text())

    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/event-stream")
    assert (
        'data: {"type": "agent_status", "payload": {"status": "started"}}\n\n' in chunks
    )
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
        response = client.post("/api/agent/chat", json={"message": "", "history": []})

        # ASSERT
        assert response.status_code == 400

    def test_chat_whitespace_only_message(self, client):
        """
        GIVEN message with only whitespace
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "   ", "history": []}
        )

        # ASSERT
        assert response.status_code == 400

    def test_chat_message_too_long(self, client):
        """
        GIVEN message exceeds MAX_MESSAGE_LENGTH (10 000 chars)
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE — MAX_MESSAGE_LENGTH is 10_000
        very_long_message = "a" * 10_001

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": very_long_message, "history": []}
        )

        # ASSERT
        assert response.status_code == 400
        assert "10000" in response.json()["error"]

    def test_chat_missing_message_field(self, client):
        """
        GIVEN request is missing message field
        WHEN POST /api/agent/chat is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post("/api/agent/chat", json={"history": []})

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
            "/api/agent/chat",
            content=b"invalid json",
            headers={"Content-Type": "application/json"},
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
        mock_settings_service = mocker.patch("api.agent.SettingsService")
        mock_settings_service.return_value.get_settings.return_value = {
            "gemini_api_key": "test-api-key"
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
        mocker.patch.dict("sys.modules", {"google.generativeai": mock_genai})

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "What is water?", "history": []}
        )

        # ASSERT
        assert response.status_code == 200
        assert response.headers["content-type"].startswith("text/event-stream")

        # Parse SSE stream
        data_str = response.content.decode("utf-8")
        lines = [line for line in data_str.split("\n") if line.startswith("data:")]

        # Should have chunk events and a done event
        assert len(lines) >= 2  # At least some chunks + done

        # Verify last event is 'done'
        last_event = json.loads(lines[-1].replace("data: ", ""))
        assert last_event["type"] == "done"

    def test_chat_no_api_key(self, client, mocker):
        """
        GIVEN Gemini API key is not configured
        WHEN POST /api/agent/chat is called
        THEN error event is returned in SSE stream
        """
        # ARRANGE
        mock_settings_service = mocker.patch("api.agent.SettingsService")
        mock_settings_service.return_value.get_settings.return_value = {
            "gemini_api_key": None
        }

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "What is water?", "history": []}
        )

        # ASSERT
        assert response.status_code == 200
        assert response.headers["content-type"].startswith("text/event-stream")

        # Parse SSE stream
        data_str = response.content.decode("utf-8")
        lines = [line for line in data_str.split("\n") if line.startswith("data:")]

        # Should contain error event about missing API key
        events = [json.loads(line.replace("data: ", "")) for line in lines]
        error_events = [e for e in events if e["type"] == "error"]
        assert len(error_events) > 0
        assert "API key" in error_events[0]["payload"]["message"]

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


class TestAgentChatExtraFieldRejection:
    """Verify that AgentChatRequest with extra='forbid' rejects unknown fields."""

    def test_agent_chat_rejects_unknown_fields(self, client):
        """
        GIVEN a valid chat payload with an extra unknown field
        WHEN POST /api/agent/chat is called
        THEN 400 is returned because AgentChatRequest has extra='forbid'
        """
        payload = {
            "message": "Hello",
            "history": [],
            "unknown_field": "should be rejected",
        }

        response = client.post("/api/agent/chat", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]


class TestAgentChatInputSizeLimits:
    """Tests for security-related input size validation on POST /api/agent/chat.

    Constants under test (from api.agent):
        MAX_MESSAGE_LENGTH          = 10_000
        MAX_HISTORY_ITEMS           = 200
        MAX_TOTAL_HISTORY_CHARS     = 500_000
        MAX_HISTORY_PART_TEXT_LENGTH = 10_000
        HistoryItem.parts max_length = 100  (Pydantic model cap)
    """

    # -- helpers ----------------------------------------------------------

    @staticmethod
    def _setup_gemini_mock(mocker):
        """Configure Gemini SDK mock so valid requests stream without network."""
        mock_settings_service = mocker.patch("api.agent.SettingsService")
        mock_settings_service.return_value.get_settings.return_value = {
            "gemini_api_key": "test-api-key"
        }

        mock_genai = mocker.MagicMock()
        mock_model = mocker.MagicMock()
        mock_chat = mocker.MagicMock()
        mock_chunk = mocker.MagicMock()
        mock_chunk.text = "OK"

        mock_chat.send_message.return_value = iter([mock_chunk])
        mock_model.start_chat.return_value = mock_chat
        mock_genai.GenerativeModel.return_value = mock_model
        mock_genai.configure = mocker.MagicMock()

        mocker.patch.dict("sys.modules", {"google.generativeai": mock_genai})

    @staticmethod
    def _make_history_item(text="hi", role="user"):
        """Return a single history item dict."""
        return {"role": role, "parts": [{"text": text}]}

    # -- 1. Message length boundary --------------------------------------

    def test_chat_message_at_max_length_accepted(self, client, mocker):
        """
        GIVEN message is exactly MAX_MESSAGE_LENGTH (10 000) characters
        WHEN POST /api/agent/chat is called
        THEN request is accepted (200 streaming response)
        """
        # ARRANGE
        self._setup_gemini_mock(mocker)
        message = "a" * 10_000

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": message, "history": []}
        )

        # ASSERT
        assert response.status_code == 200

    def test_chat_message_one_over_max_length_rejected(self, client):
        """
        GIVEN message is MAX_MESSAGE_LENGTH + 1 (10 001) characters
        WHEN POST /api/agent/chat is called
        THEN 400 is returned with an error mentioning 10000
        """
        # ARRANGE
        message = "a" * 10_001

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": message, "history": []}
        )

        # ASSERT
        assert response.status_code == 400
        assert "10000" in response.json()["error"]

    # -- 2. History item count boundary -----------------------------------

    def test_chat_history_at_max_items_accepted(self, client, mocker):
        """
        GIVEN history has exactly MAX_HISTORY_ITEMS (200) entries
        WHEN POST /api/agent/chat is called
        THEN request is accepted (200 streaming response)
        """
        # ARRANGE
        self._setup_gemini_mock(mocker)
        history = [self._make_history_item() for _ in range(200)]

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "test", "history": history}
        )

        # ASSERT
        assert response.status_code == 200

    def test_chat_history_one_over_max_items_rejected(self, client):
        """
        GIVEN history has MAX_HISTORY_ITEMS + 1 (201) entries
        WHEN POST /api/agent/chat is called
        THEN 400 is returned with an error mentioning 200
        """
        # ARRANGE
        history = [self._make_history_item() for _ in range(201)]

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "test", "history": history}
        )

        # ASSERT
        assert response.status_code == 400
        assert "200" in response.json()["error"]

    # -- 3. Per-part text length ------------------------------------------

    def test_chat_history_part_text_over_max_rejected(self, client):
        """
        GIVEN a single history part has MAX_HISTORY_PART_TEXT_LENGTH + 1 (10 001) chars
        WHEN POST /api/agent/chat is called
        THEN 400 is returned with an error mentioning 10000
        """
        # ARRANGE
        history = [self._make_history_item(text="a" * 10_001)]

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "test", "history": history}
        )

        # ASSERT
        assert response.status_code == 400
        assert "10000" in response.json()["error"]

    # -- 4. Total history characters --------------------------------------

    def test_chat_total_history_chars_over_max_rejected(self, client):
        """
        GIVEN total chars across all history parts exceeds MAX_TOTAL_HISTORY_CHARS
              (51 items * 10 000 chars each = 510 000 > 500 000)
              Each individual part is within per-part cap, item count within cap.
        WHEN POST /api/agent/chat is called
        THEN 400 is returned with an error about total history text
        """
        # ARRANGE — 51 items, each with one 10_000-char part → 510_000 total
        history = [self._make_history_item(text="a" * 10_000) for _ in range(51)]

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "test", "history": history}
        )

        # ASSERT
        assert response.status_code == 400
        assert "500000" in response.json()["error"]

    # -- 5. Parts maxItems (Pydantic model cap: 100) ----------------------

    def test_chat_history_item_parts_over_model_max_rejected(self, client):
        """
        GIVEN a single history item has 101 parts (Pydantic max_length=100)
        WHEN POST /api/agent/chat is called
        THEN 400 is returned (Pydantic validation mapped to 400 by app handler)
        """
        # ARRANGE
        parts = [{"text": "x"} for _ in range(101)]
        history = [{"role": "user", "parts": parts}]

        # ACT
        response = client.post(
            "/api/agent/chat", json={"message": "test", "history": history}
        )

        # ASSERT
        assert response.status_code == 400
        assert "100" in response.json()["error"]


class TestAgentChatDebugLogging:
    """SEC-006: Verify chat message content is not leaked into debug logs."""

    def test_debug_log_does_not_contain_message_content(self, mocker, caplog):
        """
        GIVEN a chat request with a distinctive message
        WHEN stream_chat_response is consumed
        THEN the debug log contains message length but NOT the message text
        """
        # ARRANGE
        secret_message = "TopSecretContent_XYZ_12345"
        mock_settings_service = mocker.patch("api.agent.SettingsService")
        mock_settings_service.return_value.get_settings.return_value = {
            "gemini_api_key": "test-api-key"
        }

        mock_genai = mocker.MagicMock()
        mock_model = mocker.MagicMock()
        mock_chat = mocker.MagicMock()
        mock_chunk = mocker.MagicMock()
        mock_chunk.text = "OK"
        mock_chat.send_message.return_value = iter([mock_chunk])
        mock_model.start_chat.return_value = mock_chat
        mock_genai.GenerativeModel.return_value = mock_model
        mock_genai.configure = mocker.MagicMock()
        mocker.patch.dict("sys.modules", {"google.generativeai": mock_genai})

        request = AgentChatRequest(
            message=secret_message,
            history=[],
            session_id=None,
        )

        # ACT
        with caplog.at_level(logging.DEBUG, logger="api.agent"):
            events = list(stream_chat_response(request))

        # ASSERT — stream completed
        assert any(e["type"] == "done" for e in events)
        # The message text must NOT appear in any log record
        for record in caplog.records:
            assert (
                secret_message not in record.getMessage()
            ), f"Message content leaked into log: {record.getMessage()}"
        # The message length SHOULD appear in debug logs
        debug_messages = [
            r.getMessage() for r in caplog.records if r.levelno == logging.DEBUG
        ]
        assert len(debug_messages) > 0, "No DEBUG records captured; test is vacuous"
        assert any(str(len(secret_message)) in msg for msg in debug_messages)


class TestExtractTextPartsNullSafety:
    """Verify _extract_text_parts skips None and non-str values (F1 soundness fix)."""

    def test_dict_part_with_none_text_is_skipped(self):
        """Dict part whose text value is None must not appear in the result."""

        from api.agent import _extract_text_parts

        parts = [{"text": None}, {"text": "hello"}]
        result = _extract_text_parts(parts)
        assert result == ["hello"]

    def test_dict_part_with_non_str_text_is_skipped(self):
        """Dict part whose text value is a non-str (e.g. int) must be skipped."""
        from api.agent import _extract_text_parts

        parts = [{"text": 123}, {"text": "world"}]
        result = _extract_text_parts(parts)
        assert result == ["world"]

    def test_pydantic_part_with_none_text_is_skipped(self):
        """Pydantic-like object with text=None must be skipped."""
        import types

        from api.agent import _extract_text_parts

        obj_none = types.SimpleNamespace(text=None)
        obj_valid = types.SimpleNamespace(text="valid")
        result = _extract_text_parts([obj_none, obj_valid])
        assert result == ["valid"]

    def test_convert_history_skips_none_text_parts(self):
        """End-to-end: history item with a None text part does not crash."""
        import types

        from api.agent import _convert_history_to_gemini_format

        history = [
            {
                "role": "user",
                "parts": [{"text": None}, {"text": "hi"}],
            },
            {
                "role": "model",
                "parts": [types.SimpleNamespace(text=None)],
            },
        ]
        result = _convert_history_to_gemini_format(history)
        # First item keeps the valid "hi" part; second item is empty and dropped
        assert len(result) == 1
        assert result[0] == {"role": "user", "parts": ["hi"]}


class TestAgentStreamErrorRedaction:
    """IMP-T2: Verify SSE error payload is redacted and partial responses
    are not persisted when the LLM stream raises mid-response."""

    def test_stream_error_redacts_internal_detail_and_skips_db_save(self, mocker):
        """
        GIVEN Gemini send_message yields one chunk then raises RuntimeError
        WHEN stream_chat_response is consumed to completion
        THEN (a) the SSE error event uses the generic message, not the
             internal detail, and (b) no model message is persisted.
        """
        session_id = "session-error-redact"
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

        # Yield one chunk then raise
        def _failing_stream(*args, **kwargs):
            mock_chunk = mocker.MagicMock()
            mock_chunk.text = "Partial"
            yield mock_chunk
            raise RuntimeError("secret internal detail")

        mock_chat.send_message.side_effect = _failing_stream
        mock_model.start_chat.return_value = mock_chat
        mock_genai.GenerativeModel.return_value = mock_model
        mocker.patch.dict("sys.modules", {"google.generativeai": mock_genai})

        request = AgentChatRequest(
            message="What is water?",
            history=[],
            session_id=session_id,
        )

        # ACT — consume the entire generator
        events = list(stream_chat_response(request))

        # ASSERT — error event has generic message
        error_events = [e for e in events if e.get("type") == "error"]
        assert len(error_events) == 1
        error_msg = error_events[0]["payload"]["message"]
        assert error_msg == "AI chat failed. Check settings and logs."
        assert "secret internal detail" not in error_msg

        # No "done" event (stream did not complete normally)
        assert not any(e.get("type") == "done" for e in events)

        # No model message persisted (stream_completed is False)
        model_saves = [
            call
            for call in mock_chat_service.add_message.call_args_list
            if call.args[1] == "model"
        ]
        assert model_saves == []
