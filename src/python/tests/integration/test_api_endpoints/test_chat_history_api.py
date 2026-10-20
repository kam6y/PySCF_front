from services.exceptions import ServiceError


class TestChatHistoryAPI:
    def test_get_chat_sessions_success(self, client, mocker):
        mock_data = {
            "sessions": [{"id": "session-1", "name": "First chat"}],
            "total_count": 1,
        }
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.list_sessions.return_value = mock_data

        response = client.get("/api/chat-history/sessions")

        assert response.status_code == 200
        assert response.json() == {"success": True, "data": mock_data}
        mock_service.return_value.list_sessions.assert_called_once_with()

    def test_create_chat_session_success(self, client, mocker):
        mock_session = {"id": "session-1", "name": "New chat"}
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.create_session.return_value = mock_session

        response = client.post(
            "/api/chat-history/sessions",
            json={"name": "New chat"},
        )

        assert response.status_code == 201
        assert response.json() == {
            "success": True,
            "data": {"session": mock_session},
        }
        mock_service.return_value.create_session.assert_called_once_with(
            name="New chat",
        )

    def test_get_chat_session_success(self, client, mocker):
        mock_data = {
            "session": {"id": "session-1", "name": "First chat"},
            "messages": [{"role": "user", "content": "Hello"}],
        }
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.get_session_with_messages.return_value = mock_data

        response = client.get("/api/chat-history/sessions/session-1")

        assert response.status_code == 200
        assert response.json() == {"success": True, "data": mock_data}
        mock_service.return_value.get_session_with_messages.assert_called_once_with(
            "session-1",
        )

    def test_get_chat_session_not_found(self, client, mocker):
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.get_session_with_messages.return_value = None

        response = client.get("/api/chat-history/sessions/missing-session")

        assert response.status_code == 404
        assert response.json() == {
            "success": False,
            "error": "Chat session not found: missing-session",
        }

    def test_update_chat_session_success(self, client, mocker):
        mock_session = {"id": "session-1", "name": "Renamed chat"}
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.update_session.return_value = mock_session

        response = client.patch(
            "/api/chat-history/sessions/session-1",
            json={"name": "Renamed chat"},
        )

        assert response.status_code == 200
        assert response.json() == {
            "success": True,
            "data": {"session": mock_session},
        }
        mock_service.return_value.update_session.assert_called_once_with(
            "session-1",
            "Renamed chat",
        )

    def test_update_chat_session_not_found(self, client, mocker):
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.update_session.return_value = None

        response = client.patch(
            "/api/chat-history/sessions/missing-session",
            json={"name": "Renamed chat"},
        )

        assert response.status_code == 404
        assert response.json() == {
            "success": False,
            "error": "Chat session not found: missing-session",
        }

    def test_delete_chat_session_success(self, client, mocker):
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.delete_session.return_value = True

        response = client.delete("/api/chat-history/sessions/session-1")

        assert response.status_code == 200
        assert response.json() == {
            "success": True,
            "data": {
                "message": "Chat session deleted successfully",
                "deleted_id": "session-1",
            },
        }
        mock_service.return_value.delete_session.assert_called_once_with("session-1")

    def test_delete_chat_session_not_found(self, client, mocker):
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.delete_session.return_value = False

        response = client.delete("/api/chat-history/sessions/missing-session")

        assert response.status_code == 404
        assert response.json() == {
            "success": False,
            "error": "Chat session not found: missing-session",
        }

    def test_create_chat_session_rejects_unknown_fields(self, client):
        """
        GIVEN a valid create-session payload with an extra unknown field
        WHEN POST /api/chat-history/sessions is called
        THEN 400 is returned because CreateChatSessionRequest has extra='forbid'
        """
        payload = {"name": "New chat", "unknown_field": "should be rejected"}

        response = client.post("/api/chat-history/sessions", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]

    def test_update_chat_session_rejects_unknown_fields(self, client):
        """
        GIVEN a valid update-session payload with an extra unknown field
        WHEN PATCH /api/chat-history/sessions/<id> is called
        THEN 400 is returned because UpdateChatSessionRequest has extra='forbid'
        """
        payload = {"name": "Renamed chat", "unknown_field": "should be rejected"}

        response = client.patch(
            "/api/chat-history/sessions/session-1",
            json=payload,
        )

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]

    def test_service_error_uses_global_error_envelope(self, client, mocker):
        mock_service = mocker.patch("api.chat_history.get_chat_history_service")
        mock_service.return_value.list_sessions.side_effect = ServiceError(
            "Chat history unavailable",
            status_code=503,
        )

        response = client.get("/api/chat-history/sessions")

        assert response.status_code == 503
        # 503 is a curated operational status code -- the developer-authored
        # message is safe to return to the client (unlike generic 500 which
        # is redacted to prevent leaking raw exception text).
        assert response.json() == {
            "success": False,
            "error": "Chat history unavailable",
        }
