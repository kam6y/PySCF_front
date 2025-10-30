"""
Integration tests for AI Agent API endpoints.

Tests the chat endpoint which uses Server-Sent Events (SSE) for streaming
responses from the LangGraph multi-agent dispatcher.
"""

import pytest
import json
from langchain_core.messages import HumanMessage, AIMessage


class TestAgentChatAPI:
    """Integration tests for POST /api/agent/chat endpoint."""

    def test_chat_success_with_streaming(self, client, mocker):
        """
        GIVEN LangGraph dispatcher is available and returns streaming response
        WHEN POST /api/agent/chat is called with valid message
        THEN SSE stream is returned with chunks and done event
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        # Mock the compiled graph
        mock_graph = mocker.MagicMock()

        # Mock the stream method to return state updates with AI response in supervisor format
        def mock_stream(graph_input, stream_mode="updates", subgraphs=True):
            # Simulate supervisor returning final state with AI message
            user_msg = graph_input["messages"][-1]
            ai_response = AIMessage(content="This is a test response.")
            final_state = {
                "supervisor": {
                    "messages": graph_input["messages"] + [ai_response]
                }
            }
            yield final_state

        mock_graph.stream.side_effect = mock_stream

        # Patch get_compiled_graph to return our mock
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

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
        THEN graph receives the history converted to LangChain format
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        chat_history = [
            {'role': 'user', 'parts': [{'text': 'Hello'}]},
            {'role': 'model', 'parts': [{'text': 'Hi there!'}]}
        ]

        # Mock the compiled graph
        mock_graph = mocker.MagicMock()

        def mock_stream(graph_input, stream_mode="updates", subgraphs=True):
            # Verify history was converted properly (should have 2 history messages + 1 new user message)
            num_messages = len(graph_input["messages"])
            ai_response = AIMessage(content=f"Received {num_messages} total messages (including history)")
            final_state = {
                "supervisor": {
                    "messages": graph_input["messages"] + [ai_response]
                }
            }
            yield final_state

        mock_graph.stream.side_effect = mock_stream
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Continue conversation',
            'history': chat_history
        })

        # ASSERT
        assert response.status_code == 200
        mock_graph.stream.assert_called_once()
        call_args = mock_graph.stream.call_args
        # Should have 3 messages total: 2 from history + 1 new user message
        assert len(call_args[0][0]["messages"]) == 3

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

    def test_chat_graph_error_handling(self, client, mocker):
        """
        GIVEN LangGraph raises an exception
        WHEN POST /api/agent/chat is called
        THEN error event is sent in SSE stream
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        # Mock graph to raise an exception
        mock_graph = mocker.MagicMock()
        mock_graph.stream.side_effect = Exception("Graph execution error")
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

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
        
        # Should contain error event
        assert len(lines) >= 1
        events = [json.loads(line.replace('data: ', '')) for line in lines]
        error_events = [e for e in events if e['type'] == 'error']
        assert len(error_events) > 0

    def test_chat_routing_to_molecular_agent(self, client, mocker):
        """
        GIVEN a molecular chemistry query
        WHEN POST /api/agent/chat is called
        THEN graph routes to molecular agent and returns response
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        mock_graph = mocker.MagicMock()

        def mock_stream(graph_input, stream_mode="updates", subgraphs=True):
            # Simulate routing to molecular agent via supervisor
            ai_response = AIMessage(content="HOMO energy calculation result: -0.5 Hartree")
            final_state = {
                "supervisor": {
                    "messages": graph_input["messages"] + [ai_response]
                }
            }
            yield final_state

        mock_graph.stream.side_effect = mock_stream
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Calculate HOMO energy for water',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Should contain molecular calculation response
        chunk_events = [json.loads(line.replace('data: ', '')) for line in lines if 'chunk' in line]
        assert len(chunk_events) > 0
        
        # Combine all chunks
        full_response = ''.join([e['payload']['text'] for e in chunk_events if e['type'] == 'chunk'])
        assert 'HOMO' in full_response or 'Hartree' in full_response

    def test_chat_routing_to_research_agent(self, client, mocker):
        """
        GIVEN a research/paper search query
        WHEN POST /api/agent/chat is called
        THEN graph routes to research agent and returns response
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        mock_graph = mocker.MagicMock()

        def mock_stream(graph_input, stream_mode="updates", subgraphs=True):
            # Simulate routing to research agent via supervisor
            ai_response = AIMessage(content="Found 3 papers on DFT:\n1. Paper Title 1\n2. Paper Title 2")
            final_state = {
                "supervisor": {
                    "messages": graph_input["messages"] + [ai_response]
                }
            }
            yield final_state

        mock_graph.stream.side_effect = mock_stream
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

        # ACT
        response = client.post('/api/agent/chat', json={
            'message': 'Find papers about density functional theory',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        data_str = response.data.decode('utf-8')
        lines = [line for line in data_str.split('\n') if line.startswith('data:')]
        
        # Should contain research response
        chunk_events = [json.loads(line.replace('data: ', '')) for line in lines if 'chunk' in line]
        assert len(chunk_events) > 0
        
        # Combine all chunks
        full_response = ''.join([e['payload']['text'] for e in chunk_events if e['type'] == 'chunk'])
        assert 'papers' in full_response.lower() or 'DFT' in full_response

    def test_chat_history_optional(self, client, mocker):
        """
        GIVEN history is empty list
        WHEN POST /api/agent/chat is called
        THEN graph receives only the current message
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        mock_graph = mocker.MagicMock()

        def mock_stream(graph_input, stream_mode="updates", subgraphs=True):
            ai_response = AIMessage(content="Response without history")
            final_state = {
                "supervisor": {
                    "messages": graph_input["messages"] + [ai_response]
                }
            }
            yield final_state

        mock_graph.stream.side_effect = mock_stream
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

        # ACT - empty history list
        response = client.post('/api/agent/chat', json={
            'message': 'Test message',
            'history': []
        })

        # ASSERT
        assert response.status_code == 200
        mock_graph.stream.assert_called_once()
        call_args = mock_graph.stream.call_args
        # Should have only 1 message (the current user message, no history)
        assert len(call_args[0][0]["messages"]) == 1

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
        GIVEN graph returns a long response
        WHEN POST /api/agent/chat is called
        THEN chunks are received in correct order
        """
        # ARRANGE
        # Clear the compiled graph cache
        import api.agent
        api.agent._compiled_graph_cache = None

        mock_graph = mocker.MagicMock()

        def mock_stream(graph_input, stream_mode="updates", subgraphs=True):
            # Return a long response that will be chunked
            long_response = "First Second Third"
            ai_response = AIMessage(content=long_response)
            final_state = {
                "supervisor": {
                    "messages": graph_input["messages"] + [ai_response]
                }
            }
            yield final_state

        mock_graph.stream.side_effect = mock_stream
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

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
        chunk_events = [json.loads(line.replace('data: ', '')) for line in lines if 'chunk' in line]
        
        # Verify we received chunks
        assert len(chunk_events) > 0
        # Combine all chunks to verify complete message
        full_response = ''.join([e['payload']['text'] for e in chunk_events if e['type'] == 'chunk'])
        assert full_response == "First Second Third"

    def test_chat_response_format(self, client, mocker):
        """
        GIVEN graph returns response
        WHEN POST /api/agent/chat is called
        THEN SSE events have correct format
        """
        # ARRANGE
        mock_graph = mocker.MagicMock()
        
        def mock_stream(graph_input):
            ai_response = AIMessage(content="Test response")
            final_state = {
                "molecular": {
                    "messages": graph_input["messages"] + [ai_response],
                    "next_node": "end"
                }
            }
            yield final_state
        
        mock_graph.stream.side_effect = mock_stream
        mocker.patch('api.agent.get_compiled_graph', return_value=mock_graph)

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
            # agent_status is now also a valid event type
            assert event['type'] in ['chunk', 'done', 'error', 'agent_status']

            if event['type'] == 'chunk':
                assert 'payload' in event
                assert 'text' in event['payload']
            elif event['type'] == 'agent_status':
                assert 'payload' in event
                assert 'status' in event['payload']
                assert 'agent' in event['payload']
