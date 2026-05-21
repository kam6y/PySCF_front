"""Chat history API endpoints."""

from fastapi import APIRouter
from fastapi.responses import JSONResponse

from generated_models import CreateChatSessionRequest, UpdateChatSessionRequest
from services.chat_history_service import get_chat_history_service

router = APIRouter(prefix='/api/chat-history')


@router.get('/sessions')
def get_chat_sessions() -> dict:
    data = get_chat_history_service().list_sessions()
    return {'success': True, 'data': data}


@router.post('/sessions', status_code=201)
def create_chat_session(body: CreateChatSessionRequest) -> dict:
    session = get_chat_history_service().create_session(name=body.name)
    return {'success': True, 'data': {'session': session}}


@router.get('/sessions/{session_id}', response_model=None)
def get_chat_session(session_id: str) -> dict | JSONResponse:
    session_data = get_chat_history_service().get_session_with_messages(session_id)
    if session_data is None:
        return JSONResponse(
            {'success': False, 'error': f'Chat session not found: {session_id}'},
            status_code=404,
        )
    return {'success': True, 'data': session_data}


@router.patch('/sessions/{session_id}', response_model=None)
def update_chat_session(
    session_id: str,
    body: UpdateChatSessionRequest,
) -> dict | JSONResponse:
    session = get_chat_history_service().update_session(session_id, body.name)
    if session is None:
        return JSONResponse(
            {'success': False, 'error': f'Chat session not found: {session_id}'},
            status_code=404,
        )
    return {'success': True, 'data': {'session': session}}


@router.delete('/sessions/{session_id}', response_model=None)
def delete_chat_session(session_id: str) -> dict | JSONResponse:
    deleted = get_chat_history_service().delete_session(session_id)
    if not deleted:
        return JSONResponse(
            {'success': False, 'error': f'Chat session not found: {session_id}'},
            status_code=404,
        )
    return {
        'success': True,
        'data': {
            'message': 'Chat session deleted successfully',
            'deleted_id': session_id,
        },
    }
