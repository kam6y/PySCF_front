import { components } from '../types/generated-api';
import { request, ApiError } from './core';

type ChatHistoryListResponse = components['schemas']['ChatHistoryListResponse'];
type ChatSessionResponse = components['schemas']['ChatSessionResponse'];
type ChatSessionDetailResponse =
  components['schemas']['ChatSessionDetailResponse'];

const rejectInvalidSessionId = (sessionId: string): Promise<never> => {
  return Promise.reject(
    new ApiError(
      'Invalid session ID provided.',
      400,
      'Bad Request',
      `/api/chat-history/sessions/${sessionId}`,
      null,
      false
    )
  );
};

export const getChatSessions = (): Promise<ChatHistoryListResponse['data']> => {
  return request<ChatHistoryListResponse['data']>(
    '/api/chat-history/sessions',
    {
      method: 'GET',
    }
  );
};

export const createChatSession = (
  name?: string
): Promise<ChatSessionResponse['data']> => {
  return request<ChatSessionResponse['data']>('/api/chat-history/sessions', {
    method: 'POST',
    body: JSON.stringify({ name }),
  });
};

export const getChatSessionDetail = (
  sessionId: string
): Promise<ChatSessionDetailResponse['data']> => {
  if (!sessionId || sessionId.trim() === '') {
    return rejectInvalidSessionId(sessionId);
  }

  return request<ChatSessionDetailResponse['data']>(
    `/api/chat-history/sessions/${sessionId}`,
    { method: 'GET' }
  );
};

export const updateChatSession = (
  sessionId: string,
  name: string
): Promise<ChatSessionResponse['data']> => {
  if (!sessionId || sessionId.trim() === '') {
    return rejectInvalidSessionId(sessionId);
  }

  return request<ChatSessionResponse['data']>(
    `/api/chat-history/sessions/${sessionId}`,
    {
      method: 'PATCH',
      body: JSON.stringify({ name }),
    }
  );
};

export const deleteChatSession = (sessionId: string): Promise<void> => {
  if (!sessionId || sessionId.trim() === '') {
    return rejectInvalidSessionId(sessionId);
  }

  return request<void>(`/api/chat-history/sessions/${sessionId}`, {
    method: 'DELETE',
  });
};
