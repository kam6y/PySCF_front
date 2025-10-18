import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import {
  getChatSessions,
  createChatSession,
  getChatSessionDetail,
  updateChatSession,
  deleteChatSession,
} from '../apiClient';

// Query keys
export const chatHistoryKeys = {
  all: ['chatHistory'] as const,
  sessions: () => [...chatHistoryKeys.all, 'sessions'] as const,
  sessionDetail: (sessionId: string) =>
    [...chatHistoryKeys.all, 'session', sessionId] as const,
};

/**
 * Get all chat sessions
 */
export const useGetChatSessions = () => {
  return useQuery({
    queryKey: chatHistoryKeys.sessions(),
    queryFn: getChatSessions,
    staleTime: 30000, // 30 seconds
  });
};

/**
 * Get chat session detail
 */
export const useGetChatSessionDetail = (sessionId: string | null) => {
  return useQuery({
    queryKey: chatHistoryKeys.sessionDetail(sessionId || ''),
    queryFn: () => {
      if (!sessionId) {
        return Promise.reject(new Error('No session ID provided'));
      }
      return getChatSessionDetail(sessionId);
    },
    enabled: !!sessionId,
    staleTime: 30000, // 30 seconds - balance between freshness and performance
    refetchOnWindowFocus: true, // Re-enabled for better UX - refetch when user returns to tab
    // Conditional retry: Only retry on 404 errors once (newly created sessions may not be immediately available)
    retry: (failureCount, error) => {
      // Don't retry on 404 (not found) - session doesn't exist
      if (error instanceof Error && error.message.includes('404')) {
        return false;
      }
      // Retry once for network errors
      return failureCount < 1;
    },
  });
};

/**
 * Create a new chat session
 */
export const useCreateChatSession = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: (name?: string) => createChatSession(name),
    onSuccess: () => {
      // Invalidate sessions list to refetch
      queryClient.invalidateQueries({ queryKey: chatHistoryKeys.sessions() });
    },
  });
};

/**
 * Update chat session name
 */
export const useUpdateChatSession = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: ({ sessionId, name }: { sessionId: string; name: string }) =>
      updateChatSession(sessionId, name),
    onSuccess: (_, variables) => {
      // Invalidate sessions list and specific session detail
      queryClient.invalidateQueries({ queryKey: chatHistoryKeys.sessions() });
      queryClient.invalidateQueries({
        queryKey: chatHistoryKeys.sessionDetail(variables.sessionId),
      });
    },
  });
};

/**
 * Delete a chat session
 */
export const useDeleteChatSession = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: (sessionId: string) => deleteChatSession(sessionId),
    onSuccess: (_, sessionId) => {
      // Invalidate sessions list and remove specific session from cache
      queryClient.invalidateQueries({ queryKey: chatHistoryKeys.sessions() });
      queryClient.removeQueries({
        queryKey: chatHistoryKeys.sessionDetail(sessionId),
      });
    },
  });
};
