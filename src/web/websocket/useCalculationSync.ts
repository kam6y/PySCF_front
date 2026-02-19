import { useCallback, useRef } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import { Socket } from 'socket.io-client';
import { CalculationInstance } from '../types/api-types';
import { handleError } from '../utils/errorHandler';
import { calculationQueryKeys } from '../hooks/useCalculationQueries';
import { invalidateQueriesWithRetry } from './invalidateQueriesWithRetry';
import type { SocketEventHandlers } from './useSocketTransport';

export interface UseCalculationSyncOptions {
  activeCalculationId: string | null;
  onCalculationUpdate?: (
    updated: CalculationInstance,
    previousStatus: string | undefined
  ) => void;
  onWebSocketError?: (error: string) => void;
}

export interface UseCalculationSyncReturn {
  transportHandlers: SocketEventHandlers;
  manageActiveCalculationRoom: (
    newActiveId: string | null,
    emit: (event: string, data?: unknown) => void
  ) => void;
}

export const useCalculationSync = ({
  activeCalculationId,
  onCalculationUpdate,
  onWebSocketError,
}: UseCalculationSyncOptions): UseCalculationSyncReturn => {
  const queryClient = useQueryClient();
  const activeCalculationIdRef = useRef(activeCalculationId);
  activeCalculationIdRef.current = activeCalculationId;
  const currentActiveCalculationId = useRef<string | null>(null);
  const lastNotifiedStatus = useRef<Map<string, string>>(new Map());
  const lastNotifiedUpdatedAt = useRef<Map<string, string>>(new Map());
  const syncFailureCountRef = useRef<number>(0);
  const isConnectedRef = useRef<boolean>(false);

  // コールバックを ref 化し、依存チェーンを安定化
  const onCalculationUpdateRef = useRef(onCalculationUpdate);
  onCalculationUpdateRef.current = onCalculationUpdate;
  const onWebSocketErrorRef = useRef(onWebSocketError);
  onWebSocketErrorRef.current = onWebSocketError;

  // 統一されたキャッシュ更新ロジック
  const updateCalculationCache = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;

      try {
        // キャッシュを無効化して再取得を促す（データの整合性はサーバーが保証）
        // 1. 個別計算詳細のキャッシュを無効化
        queryClient.invalidateQueries({
          queryKey: calculationQueryKeys.detail(calculationId),
        });

        // 2. 計算リストのキャッシュを無効化
        queryClient.invalidateQueries({
          queryKey: calculationQueryKeys.list(),
        });

        console.log(
          `[UnifiedWebSocket] Invalidated queries for calculation ${calculationId}`
        );
      } catch (error) {
        console.error(
          `[UnifiedWebSocket] Failed to invalidate queries for calculation ${calculationId}:`,
          error
        );
      }
    },
    [queryClient]
  );

  const notifyWebSocketError = useCallback(
    (errorMessage: string) => {
      queryClient.invalidateQueries({
        queryKey: calculationQueryKeys.list(),
      });
      onWebSocketErrorRef.current?.(errorMessage);
    },
    [queryClient]
  );

  // 統一された計算更新処理
  const handleCalculationUpdate = useCallback(
    (updatedCalculation: CalculationInstance) => {
      const calculationId = updatedCalculation.id;

      if (!updatedCalculation || !calculationId) {
        console.warn(
          '[UnifiedWebSocket] Received invalid calculation update data:',
          updatedCalculation
        );
        return;
      }

      console.log(
        `[UnifiedWebSocket] Processing update for calculation ${calculationId}: ${updatedCalculation.status}`
      );

      // updatedAt ベースで重複を判定（同一イベントの二重受信を排除）
      // status ベースの判定は廃止: running 状態中の中間結果（scf_iterations等）を取りこぼすため
      const previousStatus = lastNotifiedStatus.current.get(calculationId);
      const previousUpdatedAt =
        lastNotifiedUpdatedAt.current.get(calculationId);

      if (previousUpdatedAt === updatedCalculation.updatedAt) {
        console.log(
          `[UnifiedWebSocket] Skipping duplicate update for ${calculationId}: updatedAt unchanged (${updatedCalculation.updatedAt})`
        );
        return;
      }

      lastNotifiedStatus.current.set(calculationId, updatedCalculation.status);
      lastNotifiedUpdatedAt.current.set(
        calculationId,
        updatedCalculation.updatedAt
      );

      console.log(
        `[UnifiedWebSocket] Update detected for ${calculationId}: status=${previousStatus ?? 'none'} -> ${updatedCalculation.status}, updatedAt=${updatedCalculation.updatedAt}`
      );

      // キャッシュ更新
      updateCalculationCache(updatedCalculation);

      // 通知処理
      onCalculationUpdateRef.current?.(updatedCalculation, previousStatus);
    },
    [queryClient, updateCalculationCache]
  );

  const handleSocketError = useCallback(
    (errorData: unknown) => {
      console.error('[UnifiedWebSocket] Socket error:', errorData);

      let errorMessage = 'A problem occurred with calculation monitoring.';
      if (typeof errorData === 'object' && errorData !== null) {
        const maybeError = (errorData as Record<string, unknown>).error;
        if (typeof maybeError === 'string' && maybeError.length > 0) {
          errorMessage = maybeError;
        }
      }

      // If calculation not found (likely due to directory change), silently ignore
      if (errorMessage.includes('not found')) {
        console.log(
          '[UnifiedWebSocket] Calculation not found (likely directory changed), ignoring error'
        );
        return;
      }

      notifyWebSocketError(errorMessage);
    },
    [notifyWebSocketError]
  );

  const handleConnectError = useCallback(
    (error: Error) => {
      // その他のエラーのみ通知
      let errorMessage = 'Failed to connect to monitoring server.';

      if (error.message.includes('timeout')) {
        errorMessage = 'Connection timed out. The server may not be running.';
      } else if (error.message.includes('xhr poll error')) {
        errorMessage = 'Server communication was interrupted.';
      } else if (error.message.includes('400')) {
        errorMessage =
          'Server rejected connection (HTTP 400). Check CORS settings.';
      } else if (error.message.includes('403')) {
        errorMessage =
          'Connection forbidden (HTTP 403). Check server authentication.';
      } else if (error.message.includes('network')) {
        errorMessage = 'Unable to connect to server due to network error.';
      } else if (error.message.includes('security')) {
        errorMessage = 'Connection blocked by security settings.';
      }

      notifyWebSocketError(errorMessage);
    },
    [notifyWebSocketError]
  );

  const handleReconnect = useCallback(
    async (_socket: Socket, attemptNumber: number) => {
      console.log(`[UnifiedWebSocket] Reconnected after ${attemptNumber} attempts`);
      console.log('[UnifiedWebSocket] Syncing data after reconnection...');

      const activeId = activeCalculationIdRef.current;

      try {
        // 並行実行で効率化しつつ、各クエリにリトライロジックを適用
        await Promise.all([
          invalidateQueriesWithRetry({
            queryClient,
            queryKey: [...calculationQueryKeys.list()],
          }),
          activeId && !activeId.startsWith('new-calculation-')
            ? invalidateQueriesWithRetry({
                queryClient,
                queryKey: [...calculationQueryKeys.detail(activeId)],
              })
            : Promise.resolve(),
        ]);

        syncFailureCountRef.current = 0;
        console.log('[UnifiedWebSocket] Data sync completed successfully');
      } catch (error) {
        syncFailureCountRef.current++;
        console.error(
          `[UnifiedWebSocket] Data sync failed (attempt ${syncFailureCountRef.current}/2):`,
          error
        );

        if (syncFailureCountRef.current >= 2) {
          handleError(
            new Error(
              'Unable to sync calculation data after reconnection. Please refresh the page if calculations appear outdated.'
            ),
            'Connection synchronization failed'
          );
          syncFailureCountRef.current = 0;
        }
        // エラーでもアプリケーションは継続（次のWebSocket更新で回復可能）
      }
    },
    [queryClient]
  );

  // transportHandlers は参照を固定し、内部の処理は ref 経由で差し替える
  const handleCalculationUpdateRef = useRef(handleCalculationUpdate);
  handleCalculationUpdateRef.current = handleCalculationUpdate;
  const handleSocketErrorRef = useRef(handleSocketError);
  handleSocketErrorRef.current = handleSocketError;
  const handleConnectErrorRef = useRef(handleConnectError);
  handleConnectErrorRef.current = handleConnectError;
  const handleReconnectRef = useRef(handleReconnect);
  handleReconnectRef.current = handleReconnect;

  const transportHandlersRef = useRef<SocketEventHandlers>({
    onConnect: socket => {
      isConnectedRef.current = true;
      console.log('[UnifiedWebSocket] Connected, joining global_updates room');
      socket.emit('join_global_updates');

      const activeId = activeCalculationIdRef.current;
      if (activeId && !activeId.startsWith('new-calculation-')) {
        console.log(`[UnifiedWebSocket] Joining calculation room: ${activeId}`);
        socket.emit('join_calculation', { calculation_id: activeId });
        currentActiveCalculationId.current = activeId;
      }
    },
    onDisconnect: _reason => {
      isConnectedRef.current = false;
      currentActiveCalculationId.current = null;
    },
    onReconnect: (socket, attemptNumber) => {
      return handleReconnectRef.current(socket, attemptNumber);
    },
    onConnectError: error => {
      return handleConnectErrorRef.current(error);
    },
    onError: errorData => {
      return handleSocketErrorRef.current(errorData);
    },
    dataListeners: {
      calculation_update: data => {
        return handleCalculationUpdateRef.current(data as CalculationInstance);
      },
    },
  });

  const manageActiveCalculationRoom = useCallback(
    (newActiveId: string | null, emit: (event: string, data?: unknown) => void) => {
      if (!isConnectedRef.current) return;

      const previousId = currentActiveCalculationId.current;

      if (previousId === newActiveId) {
        activeCalculationIdRef.current = newActiveId;
        currentActiveCalculationId.current = newActiveId;
        return;
      }

      // 前のルームから退出
      if (previousId && !previousId.startsWith('new-calculation-')) {
        console.log(`[UnifiedWebSocket] Leaving calculation room: ${previousId}`);
        emit('leave_calculation', { calculation_id: previousId });
      }

      // 新しいルームに参加
      if (newActiveId && !newActiveId.startsWith('new-calculation-')) {
        console.log(`[UnifiedWebSocket] Joining calculation room: ${newActiveId}`);
        emit('join_calculation', { calculation_id: newActiveId });
      }

      activeCalculationIdRef.current = newActiveId;
      currentActiveCalculationId.current = newActiveId;
    },
    []
  );

  return {
    transportHandlers: transportHandlersRef.current,
    manageActiveCalculationRoom,
  };
};
