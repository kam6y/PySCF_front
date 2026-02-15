// src/web/pages/DrawMoleculePage.tsx

import React, { useState, useRef, useEffect, useCallback, useMemo, memo } from 'react';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Ketcher } from 'ketcher-core';
import 'ketcher-react/dist/index.css';
import styles from './DrawMoleculePage.module.css';
import { convertSmilesToXyz } from '../api/molecule';
import { useUIStore } from '../store/uiStore';
import { useCalculationStore } from '../store/calculationStore';
import { useNotificationStore } from '../store/notificationStore';
import { useActiveCalculation } from '../hooks/useActiveCalculation';
import {
  DEFAULT_CALCULATION_PARAMETERS,
  STATUS_MESSAGES,
} from '../constants/calculationDefaults';
import { isCalculationEditable } from '../utils/calculationStatus';

type KetcherWindow = Window & {
  ketcher?: {
    logging?: {
      enabled?: boolean;
      level?: number;
      showTrace?: boolean;
    };
    editor?: {
      errorHandler?: (...args: unknown[]) => void;
    };
    [key: string]: unknown;
  };
};

const noopErrorHandler = () => {};
const RESTORE_MAX_RETRIES = 10;
const RESTORE_RETRY_DELAY_MS = 300;

type PersistenceFormat = 'molfile' | 'ket' | 'smiles';

const hasMeaningfulContent = (value: string | null | undefined): boolean =>
  typeof value === 'string' && value.trim().length > 0;

const getRestoreKey = (
  calculationId: string | undefined,
  ketcherData: string | null | undefined
): string => `${calculationId ?? 'no-calculation'}::${ketcherData ?? ''}`;

const serializeForPersistence = async (
  ketcher: Ketcher,
  smilesFallback?: string
): Promise<{ data: string; format: PersistenceFormat }> => {
  const serializers: Array<{
    format: PersistenceFormat;
    getter: () => Promise<string>;
  }> = [
    { format: 'molfile', getter: () => ketcher.getMolfile() },
    { format: 'ket', getter: () => ketcher.getKet() },
    {
      format: 'smiles',
      getter: async () => smilesFallback ?? ketcher.getSmiles(),
    },
  ];

  for (const serializer of serializers) {
    try {
      const serializedData = await serializer.getter();
      if (hasMeaningfulContent(serializedData)) {
        return { data: serializedData, format: serializer.format };
      }
    } catch {
      // Try next serializer
    }
  }

  throw new Error('Failed to serialize molecule data for persistence.');
};

const ensureKetcherGlobals = () => {
  if (typeof window === 'undefined') {
    return;
  }

  const browserWindow = window as KetcherWindow;
  if (!browserWindow.ketcher || typeof browserWindow.ketcher !== 'object') {
    browserWindow.ketcher = {};
  }

  const ketcherGlobal = browserWindow.ketcher;
  const logging =
    ketcherGlobal.logging && typeof ketcherGlobal.logging === 'object'
      ? ketcherGlobal.logging
      : {};

  ketcherGlobal.logging = {
    enabled: typeof logging.enabled === 'boolean' ? logging.enabled : false,
    level: typeof logging.level === 'number' ? logging.level : 0,
    showTrace: typeof logging.showTrace === 'boolean' ? logging.showTrace : false,
  };

  if (!ketcherGlobal.editor || typeof ketcherGlobal.editor !== 'object') {
    try {
      ketcherGlobal.editor = { errorHandler: noopErrorHandler };
    } catch {
      // Ignore assignment errors for readonly getter-based editor objects.
    }
  }

  if (ketcherGlobal.editor && typeof ketcherGlobal.editor === 'object') {
    try {
      if (typeof ketcherGlobal.editor.errorHandler !== 'function') {
        ketcherGlobal.editor.errorHandler = noopErrorHandler;
      }
    } catch {
      // Ignore assignment errors for readonly error handler fields.
    }
  }
};

// KetcherLoggerがwindow.ketcherを参照するため、初期表示前から最小限のグローバルを保証する
ensureKetcherGlobals();

interface KetcherEditorErrorBoundaryProps {
  children: React.ReactNode;
  onRetry: () => void;
}

interface KetcherEditorErrorBoundaryState {
  hasError: boolean;
}

class KetcherEditorErrorBoundary extends React.Component<
  KetcherEditorErrorBoundaryProps,
  KetcherEditorErrorBoundaryState
> {
  constructor(props: KetcherEditorErrorBoundaryProps) {
    super(props);
    this.state = { hasError: false };
  }

  static getDerivedStateFromError(): KetcherEditorErrorBoundaryState {
    return { hasError: true };
  }

  componentDidCatch(error: Error) {
    console.error('[Ketcher Error Boundary] Editor crashed', error);
  }

  private handleRetry = () => {
    this.setState({ hasError: false });
    this.props.onRetry();
  };

  render() {
    if (this.state.hasError) {
      return (
        <div className={styles.editorFallback}>
          <h3 className={styles.editorFallbackTitle}>Editor crashed</h3>
          <p className={styles.editorFallbackText}>
            The molecule editor failed to initialize. Please retry.
          </p>
          <button
            type="button"
            className={styles.editorRetryButton}
            onClick={this.handleRetry}
          >
            Retry Editor
          </button>
        </div>
      );
    }

    return this.props.children;
  }
}

// Miewをwindowに設定（Ketcherが3D表示に使用）
if (typeof window !== 'undefined') {
  import('miew').then(Miew => {
    (window as any).Miew = Miew.default || Miew;
  });
}

// Ketcherエディタコンポーネント（マウントごとにStructServiceProviderを新規作成）
const KetcherEditor = memo<{
  errorHandler: (message: string) => void;
  onInit: (ketcher: Ketcher) => void;
}>(({ errorHandler, onInit }) => {
  // マウントごとに新しいProviderを作成（アンマウント時にKetcherが内部状態をクリーンアップするため再利用不可）
  const structServiceProvider = useMemo(
    () => new StandaloneStructServiceProvider(),
    []
  );

  return (
    <Editor
      staticResourcesUrl=""
      structServiceProvider={structServiceProvider}
      errorHandler={errorHandler}
      onInit={onInit}
    />
  );
});

export const DrawMoleculePage: React.FC = () => {
  const ketcherInstanceRef = useRef<Ketcher | null>(null);
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);
  const [isKetcherReady, setIsKetcherReady] = useState(false);
  const [ketcherInitVersion, setKetcherInitVersion] = useState(0);
  const [editorRetryKey, setEditorRetryKey] = useState(0);
  const hasRestoredRef = useRef<string | null>(null); // 復元済みキー（calculation ID + ketcher_data）
  const restoreTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const clearRestoreTimer = useCallback(() => {
    if (restoreTimerRef.current) {
      clearTimeout(restoreTimerRef.current);
      restoreTimerRef.current = null;
    }
  }, []);

  // コンポーネントアンマウント時のタイマークリーンアップ
  useEffect(() => {
    return () => {
      clearRestoreTimer();
    };
  }, [clearRestoreTimer]);

  // Zustandストア
  const setCurrentPage = useUIStore(state => state.setCurrentPage);
  const setStagedCalculation = useCalculationStore(
    state => state.setStagedCalculation
  );
  const setActiveCalculationId = useCalculationStore(
    state => state.setActiveCalculationId
  );

  // アクティブな計算を取得
  const { activeCalculation } = useActiveCalculation();

  // Draw Moleculeを使わずに開始された計算かどうか判定
  const isNonDrawMoleculeCalculation = Boolean(
    activeCalculation &&
      activeCalculation.status !== 'pending' &&
      !activeCalculation.parameters?.ketcher_data
  );

  // 編集可否の判定（pending/completed/error/undefined のみ編集可）
  const canEdit = isCalculationEditable(activeCalculation?.status);

  // Ketcherインスタンスの初期化（メモ化して安定化）
  const handleOnInit = useCallback((ketcher: Ketcher) => {
    if (typeof window !== 'undefined') {
      const browserWindow = window as KetcherWindow;
      browserWindow.ketcher = ketcher as unknown as KetcherWindow['ketcher'];
      ensureKetcherGlobals();
    }
    ketcherInstanceRef.current = ketcher;
    hasRestoredRef.current = null;
    setIsKetcherReady(true);
    // Ketcher内部で再初期化が発生しても復元処理を再トリガーする
    setKetcherInitVersion(prev => prev + 1);
  }, []);

  // activeCalculation変更時にKetcherデータを復元
  useEffect(() => {
    if (isNonDrawMoleculeCalculation) {
      return;
    }

    let cancelled = false;

    const restoreKetcherData = async (attempt: number = 0) => {
      if (!ketcherInstanceRef.current || !isKetcherReady) {
        return;
      }

      const currentCalcId = activeCalculation?.id;
      const ketcherData = activeCalculation?.parameters?.ketcher_data;
      const restoreKey = getRestoreKey(currentCalcId, ketcherData);
      const payload = ketcherData ?? '';

      // 既に同じ復元キーを適用済みの場合はスキップ
      if (hasRestoredRef.current === restoreKey) {
        return;
      }

      try {
        await ketcherInstanceRef.current.setMolecule(payload);
        if (cancelled) {
          return;
        }
        hasRestoredRef.current = restoreKey;
      } catch (error) {
        if (attempt < RESTORE_MAX_RETRIES - 1) {
          const nextAttempt = attempt + 1;
          clearRestoreTimer();
          restoreTimerRef.current = setTimeout(() => {
            if (cancelled) {
              return;
            }
            restoreKetcherData(nextAttempt);
          }, RESTORE_RETRY_DELAY_MS);
          return;
        }

        console.error('[Ketcher Restore] Failed to apply molecule state after retries', {
          calculationId: currentCalcId,
          attempts: RESTORE_MAX_RETRIES,
          error,
        });
        useNotificationStore.getState().addNotification({
          type: 'error',
          title: 'Restore Error',
          message: 'Failed to restore molecule structure after retries',
          autoClose: true,
          duration: 3000,
        });
      }
    };

    restoreKetcherData();
    return () => {
      cancelled = true;
      clearRestoreTimer();
    };
  }, [
    activeCalculation?.id,
    isKetcherReady,
    ketcherInitVersion,
    activeCalculation?.parameters?.ketcher_data,
    isNonDrawMoleculeCalculation,
    clearRestoreTimer,
  ]);

  useEffect(() => {
    if (!isNonDrawMoleculeCalculation) {
      return;
    }

    clearRestoreTimer();
    ketcherInstanceRef.current = null;
    setIsKetcherReady(false);
    hasRestoredRef.current = null;
  }, [isNonDrawMoleculeCalculation, clearRestoreTimer]);

  const handleRetryEditor = useCallback(() => {
    clearRestoreTimer();
    ensureKetcherGlobals();
    setConvertError(null);
    setIsKetcherReady(false);
    ketcherInstanceRef.current = null;
    hasRestoredRef.current = null;
    setEditorRetryKey(prev => prev + 1);
  }, [clearRestoreTimer]);

  // エラーハンドラー（非同期化してレンダリング中の状態更新を回避）
  const handleError = useCallback(
    (message: string) => {
      console.error('Ketcher error:', message);

      // レンダリングサイクルの外で状態更新を実行
      queueMicrotask(() => {
        setConvertError(message);
        useNotificationStore.getState().addNotification({
          type: 'error',
          title: 'Ketcher Error',
          message: message,
          autoClose: false,
          duration: 0,
        });
      });
    },
    []
  );

  // SMILESをXYZに変換してCalculation Settingsページへ遷移
  const handleConvertToXyz = async () => {
    if (!ketcherInstanceRef.current) {
      setConvertError('Ketcher editor is not initialized');
      useNotificationStore.getState().addNotification({
        type: 'error',
        title: 'Editor Not Ready',
        message: 'Ketcher editor is not initialized',
        autoClose: false,
        duration: 0,
      });
      return;
    }

    setIsConverting(true);
    setConvertError(null);

    try {
      // KetcherからSMILES形式で構造を取得
      const smiles = await ketcherInstanceRef.current.getSmiles();

      if (!smiles || smiles.trim() === '') {
        throw new Error('No molecule drawn. Please draw a molecule first.');
      }

      const serialized = await serializeForPersistence(
        ketcherInstanceRef.current,
        smiles
      );

      // SMILES → XYZ変換APIを呼び出し
      const response = await convertSmilesToXyz(smiles);

      if (response.xyz) {
        // 完了済みまたはエラー状態の計算を編集した場合は新規計算として扱う
        const isExistingCompleted =
          activeCalculation &&
          (activeCalculation.status === 'completed' ||
            activeCalculation.status === 'error');

        // 新しいIDを生成（新規計算または既存完了計算の編集の場合）
        const newId = isExistingCompleted
          ? `new-calculation-${crypto.randomUUID()}`
          : activeCalculation?.id || `new-calculation-${crypto.randomUUID()}`;

        const moleculeName = smiles.substring(0, 50);
        const calculationName = `Drawn Molecule (${moleculeName}${smiles.length > 50 ? '...' : ''})`;

        // 既存の計算パラメータを引き継ぐか、デフォルト値を使用
        const baseParams =
          activeCalculation?.parameters || DEFAULT_CALCULATION_PARAMETERS;

        const newCalculation = {
          id: newId,
          name: calculationName,
          status: 'pending' as const,
          createdAt: new Date().toISOString(),
          updatedAt: new Date().toISOString(),
          parameters: {
            ...baseParams,
            xyz: response.xyz,
            ketcher_data: serialized.data,
            name: calculationName,
          },
          results: undefined,
        };

        // Staged Calculationを設定
        setStagedCalculation(newCalculation);
        setActiveCalculationId(newId);

        // 復元状態をリセット（次回の遷移で復元できるように）
        hasRestoredRef.current = null;

        // Calculation Settingsページへ遷移
        setCurrentPage('calculation-settings');

        useNotificationStore.getState().addNotification({
          type: 'success',
          title: 'Success',
          message: isExistingCompleted
            ? 'New calculation created from edited molecule!'
            : 'Molecule converted successfully!',
          autoClose: true,
          duration: 3000,
        });
      } else {
        throw new Error('Failed to convert SMILES to XYZ');
      }
    } catch (error: any) {
      const errorMessage =
        error.message || 'An error occurred during conversion';
      setConvertError(errorMessage);
      useNotificationStore.getState().addNotification({
        type: 'error',
        title: 'Conversion Error',
        message: errorMessage,
        autoClose: false,
        duration: 0,
      });
      console.error('Conversion error:', error);
    } finally {
      setIsConverting(false);
    }
  };

  const statusMessage = activeCalculation
    ? STATUS_MESSAGES[activeCalculation.status] ?? null
    : null;

  return (
    <div className={styles.pageContainer}>
      {/* ヘッダーとボタンを横並びに */}
      <div className={styles.headerRow}>
        <div className={styles.pageHeader}>
          <h2 className={styles.pageTitle}>Draw Molecule</h2>
        </div>

        {/* アクションボタン */}
        {!isNonDrawMoleculeCalculation && (
          <div className={styles.actionsContainer}>
            <button
              className={styles.convertButton}
              onClick={handleConvertToXyz}
              disabled={isConverting || !canEdit}
            >
              {isConverting ? 'Converting...' : 'Convert to XYZ & Continue'}
            </button>
          </div>
        )}
      </div>

      <div className={styles.pageContent}>
        {/* ステータスバナー */}
        {statusMessage && (
          <div
            className={`${styles.statusBanner} ${
              activeCalculation?.status ? styles[activeCalculation.status] : ''
            }`}
          >
            <svg
              width="20"
              height="20"
              viewBox="0 0 20 20"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
              className={styles.statusIcon}
            >
              {activeCalculation?.status === 'running' ||
              activeCalculation?.status === 'waiting' ? (
                <path
                  d="M10 2C5.58172 2 2 5.58172 2 10C2 14.4183 5.58172 18 10 18C14.4183 18 18 14.4183 18 10"
                  stroke="currentColor"
                  strokeWidth="2"
                  strokeLinecap="round"
                />
              ) : activeCalculation?.status === 'completed' ? (
                <>
                  <path
                    d="M10 18C14.4183 18 18 14.4183 18 10C18 5.58172 14.4183 2 10 2C5.58172 2 2 5.58172 2 10C2 14.4183 5.58172 18 10 18Z"
                    stroke="currentColor"
                    strokeWidth="2"
                  />
                  <path
                    d="M6 10L9 13L14 7"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </>
              ) : (
                <>
                  <path
                    d="M10 18C14.4183 18 18 14.4183 18 10C18 5.58172 14.4183 2 10 2C5.58172 2 2 5.58172 2 10C2 14.4183 5.58172 18 10 18Z"
                    stroke="currentColor"
                    strokeWidth="2"
                  />
                  <path
                    d="M10 6V10"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                  />
                  <circle cx="10" cy="14" r="1" fill="currentColor" />
                </>
              )}
            </svg>
            <span className={styles.statusText}>{statusMessage}</span>
          </div>
        )}

        {/* Ketcher エディタ */}
        {isNonDrawMoleculeCalculation ? (
          <div className={styles.unavailableContainer}>
            <svg
              width="48"
              height="48"
              viewBox="0 0 24 24"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
              className={styles.unavailableIcon}
            >
              <path
                d="M12 22C17.5228 22 22 17.5228 22 12C22 6.47715 17.5228 2 12 2C6.47715 2 2 6.47715 2 12C2 17.5228 6.47715 22 12 22Z"
                stroke="currentColor"
                strokeWidth="2"
              />
              <path
                d="M12 8V12"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
              />
              <circle cx="12" cy="16" r="1" fill="currentColor" />
            </svg>
            <h3 className={styles.unavailableTitle}>Draw Molecule Unavailable</h3>
            <p className={styles.unavailableText}>
              This calculation was not started using Draw Molecule. The molecular
              editor is only available for calculations created through this
              page.
            </p>
          </div>
        ) : (
          <KetcherEditorErrorBoundary onRetry={handleRetryEditor}>
            <div
              className={`${styles.editorContainer} ${!canEdit ? styles.readOnly : ''}`}
            >
              <KetcherEditor
                key={editorRetryKey}
                errorHandler={handleError}
                onInit={handleOnInit}
              />
            </div>
          </KetcherEditorErrorBoundary>
        )}

        {/* エラー表示 */}
        {convertError && (
          <div className={styles.errorContainer}>
            <svg
              width="20"
              height="20"
              viewBox="0 0 20 20"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
              className={styles.errorIcon}
            >
              <path
                d="M10 18C14.4183 18 18 14.4183 18 10C18 5.58172 14.4183 2 10 2C5.58172 2 2 5.58172 2 10C2 14.4183 5.58172 18 10 18Z"
                stroke="currentColor"
                strokeWidth="2"
              />
              <path
                d="M10 6V10"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
              />
              <circle cx="10" cy="14" r="1" fill="currentColor" />
            </svg>
            <span className={styles.errorText}>{convertError}</span>
          </div>
        )}
      </div>
    </div>
  );
};
