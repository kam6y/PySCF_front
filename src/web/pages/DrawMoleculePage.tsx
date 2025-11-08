// src/web/pages/DrawMoleculePage.tsx

import React, { useState, useRef, useEffect, useCallback, memo } from 'react';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Ketcher } from 'ketcher-core';
import 'ketcher-react/dist/index.css';
import styles from './DrawMoleculePage.module.css';
import { convertSmilesToXyz } from '../apiClient';
import { useUIStore } from '../store/uiStore';
import { useCalculationStore } from '../store/calculationStore';
import { useNotificationStore } from '../store/notificationStore';
import { useActiveCalculation } from '../hooks/useActiveCalculation';

// Miewをwindowに設定（Ketcherが3D表示に使用）
if (typeof window !== 'undefined') {
  import('miew').then(Miew => {
    (window as any).Miew = Miew.default || Miew;
  });
}

// StandaloneStructServiceProviderのインスタンスを作成
const structServiceProvider = new StandaloneStructServiceProvider();

// Ketcherエディタコンポーネント（メモ化して不要な再レンダリングを防ぐ）
const KetcherEditor = memo<{
  errorHandler: (message: string) => void;
  onInit: (ketcher: Ketcher) => void;
}>(({ errorHandler, onInit }) => {
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
  const hasRestoredRef = useRef<string | null>(null); // 復元済みのcalculation ID

  // Zustandストア
  const setCurrentPage = useUIStore(state => state.setCurrentPage);
  const setStagedCalculation = useCalculationStore(
    state => state.setStagedCalculation
  );
  const setActiveCalculationId = useCalculationStore(
    state => state.setActiveCalculationId
  );
  const addNotification = useNotificationStore(state => state.addNotification);

  // アクティブな計算を取得
  const { activeCalculation } = useActiveCalculation();

  // 編集可否の判定（running/waitingの場合は編集不可）
  const canEdit =
    !activeCalculation ||
    activeCalculation.status === 'pending' ||
    activeCalculation.status === 'error' ||
    activeCalculation.status === 'completed';

  // ページマウント時のログ
  useEffect(() => {
    console.log('[DrawMoleculePage] Component mounted/updated', {
      activeCalculationId: activeCalculation?.id,
      hasKetcherData: !!activeCalculation?.parameters?.ketcher_data,
      ketcherDataLength: activeCalculation?.parameters?.ketcher_data?.length,
      status: activeCalculation?.status,
    });
  }, [activeCalculation?.id]);

  // Ketcherインスタンスの初期化（メモ化して安定化）
  const handleOnInit = useCallback((ketcher: Ketcher) => {
    ketcherInstanceRef.current = ketcher;

    // デバッグ用にwindowにも設定
    if (typeof window !== 'undefined') {
      (window as any).ketcher = ketcher;
    }

    // 少し遅延してから復元（複数回の初期化に対応）
    setTimeout(() => {
      setIsKetcherReady(true);
    }, 100);
  }, []);

  // activeCalculation変更時にKetcherデータを復元
  useEffect(() => {
    const restoreKetcherData = async () => {
      if (!ketcherInstanceRef.current || !isKetcherReady) {
        return;
      }

      const currentCalcId = activeCalculation?.id;
      const ketcherData = activeCalculation?.parameters?.ketcher_data;

      // 既に同じcalculationを復元済みの場合はスキップ
      if (hasRestoredRef.current === currentCalcId) {
        return;
      }

      if (ketcherData) {
        try {
          // Ketcher JSONデータから分子を復元
          await ketcherInstanceRef.current.setMolecule(ketcherData);
          hasRestoredRef.current = currentCalcId || null;
        } catch (error) {
          console.error(
            '[Ketcher Restore] Failed to restore ketcher data:',
            error
          );
          addNotification({
            type: 'error',
            title: 'Restore Error',
            message: 'Failed to restore molecule structure',
            autoClose: true,
            duration: 3000,
          });
        }
      } else {
        // ketcher_dataがない場合はクリア
        try {
          await ketcherInstanceRef.current.setMolecule('');
          hasRestoredRef.current = currentCalcId || null;
        } catch (error) {
          console.error('[Ketcher Restore] Failed to clear ketcher:', error);
        }
      }
    };

    restoreKetcherData();
  }, [
    activeCalculation?.id,
    isKetcherReady,
    activeCalculation?.parameters?.ketcher_data,
    addNotification,
  ]);

  // エラーハンドラー（非同期化してレンダリング中の状態更新を回避）
  const handleError = useCallback(
    (message: string) => {
      console.error('Ketcher error:', message);

      // レンダリングサイクルの外で状態更新を実行
      queueMicrotask(() => {
        setConvertError(message);
        addNotification({
          type: 'error',
          title: 'Ketcher Error',
          message: message,
          autoClose: false,
          duration: 0,
        });
      });
    },
    [addNotification]
  );

  // SMILESをXYZに変換してCalculation Settingsページへ遷移
  const handleConvertToXyz = async () => {
    if (!ketcherInstanceRef.current) {
      setConvertError('Ketcher editor is not initialized');
      addNotification({
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

      // Ketcher JSONデータも取得
      const ketcherData = await ketcherInstanceRef.current.getKet();

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
          ? `new-calculation-${Date.now()}`
          : activeCalculation?.id || `new-calculation-${Date.now()}`;

        const moleculeName = smiles.substring(0, 50);
        const calculationName = `Drawn Molecule (${moleculeName}${smiles.length > 50 ? '...' : ''})`;

        // 既存の計算パラメータを引き継ぐか、デフォルト値を使用
        const baseParams = activeCalculation?.parameters || {
          calculation_method: 'DFT' as const,
          basis_function: '6-31G(d)',
          exchange_correlation: 'B3LYP',
          charges: 0,
          spin: 0,
          solvent_method: 'none' as const,
          solvent: '-',
          tddft_nstates: 10,
          tddft_method: 'TDDFT' as const,
          tddft_analyze_nto: false,
          ncas: 4,
          nelecas: 4,
          max_cycle_macro: 50,
          max_cycle_micro: 4,
          natorb: true,
          conv_tol: 1e-6,
          conv_tol_grad: 1e-4,
          optimize_geometry: true,
        };

        const newCalculation = {
          id: newId,
          name: calculationName,
          status: 'pending' as const,
          createdAt: new Date().toISOString(),
          updatedAt: new Date().toISOString(),
          parameters: {
            ...baseParams,
            xyz: response.xyz,
            ketcher_data: ketcherData, // Ketcher JSONデータを保存
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

        addNotification({
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
      addNotification({
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

  // 分子をクリア
  const handleClearMolecule = async () => {
    if (!ketcherInstanceRef.current) return;

    try {
      await ketcherInstanceRef.current.setMolecule('');
      setConvertError(null);
      addNotification({
        type: 'info',
        title: 'Cleared',
        message: 'Molecule cleared',
        autoClose: true,
        duration: 2000,
      });
    } catch (error) {
      console.error('Failed to clear molecule:', error);
    }
  };

  // ステータスメッセージの生成
  const getStatusMessage = () => {
    if (!activeCalculation) return null;

    switch (activeCalculation.status) {
      case 'running':
        return 'Calculation is running. The molecule structure is read-only.';
      case 'waiting':
        return 'Calculation is waiting for resources. The molecule structure is read-only.';
      case 'completed':
        return 'Calculation completed. Edit the structure to create a new calculation.';
      case 'error':
        return 'Previous calculation had errors. Edit the structure to create a new calculation.';
      default:
        return null;
    }
  };

  const statusMessage = getStatusMessage();

  return (
    <div className={styles.pageContainer}>
      {/* ヘッダーとボタンを横並びに */}
      <div className={styles.headerRow}>
        <div className={styles.pageHeader}>
          <h2 className={styles.pageTitle}>Draw Molecule</h2>
          <p className={styles.pageDescription}>
            Draw your molecule structure and convert it to XYZ coordinates for
            quantum calculations
          </p>
        </div>

        {/* アクションボタン */}
        <div className={styles.actionsContainer}>
          <button
            className={styles.clearButton}
            onClick={handleClearMolecule}
            disabled={isConverting || !canEdit}
          >
            Clear Molecule
          </button>
          <button
            className={styles.convertButton}
            onClick={handleConvertToXyz}
            disabled={isConverting || !canEdit}
          >
            {isConverting ? 'Converting...' : 'Convert to XYZ & Continue'}
          </button>
        </div>
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
        <div
          className={`${styles.editorContainer} ${!canEdit ? styles.readOnly : ''}`}
        >
          <KetcherEditor errorHandler={handleError} onInit={handleOnInit} />
        </div>

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
