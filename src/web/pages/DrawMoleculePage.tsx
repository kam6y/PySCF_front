// src/web/pages/DrawMoleculePage.tsx

import React, { useState, useRef, useEffect } from 'react';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Ketcher } from 'ketcher-core';
import 'ketcher-react/dist/index.css';
import 'miew/dist/miew.min.css';
import styles from './DrawMoleculePage.module.css';
import { convertSmilesToXyz } from '../apiClient';
import { useUIStore } from '../store/uiStore';
import { useCalculationStore } from '../store/calculationStore';
import { useNotificationStore } from '../store/notificationStore';

// Miewをwindowに設定（Ketcherが3D表示に使用）
if (typeof window !== 'undefined') {
  import('miew').then((Miew) => {
    (window as any).Miew = Miew.default || Miew;
  });
}

// StandaloneStructServiceProviderのインスタンスを作成
const structServiceProvider = new StandaloneStructServiceProvider();

export const DrawMoleculePage: React.FC = () => {
  const ketcherInstanceRef = useRef<Ketcher | null>(null);
  const [isConverting, setIsConverting] = useState(false);
  const [convertError, setConvertError] = useState<string | null>(null);

  // Zustandストア
  const setCurrentPage = useUIStore((state) => state.setCurrentPage);
  const setStagedCalculation = useCalculationStore(
    (state) => state.setStagedCalculation
  );
  const setActiveCalculationId = useCalculationStore(
    (state) => state.setActiveCalculationId
  );
  const addNotification = useNotificationStore((state) => state.addNotification);

  // Ketcherインスタンスの初期化
  const handleOnInit = (ketcher: Ketcher) => {
    ketcherInstanceRef.current = ketcher;
    // デバッグ用にwindowにも設定
    if (typeof window !== 'undefined') {
      (window as any).ketcher = ketcher;
    }
    console.log('Ketcher initialized successfully');
  };

  // エラーハンドラー
  const handleError = (message: string) => {
    console.error('Ketcher error:', message);
    setConvertError(message);
    addNotification({
      type: 'error',
      title: 'Ketcher Error',
      message: message,
      autoClose: false,
      duration: 0,
    });
  };

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

      console.log('Generated SMILES:', smiles);

      // SMILES → XYZ変換APIを呼び出し
      const response = await convertSmilesToXyz(smiles);

      if (response.xyz) {
        // 新しいStagedCalculationを作成
        const newId = `new-calculation-${Date.now()}`;
        const moleculeName = smiles.substring(0, 50); // Use SMILES string as name
        const newCalculation = {
          id: newId,
          name: `Drawn Molecule (${moleculeName}${smiles.length > 50 ? '...' : ''})`,
          status: 'pending' as const,
          createdAt: new Date().toISOString(),
          updatedAt: new Date().toISOString(),
          parameters: {
            calculation_method: 'DFT' as const,
            basis_function: '6-31G(d)',
            exchange_correlation: 'B3LYP',
            charges: 0,
            spin: 0,
            solvent_method: 'none' as const,
            solvent: '-',
            xyz: response.xyz,
            name: `Drawn Molecule (${moleculeName}${smiles.length > 50 ? '...' : ''})`,
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
          },
          results: undefined,
        };

        // Staged Calculationを設定
        setStagedCalculation(newCalculation);
        setActiveCalculationId(newId);

        // Calculation Settingsページへ遷移
        setCurrentPage('calculation-settings');

        addNotification({
          type: 'success',
          title: 'Success',
          message: 'Molecule converted successfully!',
          autoClose: true,
          duration: 3000,
        });
      } else {
        throw new Error('Failed to convert SMILES to XYZ');
      }
    } catch (error: any) {
      const errorMessage = error.message || 'An error occurred during conversion';
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

  return (
    <div className={styles.pageContainer}>
      {/* ヘッダーとボタンを横並びに */}
      <div className={styles.headerRow}>
        <div className={styles.pageHeader}>
          <h2 className={styles.pageTitle}>Draw Molecule</h2>
          <p className={styles.pageDescription}>
            Draw your molecule structure and convert it to XYZ coordinates for quantum calculations
          </p>
        </div>

        {/* アクションボタン */}
        <div className={styles.actionsContainer}>
          <button
            className={styles.clearButton}
            onClick={handleClearMolecule}
            disabled={isConverting}
          >
            Clear Molecule
          </button>
          <button
            className={styles.convertButton}
            onClick={handleConvertToXyz}
            disabled={isConverting}
          >
            {isConverting ? 'Converting...' : 'Convert to XYZ & Continue'}
          </button>
        </div>
      </div>

      <div className={styles.pageContent}>
        {/* Ketcher エディタ */}
        <div className={styles.editorContainer}>
          <Editor
            staticResourcesUrl=""
            structServiceProvider={structServiceProvider}
            errorHandler={handleError}
            onInit={handleOnInit}
          />
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
