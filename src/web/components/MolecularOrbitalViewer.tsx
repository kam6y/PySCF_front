// src/web/components/MolecularOrbitalViewer.tsx

import React, { useEffect, useRef, useState, useCallback } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import {
  useGetOrbitals,
  useGetOrbitalCube,
} from '../hooks/useCalculationQueries';
import { OrbitalInfo } from '../types/api-types';
import styles from './MolecularOrbitalViewer.module.css';

// 3Dmol.jsはCDNから読み込まれ、グローバル変数として利用可能
declare const $3Dmol: any;

interface MolecularOrbitalViewerProps {
  calculationId: string;
  onError?: (error: string) => void;
}

interface ViewerOptions {
  gridSize: number;
  isovaluePos: number;
  isovalueNeg: number;
}

export const MolecularOrbitalViewer: React.FC<MolecularOrbitalViewerProps> =
  React.memo(({ calculationId, onError }) => {
    const queryClient = useQueryClient();
    const viewerRef = useRef<HTMLDivElement>(null);
    const resizeObserverRef = useRef<ResizeObserver | null>(null);

    // State declarations
    const [viewer, setViewer] = useState<any>(null);
    const [selectedOrbitalIndex, setSelectedOrbitalIndex] = useState<
      number | null
    >(null);
    const [viewerOptions, setViewerOptions] = useState<ViewerOptions>({
      gridSize: 80,
      isovaluePos: 0.02,
      isovalueNeg: -0.02,
    });
    const [isLoading, setIsLoading] = useState(false);
    const [retryCount, setRetryCount] = useState(0);
    const [isDomReady, setIsDomReady] = useState(false);

    // Callback ref to properly track when DOM element becomes available
    const setViewerRef = useCallback(
      (node: HTMLDivElement | null) => {
        if (viewerRef.current !== node) {
          viewerRef.current = node;

          // Update DOM ready state when ref changes
          const isReady = !!node;
          if (isReady !== isDomReady) {
            setIsDomReady(isReady);
          }
        }
      },
      [isDomReady]
    );
    const previousCalculationIdRef = useRef<string | null>(null);

    // 軌道情報を取得
    const {
      data: orbitalsData,
      isLoading: orbitalsLoading,
      error: orbitalsError,
    } = useGetOrbitals(calculationId);

    // 選択された軌道のCUBEファイルを取得
    const {
      data: cubeData,
      isLoading: cubeLoading,
      error: cubeError,
    } = useGetOrbitalCube(calculationId, selectedOrbitalIndex, viewerOptions);

    // ビューアーのリサイズ処理
    const handleViewerResize = useCallback(() => {
      if (viewer && viewerRef.current) {
        try {
          viewer.resize();
          viewer.zoomTo();
          viewer.render();
        } catch (error) {
          console.error('Failed to resize viewer:', error);
        }
      }
    }, [viewer]);

    // 3Dmol.jsビューアーの初期化（遅延実行）
    useEffect(() => {
      if (!isDomReady || !viewerRef.current) {
        return;
      }

      if (viewer) {
        return;
      }

      let animationFrameId: number;
      let timeoutId: NodeJS.Timeout;

      const initializeViewer = () => {
        if (!viewerRef.current) {
          return;
        }

        try {
          const newViewer = $3Dmol.createViewer(viewerRef.current, {
            backgroundColor: 'white',
            antialias: true,
          });

          setViewer(newViewer);
          setRetryCount(0);

          // 初期化後にリサイズを実行
          timeoutId = setTimeout(() => {
            if (newViewer && viewerRef.current) {
              try {
                newViewer.resize();
                newViewer.render();
              } catch (error) {
                console.error(
                  'Failed to resize viewer after initialization:',
                  error
                );
              }
            }
          }, 100);
        } catch (error) {
          console.error('Failed to initialize 3Dmol viewer:', error);

          // 再試行ロジック（最大3回まで）
          if (retryCount < 3) {
            setRetryCount(prev => prev + 1);
            timeoutId = setTimeout(
              () => {
                animationFrameId = requestAnimationFrame(initializeViewer);
              },
              500 * (retryCount + 1)
            ); // 指数バックオフ
          } else {
            onError?.(
              'Failed to initialize molecular viewer after multiple attempts'
            );
          }
        }
      };

      // requestAnimationFrameで遅延実行
      animationFrameId = requestAnimationFrame(initializeViewer);

      return () => {
        if (animationFrameId) {
          cancelAnimationFrame(animationFrameId);
        }
        if (timeoutId) {
          clearTimeout(timeoutId);
        }
        if (viewer) {
          viewer.clear();
        }
      };
    }, [onError, retryCount, isDomReady, viewer]);

    // ResizeObserverのセットアップ
    useEffect(() => {
      if (!viewerRef.current) {
        return;
      }

      // ResizeObserver を設定してDOM要素のサイズ変更を監視
      resizeObserverRef.current = new ResizeObserver(entries => {
        for (const entry of entries) {
          if (entry.target === viewerRef.current) {
            // サイズ変更時にビューアーをリサイズ（遅延実行）
            setTimeout(() => {
              handleViewerResize();
            }, 50);
          }
        }
      });

      resizeObserverRef.current.observe(viewerRef.current);

      return () => {
        if (resizeObserverRef.current) {
          resizeObserverRef.current.disconnect();
        }
      };
    }, [handleViewerResize]);

    // calculationIdが変更されたときに状態をリセットとキャッシュ無効化
    useEffect(() => {
      const previousCalculationId = previousCalculationIdRef.current;

      // 実際にcalculationIdが変更された場合のみ処理を実行
      if (previousCalculationId !== calculationId) {
        // 状態をリセット
        setSelectedOrbitalIndex(null);
        setIsLoading(false);
        setIsDomReady(false);

        // 3Dmol.jsビューアーをクリア
        if (viewer) {
          try {
            viewer.clear();
          } catch (error) {
            console.error('Failed to clear viewer:', error);
          }
        }
        setViewer(null);

        // 軌道関連のクエリキャッシュを無効化（新しいcalculationIdが有効な場合のみ）
        if (calculationId && !calculationId.startsWith('new-calculation-')) {
          queryClient.invalidateQueries({
            queryKey: ['orbitals', calculationId],
          });
          queryClient.invalidateQueries({
            queryKey: ['orbital-cube', calculationId],
          });
        }

        // 現在のcalculationIdを記録
        previousCalculationIdRef.current = calculationId;
      }
    }, [calculationId, queryClient, viewer]);

    // HOMOを初期選択として設定
    useEffect(() => {
      if (orbitalsData && selectedOrbitalIndex === null) {
        const homoOrbital = orbitalsData.orbitals.find(
          (orbital: OrbitalInfo) => orbital.orbital_type === 'homo'
        );
        if (homoOrbital) {
          setSelectedOrbitalIndex(homoOrbital.index);
        }
      }
    }, [orbitalsData, selectedOrbitalIndex, calculationId]);

    // CUBEデータが更新されたときに分子軌道を表示
    useEffect(() => {
      if (!viewer || !cubeData || !(cubeData as any).cube_data) {
        return;
      }

      setIsLoading(true);

      try {
        // ビューアーをクリア
        viewer.clear();

        // CUBEファイルデータを追加
        const cubeContent = (cubeData as any).cube_data;

        // デバッグ情報: CUBEデータの内容確認
        console.log('CUBE data validation:', {
          hasContent: !!cubeContent,
          contentType: typeof cubeContent,
          contentLength: cubeContent?.length || 0,
          firstChars: cubeContent?.substring?.(0, 100) || 'N/A',
          isovaluePos: viewerOptions.isovaluePos,
          isovalueNeg: viewerOptions.isovalueNeg,
        });

        // CUBEデータの詳細検証
        if (!cubeContent || typeof cubeContent !== 'string') {
          throw new Error(
            'Invalid CUBE data: content is empty or not a string'
          );
        }

        if (cubeContent.length < 100) {
          throw new Error('Invalid CUBE data: content too short');
        }

        // CUBE形式の基本的な検証（ヘッダー行の存在確認）
        const lines = cubeContent.split('\n');
        if (lines.length < 6) {
          throw new Error('Invalid CUBE data: insufficient header lines');
        }

        // 等値面設定の範囲チェック
        if (
          Math.abs(viewerOptions.isovaluePos) < 0.001 ||
          Math.abs(viewerOptions.isovaluePos) > 1.0
        ) {
          console.warn(
            'Positive isovalue may be out of optimal range:',
            viewerOptions.isovaluePos
          );
        }
        if (
          Math.abs(viewerOptions.isovalueNeg) < 0.001 ||
          Math.abs(viewerOptions.isovalueNeg) > 1.0
        ) {
          console.warn(
            'Negative isovalue may be out of optimal range:',
            viewerOptions.isovalueNeg
          );
        }

        // 分子構造を先に追加（等値面の前に）
        try {
          viewer.addModel(cubeContent, 'cube');
          viewer.setStyle(
            {},
            { stick: { radius: 0.1 }, sphere: { radius: 0.3 } }
          );
          console.log('Successfully added molecular structure');
        } catch (modelError) {
          console.error('Failed to add molecular model:', modelError);
          throw new Error('Failed to add molecular structure from CUBE data');
        }

        // 正の等値面（赤色）- 保護された呼び出し
        try {
          console.log(
            'Adding positive isosurface with value:',
            viewerOptions.isovaluePos
          );
          viewer.addVolumetricData(cubeContent, 'cube', {
            isoval: viewerOptions.isovaluePos,
            color: 'red',
            opacity: 0.75,
          });
          console.log('Successfully added positive isosurface');
        } catch (posError) {
          console.error('Error adding positive isosurface:', posError);
          throw new Error(
            `Failed to add positive isosurface: ${posError instanceof Error ? posError.message : String(posError)}`
          );
        }

        // 負の等値面（青色）- 保護された呼び出し
        try {
          console.log(
            'Adding negative isosurface with value:',
            viewerOptions.isovalueNeg
          );
          viewer.addVolumetricData(cubeContent, 'cube', {
            isoval: viewerOptions.isovalueNeg,
            color: 'blue',
            opacity: 0.75,
          });
          console.log('Successfully added negative isosurface');
        } catch (negError) {
          console.error('Error adding negative isosurface:', negError);
          throw new Error(
            `Failed to add negative isosurface: ${negError instanceof Error ? negError.message : String(negError)}`
          );
        }

        // ビューを最適化（即座に実行）
        viewer.zoomTo();
        viewer.render();

        // レンダリング完了後の遅延再描画（レイアウトが確定してから）
        setTimeout(() => {
          try {
            if (viewer && viewerRef.current) {
              viewer.resize();
              viewer.zoomTo();
              viewer.render();
              console.log('Completed delayed re-render');
            }
          } catch (delayedRenderError) {
            console.error(
              'Failed to perform delayed re-render:',
              delayedRenderError
            );
          }
        }, 200);

        setIsLoading(false);
        console.log('Molecular orbital rendering completed successfully');
      } catch (error) {
        console.error('Failed to render molecular orbital:', error);
        onError?.(
          `Failed to render molecular orbital visualization: ${error instanceof Error ? error.message : String(error)}`
        );
        setIsLoading(false);
      }
    }, [viewer, cubeData, viewerOptions, onError]);

    // エラーハンドリング
    useEffect(() => {
      if (orbitalsError) {
        console.error('Failed to load orbital information:', orbitalsError);
        onError?.(
          orbitalsError.message || 'Failed to load orbital information'
        );
      }
    }, [orbitalsError, onError]);

    useEffect(() => {
      if (cubeError) {
        console.error('Failed to load orbital visualization data:', cubeError);
        onError?.(
          cubeError.message || 'Failed to load orbital visualization data'
        );
      }
    }, [cubeError, onError]);

    const handleOrbitalChange = (
      event: React.ChangeEvent<HTMLSelectElement>
    ) => {
      const newIndex = parseInt(event.target.value, 10);
      setSelectedOrbitalIndex(newIndex);
    };

    const handleOptionsChange = (newOptions: Partial<ViewerOptions>) => {
      setViewerOptions(prev => ({ ...prev, ...newOptions }));
    };

    if (orbitalsLoading) {
      return (
        <div className={styles.loadingContainer}>
          <div>⚛️ Loading orbital information...</div>
        </div>
      );
    }

    if (
      !orbitalsData ||
      !orbitalsData.orbitals ||
      orbitalsData.orbitals.length === 0
    ) {
      return (
        <div className={styles.noDataContainer}>
          <div>📊 No orbital data available for this calculation.</div>
        </div>
      );
    }

    const selectedOrbital = orbitalsData.orbitals.find(
      (orbital: OrbitalInfo) => orbital.index === selectedOrbitalIndex
    );

    return (
      <div className={styles.viewerContainer}>
        {/* コントロールパネル */}
        <div className={styles.controlPanel}>
          <div className={styles.controlGrid}>
            {/* 軌道選択 */}
            <div>
              <label htmlFor="orbital-select" className={styles.controlLabel}>
                Molecular Orbital:
              </label>
              <select
                id="orbital-select"
                value={selectedOrbitalIndex || ''}
                onChange={handleOrbitalChange}
                className={styles.controlSelect}
              >
                {orbitalsData.orbitals.map((orbital: OrbitalInfo) => (
                  <option key={orbital.index} value={orbital.index}>
                    {orbital.label} ({orbital.energy_ev.toFixed(4)} eV)
                  </option>
                ))}
              </select>
            </div>

            {/* グリッドサイズ */}
            <div>
              <label
                htmlFor="grid-size"
                style={{
                  display: 'block',
                  marginBottom: '4px',
                  fontWeight: 'bold',
                }}
              >
                Grid Size: {viewerOptions.gridSize}
              </label>
              <input
                id="grid-size"
                type="range"
                min="40"
                max="120"
                step="10"
                value={viewerOptions.gridSize}
                onChange={e =>
                  handleOptionsChange({
                    gridSize: parseInt(e.target.value, 10),
                  })
                }
                style={{ width: '100%' }}
              />
            </div>

            {/* 正の等値面 */}
            <div className={styles.controlGroup}>
              <label
                htmlFor="isovalue-pos"
                className={`${styles.rangeLabel} ${styles.controlLabel}`}
              >
                Positive Isovalue: {viewerOptions.isovaluePos.toFixed(3)}
              </label>
              <input
                id="isovalue-pos"
                type="range"
                min="0.001"
                max="0.1"
                step="0.001"
                value={viewerOptions.isovaluePos}
                onChange={e =>
                  handleOptionsChange({
                    isovaluePos: parseFloat(e.target.value),
                  })
                }
                className={styles.rangeInput}
              />
            </div>

            {/* 負の等値面 */}
            <div className={styles.controlGroup}>
              <label
                htmlFor="isovalue-neg"
                className={`${styles.rangeLabel} ${styles.controlLabel}`}
              >
                Negative Isovalue: {viewerOptions.isovalueNeg.toFixed(3)}
              </label>
              <input
                id="isovalue-neg"
                type="range"
                min="-0.1"
                max="-0.001"
                step="0.001"
                value={viewerOptions.isovalueNeg}
                onChange={e =>
                  handleOptionsChange({
                    isovalueNeg: parseFloat(e.target.value),
                  })
                }
                className={styles.rangeInput}
              />
            </div>
          </div>

          {/* 選択された軌道の情報 */}
          {selectedOrbital && (
            <div className={styles.selectedOrbitalInfo}>
              <div>
                <strong>{selectedOrbital.label}</strong> (Orbital #
                {selectedOrbital.index})
                <br />
                Energy: <code>{selectedOrbital.energy_ev.toFixed(4)} eV</code> (
                <code>{selectedOrbital.energy_hartree.toFixed(6)} a.u.</code>)
                <br />
                Occupancy: <code>{selectedOrbital.occupancy}</code>
              </div>
            </div>
          )}
        </div>

        {/* 3Dmol.jsビューアー */}
        <div className={styles.viewer3D}>
          <div
            ref={setViewerRef}
            style={{
              width: '100%',
              height: '100%',
              backgroundColor: 'white',
            }}
          />

          {/* ローディングオーバーレイ */}
          {(isLoading || cubeLoading) && (
            <div className={styles.loadingOverlay}>
              ⚛️ Generating orbital data...
            </div>
          )}
        </div>
      </div>
    );
  });
