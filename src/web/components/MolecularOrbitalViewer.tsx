// src/web/components/MolecularOrbitalViewer.tsx

import React, { useEffect, useRef, useState, useCallback } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import {
  useGetOrbitals,
  useGetOrbitalCube,
} from '../hooks/useCalculationQueries';
import { OrbitalInfo } from '../types/api-types';
import styles from './MolecularOrbitalViewer.module.css';
import * as $3Dmol from '3dmol';
import { GLViewer } from '../../types/3dmol';

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
    const [viewer, setViewer] = useState<GLViewer | null>(null);
    const [selectedOrbitalIndex, setSelectedOrbitalIndex] = useState<
      number | null
    >(null);
    const [viewerOptions, setViewerOptions] = useState<ViewerOptions>({
      gridSize: 80,
      isovaluePos: 0.02,
      isovalueNeg: -0.02,
    });
    const [isLoading, setIsLoading] = useState(false);
    const retryCountRef = useRef(0);
    const [isDomReady, setIsDomReady] = useState(false);
    const [isViewerReady, setIsViewerReady] = useState(false);

    const previousCalculationIdRef = useRef<string | null>(null);

    // Callback ref to track when DOM element becomes available
    const setViewerRef = useCallback(
      (node: HTMLDivElement | null) => {
        if (viewerRef.current !== node) {
          viewerRef.current = node;
          const isReady = !!node;
          if (isReady !== isDomReady) {
            setIsDomReady(isReady);
          }
        }
      },
      [isDomReady]
    );

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
      let viewerInstance: GLViewer | null = null;

      const initializeViewer = () => {
        if (!viewerRef.current) {
          return;
        }

        try {
          const newViewer = $3Dmol.createViewer(viewerRef.current, {
            backgroundColor: 'white',
            antialias: true,
          }) as unknown as GLViewer;

          viewerInstance = newViewer;

          // WebGLコンテキストの準備完了を確実に待機（2フレーム待機）
          requestAnimationFrame(() => {
            requestAnimationFrame(() => {
              // ビューアーをすぐにセット（InlineOrbitalViewerと同じパターン）
              setViewer(newViewer);
              retryCountRef.current = 0;

              // パッケージ環境でのGPU初期化を考慮して待機時間を延長
              timeoutId = setTimeout(() => {
                if (newViewer && viewerRef.current) {
                  try {
                    // ビューアーの基本的な初期化を実行
                    newViewer.resize();
                    newViewer.render();

                    // ビューアーの準備完了を明示的にマーク
                    setIsViewerReady(true);
                  } catch (error) {
                    console.error(
                      '[MolecularOrbitalViewer] Failed to finalize viewer initialization:',
                      error
                    );
                    // 初期化失敗時はビューアーをクリア
                    setViewer(null);
                    setIsViewerReady(false);

                    if (retryCountRef.current < 3) {
                      retryCountRef.current += 1;
                    } else {
                      onError?.(
                        'Failed to initialize molecular viewer after multiple attempts'
                      );
                    }
                  }
                }
              }, 200); // パッケージ環境対応のため200ms
            });
          });
        } catch (error) {
          console.error('Failed to initialize 3Dmol viewer:', error);

          // 再試行ロジック（最大3回まで）
          if (retryCountRef.current < 3) {
            retryCountRef.current += 1;
            timeoutId = setTimeout(() => {
              animationFrameId = requestAnimationFrame(initializeViewer);
            }, 500 * retryCountRef.current); // 指数バックオフ
          } else {
            console.error(
              'Failed to initialize molecular viewer after 3 attempts'
            );
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
        if (viewerInstance) {
          viewerInstance.clear();
        }
      };
    }, [onError, isDomReady]);

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
        setIsViewerReady(false);
        retryCountRef.current = 0;

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
      // ビューアーが完全に準備完了するまで待機
      if (
        !viewer ||
        !isViewerReady ||
        !cubeData ||
        !(cubeData as any).cube_data
      ) {
        return;
      }

      setIsLoading(true);

      // 非同期でレンダリング処理を実行
      const renderOrbital = async () => {
        try {
          // ビューアーが有効であることを再確認
          if (!viewer || !isViewerReady) {
            setIsLoading(false);
            return;
          }

          // ビューアーをクリアして、WebGLリソースの解放を待機
          viewer.clear();
          await new Promise(resolve => setTimeout(resolve, 50));

          // CUBEファイルデータを追加
          const cubeContent = (cubeData as any).cube_data;

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

          // 分子構造を先に追加（等値面の前に）
          try {
            viewer.addModel(cubeContent, 'cube');
            viewer.setStyle(
              {},
              { stick: { radius: 0.1 }, sphere: { radius: 0.3 } }
            );
          } catch (modelError) {
            console.error('Failed to add molecular model:', modelError);
            // ビューアーをクリーンアップ
            try {
              viewer.clear();
            } catch (clearError) {
              console.error(
                'Failed to clear viewer after model error:',
                clearError
              );
            }
            throw new Error('Failed to add molecular structure from CUBE data');
          }

          // 正の等値面（赤色）- 保護された呼び出し
          try {
            viewer.addVolumetricData(cubeContent, 'cube', {
              isoval: viewerOptions.isovaluePos,
              color: 'red',
              opacity: 0.75,
            });
          } catch (posError) {
            console.error('Error adding positive isosurface:', posError);
            // ビューアーをクリーンアップ
            try {
              viewer.clear();
            } catch (clearError) {
              console.error(
                'Failed to clear viewer after positive isosurface error:',
                clearError
              );
            }
            throw new Error(
              `Failed to add positive isosurface: ${posError instanceof Error ? posError.message : String(posError)}`
            );
          }

          // 負の等値面（青色）- 保護された呼び出し
          try {
            viewer.addVolumetricData(cubeContent, 'cube', {
              isoval: viewerOptions.isovalueNeg,
              color: 'blue',
              opacity: 0.75,
            });
          } catch (negError) {
            console.error('Error adding negative isosurface:', negError);
            // ビューアーをクリーンアップ
            try {
              viewer.clear();
            } catch (clearError) {
              console.error(
                'Failed to clear viewer after negative isosurface error:',
                clearError
              );
            }
            throw new Error(
              `Failed to add negative isosurface: ${negError instanceof Error ? negError.message : String(negError)}`
            );
          }

          // ビューを最適化（即座に実行）
          try {
            viewer.zoomTo();
            viewer.render();
          } catch (renderError) {
            console.error('Failed to perform initial render:', renderError);
            throw new Error('Failed to render orbital visualization');
          }

          // レンダリング完了を待ってから再度リサイズ
          await new Promise(resolve => setTimeout(resolve, 200));

          // 最終的な調整とレンダリング
          try {
            if (viewer && viewerRef.current && isViewerReady) {
              viewer.resize();
              viewer.zoomTo();
              viewer.render();
            }
          } catch (delayedRenderError) {
            console.error(
              'Failed to perform delayed re-render:',
              delayedRenderError
            );
            // 最終レンダリングの失敗は致命的ではない（警告のみ）
          }

          setIsLoading(false);
        } catch (error) {
          console.error('Failed to render molecular orbital:', error);
          onError?.(
            `Failed to render molecular orbital visualization: ${error instanceof Error ? error.message : String(error)}`
          );
          setIsLoading(false);
        }
      };

      renderOrbital();
    }, [viewer, isViewerReady, cubeData, viewerOptions, onError]);

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
          <div>Loading orbital information...</div>
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
          <div>No orbital data available for this calculation.</div>
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
              Generating orbital data...
            </div>
          )}
        </div>
      </div>
    );
  });
