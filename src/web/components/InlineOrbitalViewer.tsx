// src/web/components/InlineOrbitalViewer.tsx

import React, { useEffect, useRef, useState, useCallback } from 'react';
import {
  useGetOrbitals,
  useGetOrbitalCube,
} from '../hooks/useCalculationQueries';
import { OrbitalInfo } from '../types/api-types';
import styles from './InlineOrbitalViewer.module.css';
import * as $3Dmol from '3dmol';
import { GLViewer } from '../../types/3dmol';

interface InlineOrbitalViewerProps {
  calculation_id: string;
  orbital_index: number;
  grid_size?: number;
  isovalue_pos?: number;
  isovalue_neg?: number;
  onError?: (error: string) => void;
}

interface ViewerOptions {
  gridSize: number;
  isovaluePos: number;
  isovalueNeg: number;
}

export const InlineOrbitalViewer: React.FC<InlineOrbitalViewerProps> =
  React.memo(
    ({
      calculation_id,
      orbital_index,
      grid_size = 80,
      isovalue_pos = 0.02,
      isovalue_neg = -0.02,
      onError,
    }) => {
      const viewerRef = useRef<HTMLDivElement>(null);
      const resizeObserverRef = useRef<ResizeObserver | null>(null);

      // State declarations
      const [viewer, setViewer] = useState<GLViewer | null>(null);
      const [viewerOptions, setViewerOptions] = useState<ViewerOptions>({
        gridSize: grid_size,
        isovaluePos: isovalue_pos,
        isovalueNeg: isovalue_neg,
      });
      const [isLoading, setIsLoading] = useState(false);
      const [retryCount, setRetryCount] = useState(0);
      const [isDomReady, setIsDomReady] = useState(false);

      // Callback ref to properly track when DOM element becomes available
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
      } = useGetOrbitals(calculation_id);

      // 選択された軌道のCUBEファイルを取得
      const {
        data: cubeData,
        isLoading: cubeLoading,
        error: cubeError,
      } = useGetOrbitalCube(calculation_id, orbital_index, viewerOptions);

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

      // 3Dmol.jsビューアーの初期化
      useEffect(() => {
        if (!isDomReady || !viewerRef.current || viewer) {
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
                setViewer(newViewer);
                setRetryCount(0);

                // ビューアーのリサイズとレンダリングを実行
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
              });
            });
          } catch (error) {
            console.error('Failed to initialize 3Dmol viewer:', error);

            if (retryCount < 3) {
              setRetryCount(prev => prev + 1);
              timeoutId = setTimeout(
                () => {
                  animationFrameId = requestAnimationFrame(initializeViewer);
                },
                500 * (retryCount + 1)
              );
            } else {
              onError?.(
                'Failed to initialize molecular viewer after multiple attempts'
              );
            }
          }
        };

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
      }, [onError, retryCount, isDomReady]);

      // ResizeObserverのセットアップ
      useEffect(() => {
        if (!viewerRef.current) {
          return;
        }

        resizeObserverRef.current = new ResizeObserver(entries => {
          for (const entry of entries) {
            if (entry.target === viewerRef.current) {
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

      // CUBEデータが更新されたときに分子軌道を表示
      useEffect(() => {
        if (!viewer || !cubeData || !(cubeData as any).cube_data) {
          return;
        }

        setIsLoading(true);

        // 非同期でレンダリング処理を実行
        const renderOrbital = async () => {
          try {
            // ビューアーをクリアして、WebGLリソースの解放を待機
            viewer.clear();
            await new Promise(resolve => setTimeout(resolve, 50));

            const cubeContent = (cubeData as any).cube_data;

            console.log('Rendering orbital visualization:', {
              calculation_id,
              orbital_index,
              hasContent: !!cubeContent,
            });

            if (
              !cubeContent ||
              typeof cubeContent !== 'string' ||
              cubeContent.length < 100
            ) {
              throw new Error('Invalid CUBE data');
            }

            // 分子構造を追加
            viewer.addModel(cubeContent, 'cube');
            viewer.setStyle(
              {},
              { stick: { radius: 0.1 }, sphere: { radius: 0.3 } }
            );

            // 正の等値面（赤色）
            viewer.addVolumetricData(cubeContent, 'cube', {
              isoval: viewerOptions.isovaluePos,
              color: 'red',
              opacity: 0.75,
            });

            // 負の等値面（青色）
            viewer.addVolumetricData(cubeContent, 'cube', {
              isoval: viewerOptions.isovalueNeg,
              color: 'blue',
              opacity: 0.75,
            });

            viewer.zoomTo();
            viewer.render();

            // レンダリング完了を待ってから再度リサイズ
            await new Promise(resolve => setTimeout(resolve, 200));

            try {
              if (viewer && viewerRef.current) {
                viewer.resize();
                viewer.zoomTo();
                viewer.render();
              }
            } catch (delayedRenderError) {
              console.error(
                'Failed to perform delayed re-render:',
                delayedRenderError
              );
            }

            setIsLoading(false);
          } catch (error) {
            console.error('Failed to render molecular orbital:', error);
            onError?.(
              `Failed to render molecular orbital: ${error instanceof Error ? error.message : String(error)}`
            );
            setIsLoading(false);
          }
        };

        renderOrbital();
      }, [
        viewer,
        cubeData,
        viewerOptions,
        calculation_id,
        orbital_index,
        onError,
      ]);

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
          console.error(
            'Failed to load orbital visualization data:',
            cubeError
          );
          onError?.(
            cubeError.message || 'Failed to load orbital visualization data'
          );
        }
      }, [cubeError, onError]);

      if (orbitalsLoading) {
        return (
          <div className={styles.loadingContainer}>
            <div className={styles.spinner}></div>
            <div>Loading orbital information...</div>
          </div>
        );
      }

      if (orbitalsError || !orbitalsData) {
        return (
          <div className={styles.errorContainer}>
            <div>❌ Failed to load orbital data</div>
            {orbitalsError && (
              <div className={styles.errorDetail}>{orbitalsError.message}</div>
            )}
          </div>
        );
      }

      const selectedOrbital = orbitalsData.orbitals.find(
        (orbital: OrbitalInfo) => orbital.index === orbital_index
      );

      if (!selectedOrbital) {
        return (
          <div className={styles.errorContainer}>
            <div>❌ Orbital index {orbital_index} not found</div>
          </div>
        );
      }

      return (
        <div className={styles.viewerContainer}>
          {/* 軌道情報ヘッダー */}
          <div className={styles.orbitalHeader}>
            <div className={styles.orbitalTitle}>
              <strong>{selectedOrbital.label}</strong> (Orbital #
              {selectedOrbital.index})
            </div>
            <div className={styles.orbitalInfo}>
              Energy: <code>{selectedOrbital.energy_ev.toFixed(4)} eV</code> (
              <code>{selectedOrbital.energy_hartree.toFixed(6)} a.u.</code>)
              {' | '}
              Occupancy: <code>{selectedOrbital.occupancy}</code>
            </div>
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
                <div className={styles.spinner}></div>
                <div>Generating orbital data...</div>
              </div>
            )}
          </div>

          {/* 簡易コントロール */}
          <div className={styles.simpleControls}>
            <div className={styles.controlRow}>
              <label>
                Positive Isovalue: {viewerOptions.isovaluePos.toFixed(3)}
              </label>
              <input
                type="range"
                min="0.001"
                max="0.1"
                step="0.001"
                value={viewerOptions.isovaluePos}
                onChange={e =>
                  setViewerOptions(prev => ({
                    ...prev,
                    isovaluePos: parseFloat(e.target.value),
                  }))
                }
              />
            </div>
            <div className={styles.controlRow}>
              <label>
                Negative Isovalue: {viewerOptions.isovalueNeg.toFixed(3)}
              </label>
              <input
                type="range"
                min="-0.1"
                max="-0.001"
                step="0.001"
                value={viewerOptions.isovalueNeg}
                onChange={e =>
                  setViewerOptions(prev => ({
                    ...prev,
                    isovalueNeg: parseFloat(e.target.value),
                  }))
                }
              />
            </div>
          </div>
        </div>
      );
    }
  );
