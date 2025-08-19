// src/web/components/MolecularOrbitalViewer.tsx

import React, { useEffect, useRef, useState, useCallback } from 'react';
import { useGetOrbitals, useGetOrbitalCube } from '../hooks/useCalculationQueries';
import { OrbitalInfo } from '../types/api-types';

// 3dmolライブラリをインポート
import * as $3Dmol from '3dmol';

interface MolecularOrbitalViewerProps {
  calculationId: string;
  onError?: (error: string) => void;
}

interface ViewerOptions {
  gridSize: number;
  isovaluePos: number;
  isovalueNeg: number;
}

export const MolecularOrbitalViewer: React.FC<MolecularOrbitalViewerProps> = ({
  calculationId,
  onError,
}) => {
  const viewerRef = useRef<HTMLDivElement>(null);
  const resizeObserverRef = useRef<ResizeObserver | null>(null);
  const [viewer, setViewer] = useState<any>(null);
  const [selectedOrbitalIndex, setSelectedOrbitalIndex] = useState<number | null>(null);
  const [viewerOptions, setViewerOptions] = useState<ViewerOptions>({
    gridSize: 80,
    isovaluePos: 0.02,
    isovalueNeg: -0.02,
  });
  const [isLoading, setIsLoading] = useState(false);
  const [retryCount, setRetryCount] = useState(0);

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
    if (!viewerRef.current) return;

    let animationFrameId: number;
    let timeoutId: NodeJS.Timeout;

    const initializeViewer = () => {
      if (!viewerRef.current) return;

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
              console.error('Failed to resize viewer after initialization:', error);
            }
          }
        }, 100);
        
      } catch (error) {
        console.error('Failed to initialize 3Dmol viewer:', error);
        
        // 再試行ロジック（最大3回まで）
        if (retryCount < 3) {
          setRetryCount(prev => prev + 1);
          timeoutId = setTimeout(() => {
            animationFrameId = requestAnimationFrame(initializeViewer);
          }, 500 * (retryCount + 1)); // 指数バックオフ
        } else {
          onError?.('Failed to initialize molecular viewer after multiple attempts');
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
  }, [onError, retryCount]);

  // ResizeObserverのセットアップ
  useEffect(() => {
    if (!viewerRef.current) return;

    // ResizeObserver を設定してDOM要素のサイズ変更を監視
    resizeObserverRef.current = new ResizeObserver((entries) => {
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
  }, [orbitalsData, selectedOrbitalIndex]);

  // CUBEデータが更新されたときに分子軌道を表示
  useEffect(() => {
    if (!viewer || !cubeData || !(cubeData as any).cube_data) return;

    setIsLoading(true);
    
    try {
      // ビューアーをクリア
      viewer.clear();

      // CUBEファイルデータを追加
      const cubeContent = (cubeData as any).cube_data;
      
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

      // 分子構造も表示（CUBEファイルから分子情報を抽出）
      viewer.addModel(cubeContent, 'cube');
      viewer.setStyle({}, { stick: { radius: 0.1 }, sphere: { radius: 0.3 } });

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
          }
        } catch (delayedRenderError) {
          console.error('Failed to perform delayed re-render:', delayedRenderError);
        }
      }, 200);

      setIsLoading(false);
    } catch (error) {
      console.error('Failed to render molecular orbital:', error);
      onError?.('Failed to render molecular orbital visualization');
      setIsLoading(false);
    }
  }, [viewer, cubeData, viewerOptions, onError]);

  // エラーハンドリング
  useEffect(() => {
    if (orbitalsError) {
      onError?.(orbitalsError.message || 'Failed to load orbital information');
    }
  }, [orbitalsError, onError]);

  useEffect(() => {
    if (cubeError) {
      onError?.(cubeError.message || 'Failed to load orbital visualization data');
    }
  }, [cubeError, onError]);

  const handleOrbitalChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    const newIndex = parseInt(event.target.value, 10);
    setSelectedOrbitalIndex(newIndex);
  };

  const handleOptionsChange = (newOptions: Partial<ViewerOptions>) => {
    setViewerOptions(prev => ({ ...prev, ...newOptions }));
  };

  if (orbitalsLoading) {
    return (
      <div style={{ textAlign: 'center', padding: '40px' }}>
        <div>⚛️ Loading orbital information...</div>
      </div>
    );
  }

  if (!orbitalsData || !orbitalsData.orbitals || orbitalsData.orbitals.length === 0) {
    return (
      <div style={{ textAlign: 'center', padding: '40px', color: '#666' }}>
        <div>📊 No orbital data available for this calculation.</div>
      </div>
    );
  }

  const selectedOrbital = orbitalsData.orbitals.find(
    (orbital: OrbitalInfo) => orbital.index === selectedOrbitalIndex
  );

  return (
    <div style={{ border: '1px solid #ddd', borderRadius: '8px', overflow: 'hidden' }}>
      {/* コントロールパネル */}
      <div
        style={{
          padding: '16px',
          backgroundColor: '#f8f9fa',
          borderBottom: '1px solid #ddd',
        }}
      >
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
            gap: '16px',
            alignItems: 'center',
          }}
        >
          {/* 軌道選択 */}
          <div>
            <label htmlFor="orbital-select" style={{ display: 'block', marginBottom: '4px', fontWeight: 'bold' }}>
              分子軌道:
            </label>
            <select
              id="orbital-select"
              value={selectedOrbitalIndex || ''}
              onChange={handleOrbitalChange}
              style={{
                width: '100%',
                padding: '6px',
                border: '1px solid #ccc',
                borderRadius: '4px',
              }}
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
            <label htmlFor="grid-size" style={{ display: 'block', marginBottom: '4px', fontWeight: 'bold' }}>
              グリッドサイズ: {viewerOptions.gridSize}
            </label>
            <input
              id="grid-size"
              type="range"
              min="40"
              max="120"
              step="10"
              value={viewerOptions.gridSize}
              onChange={(e) =>
                handleOptionsChange({ gridSize: parseInt(e.target.value, 10) })
              }
              style={{ width: '100%' }}
            />
          </div>

          {/* 正の等値面 */}
          <div>
            <label htmlFor="isovalue-pos" style={{ display: 'block', marginBottom: '4px', fontWeight: 'bold' }}>
              正等値面: {viewerOptions.isovaluePos.toFixed(3)}
            </label>
            <input
              id="isovalue-pos"
              type="range"
              min="0.001"
              max="0.1"
              step="0.001"
              value={viewerOptions.isovaluePos}
              onChange={(e) =>
                handleOptionsChange({ isovaluePos: parseFloat(e.target.value) })
              }
              style={{ width: '100%' }}
            />
          </div>

          {/* 負の等値面 */}
          <div>
            <label htmlFor="isovalue-neg" style={{ display: 'block', marginBottom: '4px', fontWeight: 'bold' }}>
              負等値面: {viewerOptions.isovalueNeg.toFixed(3)}
            </label>
            <input
              id="isovalue-neg"
              type="range"
              min="-0.1"
              max="-0.001"
              step="0.001"
              value={viewerOptions.isovalueNeg}
              onChange={(e) =>
                handleOptionsChange({ isovalueNeg: parseFloat(e.target.value) })
              }
              style={{ width: '100%' }}
            />
          </div>
        </div>

        {/* 選択された軌道の情報 */}
        {selectedOrbital && (
          <div
            style={{
              marginTop: '12px',
              padding: '12px',
              backgroundColor: 'white',
              borderRadius: '4px',
              border: '1px solid #e0e0e0',
            }}
          >
            <div style={{ fontSize: '14px', color: '#333' }}>
              <strong>{selectedOrbital.label}</strong> (軌道 #{selectedOrbital.index})
              <br />
              エネルギー: <code>{selectedOrbital.energy_ev.toFixed(4)} eV</code> (
              <code>{selectedOrbital.energy_hartree.toFixed(6)} a.u.</code>)
              <br />
              占有: <code>{selectedOrbital.occupancy}</code>
            </div>
          </div>
        )}
      </div>

      {/* 3Dmol.jsビューアー */}
      <div style={{ position: 'relative' }}>
        <div
          ref={viewerRef}
          style={{
            width: '100%',
            height: '500px',
            backgroundColor: 'white',
          }}
        />
        
        {/* ローディングオーバーレイ */}
        {(isLoading || cubeLoading) && (
          <div
            style={{
              position: 'absolute',
              top: 0,
              left: 0,
              right: 0,
              bottom: 0,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              backgroundColor: 'rgba(255, 255, 255, 0.8)',
              fontSize: '16px',
            }}
          >
            ⚛️ 軌道データを生成中...
          </div>
        )}
      </div>
    </div>
  );
};