// src/web/components/MolecularOrbitalEnergyDiagram.tsx

import React, { useMemo, useState, useRef, useEffect } from 'react';
import {
  ReactSVGPanZoom,
  TOOL_AUTO,
  fitToViewer,
  Tool,
  Value,
} from 'react-svg-pan-zoom';
import { useGetOrbitals } from '../hooks/useCalculationQueries';
import { OrbitalInfo } from '../types/api-types';
import styles from './MolecularOrbitalEnergyDiagram.module.css';

interface MolecularOrbitalEnergyDiagramProps {
  calculationId: string;
  onError?: (error: string) => void;
  onOrbitalSelect?: (orbitalIndex: number) => void;
  selectedOrbitalIndex?: number | null;
}

interface ProcessedOrbital extends OrbitalInfo {
  yPosition: number;
  displayLevel: number;
}

const DIAGRAM_CONFIG = {
  width: 800,
  height: 600,
  margin: { top: 40, right: 100, bottom: 60, left: 80 },
  orbitalWidth: 120,
  orbitalHeight: 3,
  gapThreshold: 0.5, // eV - threshold for highlighting HOMO-LUMO gap
};

export const MolecularOrbitalEnergyDiagram: React.FC<MolecularOrbitalEnergyDiagramProps> =
  React.memo(
    ({ calculationId, onError, onOrbitalSelect, selectedOrbitalIndex }) => {
      // Viewer state for zoom/pan
      const viewerRef = useRef<any>(null);
      const [tool, setTool] = useState<Tool>(TOOL_AUTO);
      const [value, setValue] = useState<Value | null>(null);
      const [viewerSize, setViewerSize] = useState({ width: 800, height: 600 });
      const containerRef = useRef<HTMLDivElement>(null);

      // Get orbital information
      const {
        data: orbitalsData,
        isLoading: orbitalsLoading,
        error: orbitalsError,
      } = useGetOrbitals(calculationId);

      // レスポンシブ対応: コンテナサイズの監視とViewerサイズの更新（ResizeObserver使用）
      useEffect(() => {
        if (!containerRef.current) return;

        const updateSize = () => {
          if (containerRef.current) {
            const { width } = containerRef.current.getBoundingClientRect();
            const newWidth = Math.max(width, 600);
            setViewerSize({
              width: newWidth,
              height: 600,
            });

            // Set default value if initial value is not set
            if (value === null) {
              setValue({
                version: 2,
                mode: 'idle',
                focus: false,
                a: 1,
                b: 0,
                c: 0,
                d: 1,
                e: 0,
                f: 0,
                viewerWidth: newWidth,
                viewerHeight: 600,
                SVGWidth: DIAGRAM_CONFIG.width,
                SVGHeight: DIAGRAM_CONFIG.height,
                startX: null,
                startY: null,
                endX: null,
                endY: null,
                miniatureOpen: false,
                focus_miniature: false,
              } as Value);
            }
          }
        };

        // Initial size setup
        updateSize();

        // Monitor container size changes with ResizeObserver
        const resizeObserver = new ResizeObserver(() => {
          updateSize();
        });

        resizeObserver.observe(containerRef.current);

        return () => {
          resizeObserver.disconnect();
        };
      }, [value]);

      // エラーハンドリング
      React.useEffect(() => {
        if (orbitalsError) {
          console.error('Failed to load orbital information:', orbitalsError);
          onError?.(
            orbitalsError.message || 'Failed to load orbital information'
          );
        }
      }, [orbitalsError, onError]);

      // Process and sort orbital data
      const processedOrbitals: ProcessedOrbital[] = useMemo(() => {
        if (!orbitalsData?.orbitals) return [];

        // Sort by energy
        const sortedOrbitals = [...orbitalsData.orbitals].sort(
          (a, b) => a.energy_ev - b.energy_ev
        );

        // Calculate Y-axis position
        const energyRange = Math.max(
          sortedOrbitals[sortedOrbitals.length - 1].energy_ev -
            sortedOrbitals[0].energy_ev,
          10 // Set minimum range
        );
        const minEnergy = sortedOrbitals[0].energy_ev;
        const drawableHeight =
          DIAGRAM_CONFIG.height -
          DIAGRAM_CONFIG.margin.top -
          DIAGRAM_CONFIG.margin.bottom;

        return sortedOrbitals.map((orbital, index) => ({
          ...orbital,
          yPosition:
            DIAGRAM_CONFIG.margin.top +
            drawableHeight *
              (1 - (orbital.energy_ev - minEnergy) / energyRange),
          displayLevel: index,
        }));
      }, [orbitalsData]);

      // Calculate HOMO-LUMO information
      const orbitalSummary = useMemo(() => {
        const homoOrbital = processedOrbitals.find(
          o => o.orbital_type === 'homo'
        );
        const lumoOrbital = processedOrbitals.find(
          o => o.orbital_type === 'lumo'
        );

        return {
          homoOrbital,
          lumoOrbital,
          homoLumoGap:
            homoOrbital && lumoOrbital
              ? lumoOrbital.energy_ev - homoOrbital.energy_ev
              : null,
        };
      }, [processedOrbitals]);

      // Calculate zoom scale and smart scaling
      const scalingFactors = useMemo(() => {
        const scale = value?.a || 1;

        // ベースサイズ
        const baseFontSize = 12;
        const baseSmallFontSize = 10;
        const baseLabelFontSize = 14;
        const baseStrokeWidth = 2;
        const baseThinStrokeWidth = 1;

        // スケールに応じて逆スケーリング（拡大時にテキストが相対的に小さくなる）
        const invScale = 1 / scale;

        // ズームレベルに応じた表示密度の決定
        let labelDensity: 'sparse' | 'medium' | 'dense' | 'all';
        if (scale < 1.5) {
          labelDensity = 'sparse'; // 主要な軌道のみ
        } else if (scale < 3) {
          labelDensity = 'medium'; // 一部の軌道
        } else if (scale < 5) {
          labelDensity = 'dense'; // 多くの軌道
        } else {
          labelDensity = 'all'; // 全ての軌道
        }

        return {
          scale,
          invScale,
          fontSize: baseFontSize * invScale,
          smallFontSize: baseSmallFontSize * invScale,
          labelFontSize: baseLabelFontSize * invScale,
          strokeWidth: baseStrokeWidth * invScale,
          thinStrokeWidth: baseThinStrokeWidth * invScale,
          labelDensity,
        };
      }, [value]);

      // 初期ビューの設定とデータ変更時のリセット
      useEffect(() => {
        if (
          viewerRef.current &&
          processedOrbitals.length > 0 &&
          viewerSize.width > 0
        ) {
          // 二重のrequestAnimationFrameでDOMの準備完了を確実に待つ
          requestAnimationFrame(() => {
            requestAnimationFrame(() => {
              if (viewerRef.current) {
                const viewer = viewerRef.current;
                const initialValue = viewer.getValue();
                const fittedValue = fitToViewer(initialValue);
                setValue(fittedValue);
              }
            });
          });
        }
      }, [processedOrbitals.length, viewerSize.width]);

      // 軌道の色を決定
      const getOrbitalColor = (orbital: OrbitalInfo): string => {
        switch (orbital.orbital_type) {
          case 'homo':
            return '#e74c3c'; // 赤
          case 'lumo':
            return '#3498db'; // 青
          case 'core':
            return '#2ecc71'; // 緑
          case 'virtual':
            return '#95a5a6'; // 灰色
          default:
            if (orbital.occupancy > 0) {
              return '#2ecc71'; // 占有軌道は緑
            }
            return '#95a5a6'; // 仮想軌道は灰色
        }
      };

      // 軌道クリックハンドラー
      const handleOrbitalClick = (orbital: ProcessedOrbital) => {
        onOrbitalSelect?.(orbital.index);
      };

      // ズームレベルに応じて軌道ラベルを表示すべきか判定
      const shouldShowLabel = (
        orbital: ProcessedOrbital,
        index: number
      ): boolean => {
        const { labelDensity } = scalingFactors;

        // 常に表示する軌道（HOMO, LUMO）
        if (
          orbital.orbital_type === 'homo' ||
          orbital.orbital_type === 'lumo' ||
          orbital.label?.includes('HOMO') ||
          orbital.label?.includes('LUMO')
        ) {
          return true;
        }

        // ズームレベルに応じた表示密度
        switch (labelDensity) {
          case 'all':
            return true; // 全ての軌道を表示
          case 'dense':
            return index % 2 === 0; // 2つおきに表示
          case 'medium':
            return index % 4 === 0; // 4つおきに表示
          case 'sparse':
            return false; // 主要な軌道のみ（上記のHOMO/LUMOチェックで既に処理済み）
          default:
            return false;
        }
      };

      if (orbitalsLoading) {
        return (
          <div className={styles.loadingContainer}>
            <div className={styles.loadingText}>
              ⚛️ Loading energy level data...
            </div>
          </div>
        );
      }

      if (orbitalsError) {
        return (
          <div className={styles.errorContainer}>
            <div>❌ Failed to load energy level data</div>
            <div className={styles.errorMessage}>
              {orbitalsError.message || 'An unknown error occurred'}
            </div>
          </div>
        );
      }

      if (processedOrbitals.length === 0) {
        return (
          <div className={styles.noDataContainer}>
            <div>📊 No orbital energy information available.</div>
            <div className={styles.noDataMessage}>
              Calculation is not complete or orbital data has not been
              generated.
            </div>
          </div>
        );
      }

      const chartWidth =
        DIAGRAM_CONFIG.width -
        DIAGRAM_CONFIG.margin.left -
        DIAGRAM_CONFIG.margin.right;

      return (
        <div className={styles.diagramContainer}>
          {/* ヘッダー情報 */}
          <div className={styles.diagramHeader}>
            <h3 className={styles.diagramTitle}>
              Molecular Orbital Energy Level Diagram
            </h3>
            {orbitalSummary.homoLumoGap && (
              <div className={styles.homoLumoGap}>
                <strong>HOMO-LUMO Gap:</strong>{' '}
                <span className={styles.gapValue}>
                  {orbitalSummary.homoLumoGap.toFixed(4)} eV
                </span>
              </div>
            )}
          </div>

          {/* SVGエネルギー準位図 with Zoom/Pan */}
          <div ref={containerRef} className={styles.viewerContainer}>
            <ReactSVGPanZoom
              ref={viewerRef}
              width={viewerSize.width}
              height={viewerSize.height}
              tool={tool}
              onChangeTool={setTool}
              value={value}
              onChangeValue={setValue}
              detectAutoPan={false}
              preventPanOutside={true}
              background="#fafafa"
              SVGBackground="#fafafa"
              toolbarProps={{
                position: 'right',
                SVGAlignX: 'center',
                SVGAlignY: 'top',
              }}
              miniatureProps={{
                position: 'none',
                background: '#fafafa',
                width: 100,
                height: 80,
              }}
              scaleFactorMin={0.8}
              scaleFactorMax={10}
              scaleFactorOnWheel={1.06}
            >
              <svg
                width={DIAGRAM_CONFIG.width}
                height={DIAGRAM_CONFIG.height}
                className={styles.diagramSvg}
              >
                {/* 背景グリッド */}
                <defs>
                  <pattern
                    id="grid"
                    width="40"
                    height="40"
                    patternUnits="userSpaceOnUse"
                  >
                    <path
                      d="M 40 0 L 0 0 0 40"
                      fill="none"
                      stroke="#f5f5f5"
                      strokeWidth="1"
                    />
                  </pattern>
                </defs>
                <rect
                  width={DIAGRAM_CONFIG.width}
                  height={DIAGRAM_CONFIG.height}
                  fill="url(#grid)"
                />

                {/* Y軸 */}
                <line
                  x1={DIAGRAM_CONFIG.margin.left}
                  y1={DIAGRAM_CONFIG.margin.top}
                  x2={DIAGRAM_CONFIG.margin.left}
                  y2={DIAGRAM_CONFIG.height - DIAGRAM_CONFIG.margin.bottom}
                  stroke="#333"
                  strokeWidth={scalingFactors.strokeWidth}
                />

                {/* Y軸ラベル */}
                <text
                  x={20}
                  y={DIAGRAM_CONFIG.height / 2}
                  textAnchor="middle"
                  fontSize={scalingFactors.labelFontSize}
                  fill="#666"
                  transform={`rotate(-90, 20, ${DIAGRAM_CONFIG.height / 2})`}
                >
                  Energy (eV)
                </text>

                {/* HOMO-LUMOギャップの強調表示 */}
                {orbitalSummary.homoOrbital && orbitalSummary.lumoOrbital && (
                  <rect
                    x={DIAGRAM_CONFIG.margin.left}
                    y={orbitalSummary.lumoOrbital.yPosition}
                    width={chartWidth}
                    height={
                      orbitalSummary.homoOrbital.yPosition -
                      orbitalSummary.lumoOrbital.yPosition
                    }
                    fill="rgba(255, 193, 7, 0.1)"
                    stroke="rgba(255, 193, 7, 0.3)"
                    strokeWidth="1"
                    strokeDasharray="5,5"
                  />
                )}

                {/* 軌道レベル */}
                {processedOrbitals.map((orbital, index) => {
                  const isSelected = selectedOrbitalIndex === orbital.index;
                  const color = getOrbitalColor(orbital);
                  const x =
                    DIAGRAM_CONFIG.margin.left +
                    chartWidth / 2 -
                    DIAGRAM_CONFIG.orbitalWidth / 2;

                  // ズームスケールに応じた軌道の高さ
                  const scaledOrbitalHeight =
                    DIAGRAM_CONFIG.orbitalHeight * scalingFactors.invScale;

                  // 電子の半径もスケール調整
                  const electronRadius = 3 * scalingFactors.invScale;

                  // ラベルを表示するか判定
                  const showLabel = shouldShowLabel(orbital, index);

                  return (
                    <g key={orbital.index}>
                      {/* 軌道線 */}
                      <rect
                        x={x}
                        y={orbital.yPosition - scaledOrbitalHeight / 2}
                        width={DIAGRAM_CONFIG.orbitalWidth}
                        height={scaledOrbitalHeight}
                        fill={color}
                        stroke={isSelected ? '#f39c12' : color}
                        strokeWidth={
                          isSelected
                            ? scalingFactors.strokeWidth * 1.5
                            : scalingFactors.thinStrokeWidth
                        }
                        style={{
                          cursor: 'pointer',
                          transition: 'all 0.2s ease',
                          filter: isSelected ? 'brightness(1.1)' : 'none',
                        }}
                        onClick={() => handleOrbitalClick(orbital)}
                        onMouseEnter={e => {
                          e.currentTarget.style.filter = 'brightness(1.2)';
                        }}
                        onMouseLeave={e => {
                          e.currentTarget.style.filter = isSelected
                            ? 'brightness(1.1)'
                            : 'none';
                        }}
                      >
                        <title>
                          #{orbital.index}: {orbital.energy_ev.toFixed(4)} eV
                          {orbital.label ? ` (${orbital.label})` : ''}
                        </title>
                      </rect>

                      {/* 軌道ラベル：番号とエネルギーを1行で表示 */}
                      {showLabel && (
                        <text
                          x={x + DIAGRAM_CONFIG.orbitalWidth + 10}
                          y={orbital.yPosition + scalingFactors.fontSize / 3}
                          fontSize={scalingFactors.fontSize}
                          fill={
                            orbital.orbital_type === 'homo' ||
                            orbital.orbital_type === 'lumo'
                              ? '#333'
                              : '#666'
                          }
                          fontFamily="monospace"
                          style={{
                            fontWeight:
                              orbital.orbital_type === 'homo' ||
                              orbital.orbital_type === 'lumo'
                                ? 'bold'
                                : 'normal',
                          }}
                        >
                          #{orbital.index}: {orbital.energy_ev.toFixed(3)} eV
                          {orbital.label &&
                            (orbital.label.includes('HOMO') ||
                              orbital.label.includes('LUMO')) &&
                            ` (${orbital.label})`}
                        </text>
                      )}

                      {/* 占有を示す電子（占有軌道の場合） */}
                      {orbital.occupancy > 0 && (
                        <circle
                          cx={x + DIAGRAM_CONFIG.orbitalWidth / 4}
                          cy={orbital.yPosition}
                          r={electronRadius}
                          fill="#34495e"
                        />
                      )}
                      {orbital.occupancy > 1 && (
                        <circle
                          cx={x + (3 * DIAGRAM_CONFIG.orbitalWidth) / 4}
                          cy={orbital.yPosition}
                          r={electronRadius}
                          fill="#34495e"
                        />
                      )}
                    </g>
                  );
                })}

                {/* エネルギー軸の目盛り */}
                {processedOrbitals
                  .filter(
                    (_, index) =>
                      index %
                        Math.max(
                          1,
                          Math.floor(processedOrbitals.length / 10)
                        ) ===
                      0
                  )
                  .map(orbital => (
                    <g key={`tick-${orbital.index}`}>
                      <line
                        x1={DIAGRAM_CONFIG.margin.left - 5}
                        y1={orbital.yPosition}
                        x2={DIAGRAM_CONFIG.margin.left}
                        y2={orbital.yPosition}
                        stroke="#666"
                        strokeWidth={scalingFactors.thinStrokeWidth}
                      />
                      <text
                        x={DIAGRAM_CONFIG.margin.left - 10}
                        y={orbital.yPosition + scalingFactors.smallFontSize / 3}
                        textAnchor="end"
                        fontSize={scalingFactors.smallFontSize}
                        fill="#666"
                        fontFamily="monospace"
                      >
                        {orbital.energy_ev.toFixed(1)}
                      </text>
                    </g>
                  ))}
              </svg>
            </ReactSVGPanZoom>
          </div>

          {/* 凡例 */}
          <div
            style={{
              marginTop: '20px',
              display: 'flex',
              gap: '20px',
              fontSize: '14px',
            }}
          >
            <div className={styles.legendItem}>
              <div
                className={`${styles.legendColorBox} ${styles.legendHomo}`}
              />
              <span>HOMO</span>
            </div>
            <div className={styles.legendItem}>
              <div
                className={`${styles.legendColorBox} ${styles.legendLumo}`}
              />
              <span>LUMO</span>
            </div>
            <div className={styles.legendItem}>
              <div
                className={`${styles.legendColorBox} ${styles.legendOccupied}`}
              />
              <span>Occupied Orbitals</span>
            </div>
            <div className={styles.legendItem}>
              <div
                className={`${styles.legendColorBox} ${styles.legendVirtual}`}
              />
              <span>Virtual Orbitals</span>
            </div>
            <div className={styles.legendItem}>
              <div className={styles.legendElectron} />
              <span>Electrons</span>
            </div>
          </div>

          {/* 選択された軌道の詳細情報 */}
          {selectedOrbitalIndex !== null && (
            <div className={styles.selectedOrbitalDetails}>
              {(() => {
                const selectedOrbital = processedOrbitals.find(
                  o => o.index === selectedOrbitalIndex
                );
                if (!selectedOrbital) return null;

                return (
                  <div>
                    <div className={styles.selectedOrbitalTitle}>
                      Selected Orbital:{' '}
                      {selectedOrbital.label ||
                        `Orbital ${selectedOrbital.index}`}
                    </div>
                    <div className={styles.selectedOrbitalGrid}>
                      <div>
                        <strong>Energy:</strong>{' '}
                        {selectedOrbital.energy_ev.toFixed(4)} eV
                      </div>
                      <div>
                        <strong>Energy (a.u.):</strong>{' '}
                        {selectedOrbital.energy_hartree.toFixed(6)}
                      </div>
                      <div>
                        <strong>Occupancy:</strong> {selectedOrbital.occupancy}
                      </div>
                      <div>
                        <strong>Orbital type:</strong>{' '}
                        {selectedOrbital.orbital_type}
                      </div>
                    </div>
                  </div>
                );
              })()}
            </div>
          )}
        </div>
      );
    }
  );
