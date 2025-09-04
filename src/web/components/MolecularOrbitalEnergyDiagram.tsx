// src/web/components/MolecularOrbitalEnergyDiagram.tsx

import React, { useMemo } from 'react';
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

export const MolecularOrbitalEnergyDiagram: React.FC<
  MolecularOrbitalEnergyDiagramProps
> = ({ calculationId, onError, onOrbitalSelect, selectedOrbitalIndex }) => {
  // 軌道情報を取得
  const {
    data: orbitalsData,
    isLoading: orbitalsLoading,
    error: orbitalsError,
  } = useGetOrbitals(calculationId);

  // エラーハンドリング
  React.useEffect(() => {
    if (orbitalsError) {
      console.error('Failed to load orbital information:', orbitalsError);
      onError?.(orbitalsError.message || 'Failed to load orbital information');
    }
  }, [orbitalsError, onError]);

  // 軌道データの処理とソート
  const processedOrbitals: ProcessedOrbital[] = useMemo(() => {
    if (!orbitalsData?.orbitals) return [];

    // エネルギー順にソート
    const sortedOrbitals = [...orbitalsData.orbitals].sort(
      (a, b) => a.energy_ev - b.energy_ev
    );

    // Y軸位置を計算
    const energyRange = Math.max(
      sortedOrbitals[sortedOrbitals.length - 1].energy_ev -
        sortedOrbitals[0].energy_ev,
      10 // 最小範囲を設定
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
        drawableHeight * (1 - (orbital.energy_ev - minEnergy) / energyRange),
      displayLevel: index,
    }));
  }, [orbitalsData]);

  // HOMO-LUMO情報の計算
  const orbitalSummary = useMemo(() => {
    const homoOrbital = processedOrbitals.find(o => o.orbital_type === 'homo');
    const lumoOrbital = processedOrbitals.find(o => o.orbital_type === 'lumo');

    return {
      homoOrbital,
      lumoOrbital,
      homoLumoGap:
        homoOrbital && lumoOrbital
          ? lumoOrbital.energy_ev - homoOrbital.energy_ev
          : null,
    };
  }, [processedOrbitals]);

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
          Calculation is not complete or orbital data has not been generated.
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
        <h3 className={styles.diagramTitle}>Molecular Orbital Energy Level Diagram</h3>
        {orbitalSummary.homoLumoGap && (
          <div className={styles.homoLumoGap}>
            <strong>HOMO-LUMO Gap:</strong>{' '}
            <span className={styles.gapValue}>
              {orbitalSummary.homoLumoGap.toFixed(4)} eV
            </span>
          </div>
        )}
      </div>

      {/* SVGエネルギー準位図 */}
      <div className={styles.svgContainer}>
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
            strokeWidth="2"
          />

          {/* Y軸ラベル */}
          <text
            x={20}
            y={DIAGRAM_CONFIG.height / 2}
            textAnchor="middle"
            fontSize="14"
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
          {processedOrbitals.map(orbital => {
            const isSelected = selectedOrbitalIndex === orbital.index;
            const color = getOrbitalColor(orbital);
            const x =
              DIAGRAM_CONFIG.margin.left +
              chartWidth / 2 -
              DIAGRAM_CONFIG.orbitalWidth / 2;

            return (
              <g key={orbital.index}>
                {/* 軌道線 */}
                <rect
                  x={x}
                  y={orbital.yPosition - DIAGRAM_CONFIG.orbitalHeight / 2}
                  width={DIAGRAM_CONFIG.orbitalWidth}
                  height={DIAGRAM_CONFIG.orbitalHeight}
                  fill={color}
                  stroke={isSelected ? '#f39c12' : color}
                  strokeWidth={isSelected ? 3 : 1}
                  style={{
                    cursor: 'pointer',
                    transition: 'all 0.2s ease',
                    filter: isSelected ? 'brightness(1.1)' : 'none',
                  }}
                  onClick={() => handleOrbitalClick(orbital)}
                  onMouseEnter={e => {
                    e.currentTarget.style.filter = 'brightness(1.2)';
                    e.currentTarget.style.strokeWidth = '2';
                  }}
                  onMouseLeave={e => {
                    e.currentTarget.style.filter = isSelected
                      ? 'brightness(1.1)'
                      : 'none';
                    e.currentTarget.style.strokeWidth = isSelected ? '3' : '1';
                  }}
                >
                  <title>
                    {orbital.label || `Orbital ${orbital.index}`}:{' '}
                    {orbital.energy_ev.toFixed(4)} eV
                  </title>
                </rect>

                {/* 軌道ラベル（主要な軌道のみ表示） */}
                {(orbital.orbital_type === 'homo' ||
                  orbital.orbital_type === 'lumo' ||
                  orbital.label?.includes('HOMO') ||
                  orbital.label?.includes('LUMO')) && (
                  <>
                    {/* ラベル */}
                    <text
                      x={x + DIAGRAM_CONFIG.orbitalWidth + 10}
                      y={orbital.yPosition + 4}
                      fontSize="12"
                      fill="#333"
                      style={{ fontWeight: 'bold' }}
                    >
                      {orbital.label || `Orbital ${orbital.index}`}
                    </text>
                    {/* エネルギー値 */}
                    <text
                      x={x + DIAGRAM_CONFIG.orbitalWidth + 10}
                      y={orbital.yPosition + 18}
                      fontSize="10"
                      fill="#666"
                      fontFamily="monospace"
                    >
                      {orbital.energy_ev.toFixed(3)} eV
                    </text>
                  </>
                )}

                {/* 占有を示す電子（占有軌道の場合） */}
                {orbital.occupancy > 0 && (
                  <circle
                    cx={x + DIAGRAM_CONFIG.orbitalWidth / 4}
                    cy={orbital.yPosition}
                    r="3"
                    fill="#34495e"
                  />
                )}
                {orbital.occupancy > 1 && (
                  <circle
                    cx={x + (3 * DIAGRAM_CONFIG.orbitalWidth) / 4}
                    cy={orbital.yPosition}
                    r="3"
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
                  Math.max(1, Math.floor(processedOrbitals.length / 10)) ===
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
                  strokeWidth="1"
                />
                <text
                  x={DIAGRAM_CONFIG.margin.left - 10}
                  y={orbital.yPosition + 4}
                  textAnchor="end"
                  fontSize="10"
                  fill="#666"
                  fontFamily="monospace"
                >
                  {orbital.energy_ev.toFixed(1)}
                </text>
              </g>
            ))}
        </svg>
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
          <div className={`${styles.legendColorBox} ${styles.legendHomo}`} />
          <span>HOMO</span>
        </div>
        <div className={styles.legendItem}>
          <div className={`${styles.legendColorBox} ${styles.legendLumo}`} />
          <span>LUMO</span>
        </div>
        <div className={styles.legendItem}>
          <div
            className={`${styles.legendColorBox} ${styles.legendOccupied}`}
          />
          <span>Occupied Orbitals</span>
        </div>
        <div className={styles.legendItem}>
          <div className={`${styles.legendColorBox} ${styles.legendVirtual}`} />
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
                  {selectedOrbital.label || `Orbital ${selectedOrbital.index}`}
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
                    <strong>Orbital type:</strong> {selectedOrbital.orbital_type}
                  </div>
                </div>
              </div>
            );
          })()}
        </div>
      )}

      {/* 操作説明 */}
      <div className={styles.helpSection}>
        💡
        Click on orbitals to view details. The yellow gap area indicates the HOMO-LUMO gap.
      </div>
    </div>
  );
};
