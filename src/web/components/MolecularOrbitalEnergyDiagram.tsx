// src/web/components/MolecularOrbitalEnergyDiagram.tsx
// D3.js-based implementation for better zoom/pan and label handling

import React, { useMemo, useState, useRef, useEffect, useCallback } from 'react';
import * as d3 from 'd3';
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

interface LabelInfo {
  orbital: ProcessedOrbital;
  x: number;
  y: number;
  visible: boolean;
  offsetY: number;
}

const DIAGRAM_CONFIG = {
  width: 800,
  height: 600,
  margin: { top: 40, right: 150, bottom: 60, left: 160 },
  orbitalWidth: 80, // Half width for better appearance
  orbitalHeight: 3,
  gapThreshold: 0.5,
  minLabelSpacing: 14,
};

export const MolecularOrbitalEnergyDiagram: React.FC<MolecularOrbitalEnergyDiagramProps> =
  React.memo(
    ({ calculationId, onError, onOrbitalSelect, selectedOrbitalIndex }) => {
      const svgRef = useRef<SVGSVGElement>(null);
      const containerRef = useRef<HTMLDivElement>(null);
      const [viewerSize, setViewerSize] = useState({ width: 800, height: 600 });

      // Get orbital information
      const {
        data: orbitalsData,
        isLoading: orbitalsLoading,
        error: orbitalsError,
      } = useGetOrbitals(calculationId);

      // Responsive sizing
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
          }
        };

        updateSize();
        const resizeObserver = new ResizeObserver(updateSize);
        resizeObserver.observe(containerRef.current);

        return () => resizeObserver.disconnect();
      }, []);

      // Error handling
      useEffect(() => {
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

        const sortedOrbitals = [...orbitalsData.orbitals].sort(
          (a, b) => a.energy_ev - b.energy_ev
        );

        const energyRange = Math.max(
          sortedOrbitals[sortedOrbitals.length - 1].energy_ev -
            sortedOrbitals[0].energy_ev,
          10
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

      // Get orbital color
      const getOrbitalColor = useCallback((orbital: OrbitalInfo): string => {
        // Check for SOMO first (occupancy = 1 or label contains SOMO)
        const isSomo = orbital.label?.includes('SOMO') || 
          (orbital.occupancy > 0.5 && orbital.occupancy < 1.5);
        
        if (isSomo) {
          return '#f39c12'; // Orange for SOMO
        }
        
        switch (orbital.orbital_type) {
          case 'homo':
            return '#e74c3c';
          case 'lumo':
            return '#3498db';
          case 'core':
            return '#2ecc71';
          case 'virtual':
            return '#95a5a6';
          default:
            return orbital.occupancy > 0 ? '#2ecc71' : '#95a5a6';
        }
      }, []);

      // Calculate label positions with collision detection
      const calculateLabelPositions = useCallback(
        (orbitals: ProcessedOrbital[], scale: number): LabelInfo[] => {
          const invScale = 1 / scale;
          // Label X position: right after orbital end (scales with zoom)
          const labelX =
            DIAGRAM_CONFIG.margin.left +
            10 * invScale +
            DIAGRAM_CONFIG.orbitalWidth * invScale +
            10 * invScale;

          // Determine which labels to show based on zoom level
          const labelDensity =
            scale < 1.5
              ? 'sparse'
              : scale < 3
                ? 'medium'
                : scale < 5
                  ? 'dense'
                  : 'all';

          const labels: LabelInfo[] = orbitals.map((orbital, index) => {
            const isImportant =
              orbital.orbital_type === 'homo' ||
              orbital.orbital_type === 'lumo' ||
              orbital.label?.includes('HOMO') ||
              orbital.label?.includes('LUMO') ||
              orbital.label?.includes('SOMO');

            let visible = false;
            if (isImportant) {
              visible = true;
            } else {
              switch (labelDensity) {
                case 'all':
                  visible = true;
                  break;
                case 'dense':
                  visible = index % 2 === 0;
                  break;
                case 'medium':
                  visible = index % 4 === 0;
                  break;
                case 'sparse':
                  visible = false;
                  break;
              }
            }

            return {
              orbital,
              x: labelX,
              y: orbital.yPosition,
              visible,
              offsetY: 0,
            };
          });

          // Collision detection and resolution for visible labels
          const visibleLabels = labels.filter(l => l.visible);
          const minSpacing = DIAGRAM_CONFIG.minLabelSpacing / scale;

          // Sort by y position
          visibleLabels.sort((a, b) => a.y - b.y);

          // Resolve collisions
          for (let i = 1; i < visibleLabels.length; i++) {
            const prev = visibleLabels[i - 1];
            const curr = visibleLabels[i];
            const actualPrevY = prev.y + prev.offsetY;
            const actualCurrY = curr.y + curr.offsetY;
            const spacing = actualCurrY - actualPrevY;

            if (spacing < minSpacing) {
              curr.offsetY += minSpacing - spacing;
            }
          }

          return labels;
        },
        []
      );

      // D3 zoom behavior and rendering
      useEffect(() => {
        if (!svgRef.current || processedOrbitals.length === 0) return;

        const svg = d3.select(svgRef.current);
        const chartWidth =
          DIAGRAM_CONFIG.width -
          DIAGRAM_CONFIG.margin.left -
          DIAGRAM_CONFIG.margin.right;

        // Clear previous content
        svg.selectAll('*').remove();

        // Create main group for zooming
        const mainGroup = svg.append('g').attr('class', 'main-group');

        // Background
        mainGroup
          .append('rect')
          .attr('width', DIAGRAM_CONFIG.width)
          .attr('height', DIAGRAM_CONFIG.height)
          .attr('fill', '#fafafa');

        // Grid pattern
        const defs = svg.append('defs');
        const pattern = defs
          .append('pattern')
          .attr('id', 'grid-pattern')
          .attr('width', 40)
          .attr('height', 40)
          .attr('patternUnits', 'userSpaceOnUse');
        pattern
          .append('path')
          .attr('d', 'M 40 0 L 0 0 0 40')
          .attr('fill', 'none')
          .attr('stroke', '#f0f0f0')
          .attr('stroke-width', 1);

        mainGroup
          .append('rect')
          .attr('width', DIAGRAM_CONFIG.width)
          .attr('height', DIAGRAM_CONFIG.height)
          .attr('fill', 'url(#grid-pattern)');

        // Content group (this will be transformed)
        const contentGroup = mainGroup.append('g').attr('class', 'content-group');

        // UI group (labels, axes - will have inverse scaling)
        const uiGroup = mainGroup.append('g').attr('class', 'ui-group');

        // Function to render content at current transform
        const render = (transform: d3.ZoomTransform) => {
          const scale = transform.k;
          const invScale = 1 / scale;

          // Clear content
          contentGroup.selectAll('*').remove();
          uiGroup.selectAll('*').remove();

          // Apply transform to content
          contentGroup.attr('transform', transform.toString());

          // Y axis (in UI group, no scaling)
          uiGroup
            .append('line')
            .attr('x1', DIAGRAM_CONFIG.margin.left * scale + transform.x)
            .attr('y1', DIAGRAM_CONFIG.margin.top * scale + transform.y)
            .attr('x2', DIAGRAM_CONFIG.margin.left * scale + transform.x)
            .attr(
              'y2',
              (DIAGRAM_CONFIG.height - DIAGRAM_CONFIG.margin.bottom) * scale +
                transform.y
            )
            .attr('stroke', '#333')
            .attr('stroke-width', 2);

          // Y axis label
          uiGroup
            .append('text')
            .attr('x', 20)
            .attr('y', viewerSize.height / 2)
            .attr('text-anchor', 'middle')
            .attr('font-size', 14)
            .attr('fill', '#666')
            .attr(
              'transform',
              `rotate(-90, 20, ${viewerSize.height / 2})`
            )
            .text('Energy (eV)');

          // HOMO-LUMO gap highlight
          if (orbitalSummary.homoOrbital && orbitalSummary.lumoOrbital) {
            contentGroup
              .append('rect')
              .attr('x', DIAGRAM_CONFIG.margin.left)
              .attr('y', orbitalSummary.lumoOrbital.yPosition)
              .attr('width', chartWidth)
              .attr(
                'height',
                orbitalSummary.homoOrbital.yPosition -
                  orbitalSummary.lumoOrbital.yPosition
              )
              .attr('fill', 'rgba(255, 193, 7, 0.1)')
              .attr('stroke', 'rgba(255, 193, 7, 0.3)')
              .attr('stroke-width', invScale)
              .attr('stroke-dasharray', `${5 * invScale},${5 * invScale}`);
          }

          // Orbital lines positioned close to Y-axis (distance scales inversely with zoom)
          const orbitalX = DIAGRAM_CONFIG.margin.left + 10 * invScale;
          const scaledOrbitalWidth = DIAGRAM_CONFIG.orbitalWidth * invScale;

          const orbitalGroup = contentGroup
            .selectAll('.orbital')
            .data(processedOrbitals)
            .enter()
            .append('g')
            .attr('class', 'orbital');

          // Orbital rectangles with inverse scaling to maintain constant visual size
          orbitalGroup
            .append('rect')
            .attr('x', orbitalX)
            .attr('y', d => d.yPosition - (DIAGRAM_CONFIG.orbitalHeight * invScale) / 2)
            .attr('width', scaledOrbitalWidth)
            .attr('height', DIAGRAM_CONFIG.orbitalHeight * invScale)
            .attr('fill', d => getOrbitalColor(d))
            .attr('stroke', d => getOrbitalColor(d))
            .attr('stroke-width', invScale);

          // Electrons
          orbitalGroup
            .filter(d => d.occupancy > 0)
            .append('circle')
            .attr('cx', orbitalX + scaledOrbitalWidth / 4)
            .attr('cy', d => d.yPosition)
            .attr('r', 3 * invScale)
            .attr('fill', '#34495e');

          orbitalGroup
            .filter(d => d.occupancy > 1)
            .append('circle')
            .attr('cx', orbitalX + (3 * scaledOrbitalWidth) / 4)
            .attr('cy', d => d.yPosition)
            .attr('r', 3 * invScale)
            .attr('fill', '#34495e');

          // Calculate and render labels with collision detection
          const labelPositions = calculateLabelPositions(
            processedOrbitals,
            scale
          );

          labelPositions
            .filter(l => l.visible)
            .forEach(labelInfo => {
              const { orbital, x, y, offsetY } = labelInfo;
              const transformedX = x * scale + transform.x;
              const transformedY = (y + offsetY) * scale + transform.y;

              const isImportant =
                orbital.orbital_type === 'homo' ||
                orbital.orbital_type === 'lumo';

              // Connection line from orbital to label (solid line, starting at orbital edge)
              const orbitalEndX =
                (orbitalX + scaledOrbitalWidth) * scale + transform.x;
              const orbitalY = orbital.yPosition * scale + transform.y;

              uiGroup
                .append('line')
                .attr('x1', orbitalEndX)
                .attr('y1', orbitalY)
                .attr('x2', transformedX - 5 * invScale)
                .attr('y2', transformedY)
                .attr('stroke', '#ccc')
                .attr('stroke-width', 1);

              uiGroup
                .append('text')
                .attr('x', transformedX)
                .attr('y', transformedY + 4)
                .attr('font-size', 12)
                .attr('font-family', 'monospace')
                .attr('fill', isImportant ? '#333' : '#666')
                .attr('font-weight', isImportant ? 'bold' : 'normal')
                .text(
                  `#${orbital.index}: ${orbital.energy_ev.toFixed(3)} eV${
                    orbital.label &&
                    (orbital.label.includes('HOMO') ||
                      orbital.label.includes('LUMO') ||
                      orbital.label.includes('SOMO'))
                      ? ` (${orbital.label})`
                      : ''
                  }`
                );
            });

          // Energy axis ticks
          const tickOrbitals = processedOrbitals.filter(
            (_, index) =>
              index %
                Math.max(1, Math.floor(processedOrbitals.length / 10)) ===
              0
          );

          tickOrbitals.forEach(orbital => {
            const tickY = orbital.yPosition * scale + transform.y;
            const tickX = DIAGRAM_CONFIG.margin.left * scale + transform.x;

            uiGroup
              .append('line')
              .attr('x1', tickX - 5)
              .attr('y1', tickY)
              .attr('x2', tickX)
              .attr('y2', tickY)
              .attr('stroke', '#666')
              .attr('stroke-width', 1);

            uiGroup
              .append('text')
              .attr('x', tickX - 10)
              .attr('y', tickY + 4)
              .attr('text-anchor', 'end')
              .attr('font-size', 10)
              .attr('font-family', 'monospace')
              .attr('fill', '#666')
              .text(orbital.energy_ev.toFixed(1));
          });
        };

        // Initial render
        render(d3.zoomIdentity);

        // Zoom behavior with translation limits
        const zoom = d3
          .zoom<SVGSVGElement, unknown>()
          .scaleExtent([1, 100])
          .translateExtent([
            [0, 0],
            [DIAGRAM_CONFIG.width, DIAGRAM_CONFIG.height],
          ])
          .extent([
            [0, 0],
            [viewerSize.width, viewerSize.height],
          ])
          .on('zoom', event => {
            render(event.transform);
          });

        svg.call(zoom);

        // Double-click to reset
        svg.on('dblclick.zoom', () => {
          svg.transition().duration(300).call(zoom.transform, d3.zoomIdentity);
        });
      }, [
        processedOrbitals,
        orbitalSummary,
        selectedOrbitalIndex,
        viewerSize,
        getOrbitalColor,
        calculateLabelPositions,
        onOrbitalSelect,
      ]);

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

      return (
        <div className={styles.diagramContainer}>
          {/* SVG Diagram */}
          <div ref={containerRef} className={styles.viewerContainer}>
            <svg
              ref={svgRef}
              width={viewerSize.width}
              height={viewerSize.height}
              className={styles.diagramSvg}
            />
          </div>
        </div>
      );
    }
  );
