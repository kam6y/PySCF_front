import React, { useRef, useEffect, useImperativeHandle, forwardRef } from "react";
import * as $3Dmol from '3dmol'; // 変更点: 型付きでインポート
import { GLViewer, GLModel, StyleSpec } from "../../types/3dmol";

export interface MoleculeViewerProps {
  width?: number;
  height?: number;
  backgroundColor?: string;
  className?: string;
}

export interface MoleculeViewerRef {
  loadXYZ: (xyzData: string) => void;
  setStyle: (style: StyleSpec) => void;
  clearModels: () => void;
  zoomToFit: () => void;
  takeScreenshot: () => string;
}

export const MoleculeViewer = forwardRef<MoleculeViewerRef, MoleculeViewerProps>(
  ({ width = 600, height = 600, backgroundColor = "white", className = "" }, ref) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const viewerRef = useRef<GLViewer | null>(null);
    const currentModelRef = useRef<GLModel | null>(null);

    useEffect(() => {
      if (!containerRef.current) return;

      try {
        // 変更点: 'as any' を削除
        const viewer = $3Dmol.createViewer(containerRef.current, {
          backgroundColor,
        });
        
        viewer.setWidth(width);
        viewer.setHeight(height);
        
        viewerRef.current = viewer;
      } catch (error) {
        console.error("Failed to initialize 3Dmol viewer:", error);
      }

      return () => {
        if (viewerRef.current) {
          // You might need to add a proper cleanup method if available in 3Dmol.js
          // For now, clear() is a good approach.
        }
      };
    }, [backgroundColor, width, height]);

    useImperativeHandle(ref, () => ({
      loadXYZ: (xyzData: string) => {
        const viewer = viewerRef.current;
        if (!viewer) {
          console.error("Viewer not initialized");
          return;
        }

        try {
          viewer.removeAllModels();
          const model = viewer.addModel(xyzData, "xyz");
          currentModelRef.current = model;
          
          // Set default style
          viewer.setStyle({}, {
            stick: { radius: 0.2 },
            sphere: { radius: 0.3 }
          });
          
          viewer.zoomTo();
          viewer.render();
        } catch (error) {
          console.error("Failed to load XYZ data:", error);
        }
      },

      setStyle: (style: StyleSpec) => {
        const viewer = viewerRef.current;
        if (!viewer) {
          console.error("Viewer not initialized");
          return;
        }

        try {
          viewer.setStyle({}, style);
          viewer.render();
        } catch (error) {
          console.error("Failed to set style:", error);
        }
      },

      clearModels: () => {
        if (viewerRef.current) {
          viewerRef.current.removeAllModels();
          viewerRef.current.render();
          currentModelRef.current = null;
        }
      },

      zoomToFit: () => {
        if (viewerRef.current) {
          viewerRef.current.zoomTo();
          viewerRef.current.render();
        }
      },

      takeScreenshot: () => {
        if (viewerRef.current) {
          return viewerRef.current.screenshot(width, height);
        }
        return "";
      }
    }));

    return (
      <div 
        className={`molecule-viewer ${className}`}
        style={{
          width,
          height,
          border: "1px solid #ccc",
          borderRadius: "4px",
          overflow: "hidden",
          position: "relative"
        }}
      >
        <div 
          ref={containerRef}
          style={{
            width: "100%",
            height: "100%",
          }}
        />
      </div>
    );
  }
);