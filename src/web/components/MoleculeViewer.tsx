import React, { useRef, useEffect, useImperativeHandle, forwardRef } from "react";
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
      const initViewer = async () => {
        if (!containerRef.current) return;

        try {
          // Dynamic import to ensure 3Dmol.js is loaded properly
          const $3Dmol = await import("3dmol/build/3Dmol.js") as any;
          const viewer = $3Dmol.createViewer(containerRef.current, {
            backgroundColor,
          });
          
          viewer.setWidth(width);
          viewer.setHeight(height);
          
          viewerRef.current = viewer;
        } catch (error) {
          console.error("Failed to initialize 3Dmol viewer:", error);
        }
      };

      initViewer();

      return () => {
        if (viewerRef.current) {
          viewerRef.current.clear();
        }
      };
    }, [backgroundColor, width, height]);

    useImperativeHandle(ref, () => ({
      loadXYZ: (xyzData: string) => {
        if (!viewerRef.current) {
          console.error("Viewer not initialized");
          return;
        }

        try {
          // Clear existing models
          viewerRef.current.removeAllModels();
          
          // Add new model
          const model = viewerRef.current.addModel(xyzData, "xyz");
          currentModelRef.current = model;
          
          // Set default style
          viewerRef.current.setStyle({}, { 
            stick: { radius: 0.2 },
            sphere: { radius: 0.3 }
          });
          
          // Zoom to fit and render
          viewerRef.current.zoomTo();
          viewerRef.current.render();
        } catch (error) {
          console.error("Failed to load XYZ data:", error);
        }
      },

      setStyle: (style: StyleSpec) => {
        if (!viewerRef.current) {
          console.error("Viewer not initialized");
          return;
        }

        try {
          viewerRef.current.setStyle({}, style);
          viewerRef.current.render();
        } catch (error) {
          console.error("Failed to set style:", error);
        }
      },

      clearModels: () => {
        if (!viewerRef.current) return;
        
        viewerRef.current.removeAllModels();
        viewerRef.current.render();
        currentModelRef.current = null;
      },

      zoomToFit: () => {
        if (!viewerRef.current) return;
        
        viewerRef.current.zoomTo();
        viewerRef.current.render();
      },

      takeScreenshot: () => {
        if (!viewerRef.current) return "";
        
        return viewerRef.current.screenshot(width, height);
      }
    }));

    return (
      <div 
        className={`molecule-viewer ${className}`}
        style={{
          width: width,
          height: height,
          border: "1px solid #ccc",
          borderRadius: "4px",
          overflow: "hidden"
        }}
      >
        <div 
          ref={containerRef}
          style={{
            width: "100%",
            height: "100%",
            position: "relative"
          }}
        />
      </div>
    );
  }
);