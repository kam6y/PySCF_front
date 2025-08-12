import React, {
  useRef,
  useEffect,
  useImperativeHandle,
  forwardRef,
} from "react";
import * as $3Dmol from "3dmol";

// 修正点: プロジェクト構成に合わせて型定義ファイルのパスを修正
import { GLViewer, GLModel, StyleSpec } from "../../types/3dmol";

//================================================
// PropsとRefのインターフェース定義
//================================================

export interface MoleculeViewerProps {
  /**
   * ビューワーの幅。'500px'や'80dvw'のようなCSS単位が使えます。
   * @default '600px'
   */
  width?: number | string;
  /**
   * ビューワーの高さ。'500px'や'80dvh'のようなCSS単位が使えます。
   * @default '600px'
   */
  height?: number | string;
  /**
   * 背景色
   * @default 'white'
   */
  backgroundColor?: string;
  /**
   * コンポーネントのルート要素に適用するCSSクラス名
   */
  className?: string;
}

export interface MoleculeViewerRef {
  /** XYZ形式の分子データを読み込みます。 */
  loadXYZ: (xyzData: string) => void;
  /** 表示スタイルを適用します。 */
  setStyle: (style: StyleSpec) => void;
  /** 表示されているモデルをすべて消去します。 */
  clearModels: () => void;
  /** モデル全体が収まるようにズームします。 */
  zoomToFit: () => void;
  /** 現在のビューを画像データとして取得します。 */
  takeScreenshot: () => string;
}

//================================================
// MoleculeViewer コンポーネント本体
//================================================

export const MoleculeViewer = forwardRef<MoleculeViewerRef, MoleculeViewerProps>(
  (
    {
      width = "600px",
      height = "600px",
      backgroundColor = "white",
      className = "",
    },
    ref
  ) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const viewerRef = useRef<GLViewer | null>(null);

    /**
     * 3Dmol.jsのビューワーを初期化し、リサイズを監視するEffect
     */
    useEffect(() => {
      if (!containerRef.current) return;

      let viewer: GLViewer | null = null;
      try {
        viewer = $3Dmol.createViewer(containerRef.current, {
          backgroundColor,
        });
        viewerRef.current = viewer;

        // ResizeObserverでコンテナのサイズ変更を監視し、ビューワーをリサイズする
        const resizeObserver = new ResizeObserver(() => {
          viewer?.resize();
        });
        resizeObserver.observe(containerRef.current);

        // コンポーネントのアンマウント時に監視を停止
        return () => {
          resizeObserver.disconnect();
        };
      } catch (error) {
        console.error("Failed to initialize 3Dmol viewer:", error);
      }
    }, [backgroundColor]); // 背景色の変更時のみ再初期化

    /**
     * 親コンポーネントに公開するメソッドを定義
     */
    useImperativeHandle(ref, () => ({
      loadXYZ: (xyzData: string) => {
        const viewer = viewerRef.current;
        if (!viewer) return;

        try {
          viewer.removeAllModels();
          viewer.addModel(xyzData, "xyz");
          // デフォルトスタイルを設定
          viewer.setStyle(
            {},
            {
              stick: { radius: 0.2 },
              sphere: { radius: 0.3 },
            }
          );
          viewer.zoomTo();
          viewer.render();
        } catch (error) {
          console.error("Failed to load XYZ data:", error);
        }
      },

      setStyle: (style: StyleSpec) => {
        viewerRef.current?.setStyle({}, style);
        viewerRef.current?.render();
      },

      clearModels: () => {
        viewerRef.current?.removeAllModels();
        viewerRef.current?.render();
      },

      zoomToFit: () => {
        viewerRef.current?.zoomTo();
        viewerRef.current?.render();
      },

      takeScreenshot: () => {
        if (viewerRef.current && containerRef.current) {
          const w = containerRef.current.clientWidth;
          const h = containerRef.current.clientHeight;
          return viewerRef.current.screenshot(w, h);
        }
        return "";
      },
    }));

    /**
     * レンダリングするDOM
     */
    return (
      <div
        className={`molecule-viewer ${className}`}
        style={{
          width,
          height,
          border: "1px solid #ccc",
          borderRadius: "8px",
          overflow: "hidden",
          position: "relative",
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