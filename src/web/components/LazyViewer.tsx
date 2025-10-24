import React, { useRef, useEffect, useState, ReactNode } from 'react';

interface LazyViewerProps {
  children: ReactNode;
  placeholder?: ReactNode;
  /**
   * IntersectionObserverのrootMarginオプション
   * デフォルト: '200px' (ビューポートから200px手前で初期化開始)
   */
  rootMargin?: string;
}

/**
 * LazyViewer - 遅延レンダリングラッパーコンポーネント
 *
 * IntersectionObserverを使用して、コンポーネントが画面内に入った（または近づいた）
 * ときのみ子コンポーネントをレンダリングします。
 *
 * これにより、重い3Dビューアーコンポーネントの初期化を遅延させ、
 * ページの初期レンダリング時間を大幅に短縮できます。
 */
export const LazyViewer: React.FC<LazyViewerProps> = ({
  children,
  placeholder,
  rootMargin = '200px',
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const [isVisible, setIsVisible] = useState(false);
  const [hasBeenVisible, setHasBeenVisible] = useState(false);

  useEffect(() => {
    const currentContainer = containerRef.current;
    if (!currentContainer) return;

    // IntersectionObserverのコールバック
    const handleIntersection = (entries: IntersectionObserverEntry[]) => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          setIsVisible(true);
          // 一度表示されたら、その状態を保持
          setHasBeenVisible(true);
        } else {
          setIsVisible(false);
        }
      });
    };

    // IntersectionObserverを作成
    const observer = new IntersectionObserver(handleIntersection, {
      root: null, // ビューポートを基準
      rootMargin, // 指定されたマージン
      threshold: 0, // 1pxでも見えたら
    });

    observer.observe(currentContainer);

    return () => {
      observer.disconnect();
    };
  }, [rootMargin]);

  // 一度でも表示されたら、その後は常にレンダリング
  // （スクロールで外れてもアンマウントしない）
  const shouldRenderChildren = hasBeenVisible;

  return (
    <div
      ref={containerRef}
      style={{ minHeight: shouldRenderChildren ? 'auto' : '400px' }}
    >
      {shouldRenderChildren
        ? children
        : placeholder || (
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                minHeight: '400px',
                backgroundColor: '#f5f5f5',
                borderRadius: '8px',
                color: '#666',
              }}
            >
              <div>
                <div style={{ fontSize: '24px', marginBottom: '8px' }}>⚛️</div>
                <div>Loading visualization...</div>
              </div>
            </div>
          )}
    </div>
  );
};
