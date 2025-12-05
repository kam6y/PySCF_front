/**
 * スプラッシュスクリーン関連の型定義
 */

/**
 * スプラッシュスクリーンのステージ
 */
export type SplashStage =
  | 'initializing'
  | 'detecting-env'
  | 'finding-port'
  | 'starting-server'
  | 'health-check'
  | 'creating-window'
  | 'closing';

/**
 * スプラッシュスクリーンの状態更新データ
 */
export interface SplashStatusUpdate {
  /** 現在のステージ */
  stage: SplashStage;
  /** 表示するメッセージ */
  message: string;
  /** リトライカウント（ヘルスチェック時のみ） */
  retryCount?: number;
}

/**
 * Splash Renderer Process用のAPI
 * preloadスクリプトからcontextBridgeで公開される
 */
export interface SplashAPI {
  /**
   * 進捗状態の更新を受信
   * @param callback 状態更新時に呼び出されるコールバック関数
   * @returns クリーンアップ関数
   */
  onUpdateStatus: (callback: (update: SplashStatusUpdate) => void) => () => void;

  /**
   * エラーメッセージの受信
   * @param callback エラー発生時に呼び出されるコールバック関数
   * @returns クリーンアップ関数
   */
  onShowError: (callback: (message: string) => void) => () => void;

  /**
   * スプラッシュクローズ指示の受信
   * @param callback クローズ指示時に呼び出されるコールバック関数
   * @returns クリーンアップ関数
   */
  onClose: (callback: () => void) => () => void;
}

/**
 * グローバルにsplashAPIを公開
 */
declare global {
  interface Window {
    splashAPI: SplashAPI;
  }
}
