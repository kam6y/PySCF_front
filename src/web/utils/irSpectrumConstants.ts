/**
 * IRスペクトル表示とデータ生成に関する共有定数
 *
 * このファイルは IRSpectrumChart と CalculationResultsPage で共有される
 * デフォルト設定値と制約条件を一元管理します。
 */

/**
 * IRスペクトルのデフォルト設定値
 */
export const IR_SPECTRUM_DEFAULTS = {
  /** ローレンツ型ブロードニング関数の半値全幅 (cm⁻¹) */
  broadening_fwhm: 20.0,
  /** X軸（波数）の最小値 (cm⁻¹) */
  x_min: 400.0,
  /** X軸（波数）の最大値 (cm⁻¹) */
  x_max: 4000.0,
  /**
   * チャート上のピークマーカーの表示/非表示
   * 注: この設定はIRスペクトルチャート上のマーカーのみを制御します。
   * VibrationModeViewerのピークテーブルは常に表示されます（振動モード選択のため）。
   */
  show_peaks: false,
};

/**
 * サーバーにリクエストする固定範囲
 *
 * ユーザーの設定に関わらず常にこの範囲でデータを取得し、
 * クライアント側でフィルタリングすることでパフォーマンスを向上させる
 */
export const IR_SPECTRUM_API_RANGE = {
  /** APIリクエストの最小波数 (cm⁻¹) */
  x_min: 0.0,
  /** APIリクエストの最大波数 (cm⁻¹) */
  x_max: 4500.0,
} as const;

/**
 * IRスペクトル設定の制約条件
 */
export const IR_SPECTRUM_CONSTRAINTS = {
  broadening_fwhm: {
    min: 0.1,
    max: 1000,
    step: 10,
  },
  x_min: {
    min: 0,
    max: 4500,
    step: 100,
  },
  x_max: {
    min: 0,
    max: 4500,
    step: 100,
  },
} as const;

/**
 * IRスペクトル設定の型定義
 *
 * すべてのIR設定を一元管理するための型
 */
export type IRSettings = {
  broadening_fwhm: number;
  x_min: number;
  x_max: number;
  show_peaks: boolean;
};
