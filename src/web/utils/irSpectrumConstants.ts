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
  /** ピークマーカーの表示/非表示 */
  show_peaks: true,
};

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
    max: 10000,
    step: 100,
  },
  x_max: {
    min: 0,
    max: 10000,
    step: 100,
  },
} as const;

/**
 * IRスペクトル設定の型定義
 */
export type IRSettings = {
  broadening_fwhm: number;
  x_min: number;
  x_max: number;
  show_peaks: boolean;
};

/**
 * 部分的なIR設定の型（VibrationModeViewer と CalculationResultsPage で使用）
 */
export type PartialIRSettings = Pick<
  IRSettings,
  'x_min' | 'x_max' | 'show_peaks'
>;
