import { QuantumCalculationRequest } from '../types/api-types';

// DrawMoleculePage / calculationStore で使用するデフォルト計算パラメータ
export const DEFAULT_CALCULATION_PARAMETERS: QuantumCalculationRequest = {
  calculation_method: 'DFT',
  basis_function: '6-31G(d)',
  exchange_correlation: 'B3LYP',
  charges: 0,
  spin: 0,
  solvent_method: 'none',
  solvent: '-',
  xyz: '',
  name: '',
  tddft_nstates: 10,
  tddft_method: 'TDDFT',
  tddft_analyze_nto: false,
  ncas: 4,
  nelecas: 4,
  max_cycle_macro: 50,
  max_cycle_micro: 4,
  natorb: true,
  conv_tol: 1e-6,
  conv_tol_grad: 1e-4,
  density_fitting: false,
  optimize_geometry: true,
  geomopt_maxsteps: 100,
  geomopt_conv_energy: 1e-6,
} as QuantumCalculationRequest;

// DrawMoleculePage のステータスメッセージ
export const STATUS_MESSAGES: Record<string, string> = {
  running: 'Calculation is running. The molecule structure is read-only.',
  waiting:
    'Calculation is waiting for resources. The molecule structure is read-only.',
  completed:
    'Calculation completed. Edit the structure to create a new calculation.',
  error:
    'Previous calculation had errors. Edit the structure to create a new calculation.',
};

// CalculationSettingsPage の計算ボタンテキスト
export const CALCULATION_BUTTON_TEXT: Record<string, string> = {
  running: 'Running...',
  waiting: 'Waiting...',
  completed: 'Completed!',
  error: 'Error',
};

// CalculationSettingsPage の入力プレースホルダー
export const INPUT_PLACEHOLDERS: Record<string, string> = {
  smiles: 'e.g., CCO for ethanol',
  pubchem: 'e.g., aspirin, or 2244',
};

// CalculationSettingsPage の定義済み溶媒リスト
export const PREDEFINED_SOLVENTS = [
  'water',
  'dimethylsulfoxide',
  'n,n-dimethylformamide',
  'nitromethane',
  'methanol',
  'ethanol',
  'acetone',
  'dichloroethane',
  'dichloromethane',
  'tetrahydrofuran',
  'chlorobenzene',
  'chloroform',
  'diethylether',
  'toluene',
  'benzene',
  '1,4-dioxane',
  'cyclohexane',
  'custom',
] as const;

// 溶媒がカスタム誘電率かどうかを判定する純粋関数
export const isCustomDielectricConstant = (
  solventValue: string | undefined
): boolean => {
  if (!solventValue || solventValue === '-') return false;
  if (
    PREDEFINED_SOLVENTS.includes(
      solventValue.toLowerCase() as (typeof PREDEFINED_SOLVENTS)[number]
    )
  )
    return false;
  const numValue = parseFloat(solventValue);
  return !isNaN(numValue) && numValue > 0;
};
