import { create } from 'zustand';
import {
  CalculationInstance,
  QuantumCalculationRequest,
} from '../types/api-types';
import { useUIStore } from './uiStore';

interface CalculationState {
  // アクティブ計算ID管理
  activeCalculationId: string | null;
  setActiveCalculationId: (id: string | null) => void;

  // 新規計算作成時の一時状態管理
  stagedCalculation: CalculationInstance | null;
  setStagedCalculation: (calc: CalculationInstance | null) => void;
  clearStagedCalculation: () => void;

  // 統合されたアクション
  selectCalculation: (calculationId: string | null) => void;
  createNewCalculation: () => void;
  createNewFromExisting: (
    originalCalc: CalculationInstance,
    newParams: QuantumCalculationRequest
  ) => void;
  updateStagedCalculation: (updatedCalculation: CalculationInstance) => void;
  clearStaged: () => void;
  isStagedCalculation: (id: string | null) => boolean;
}

export const useCalculationStore = create<CalculationState>((set, get) => ({
  // アクティブ計算ID
  activeCalculationId: null,
  setActiveCalculationId: (id: string | null) =>
    set(state => {
      // IDが切り替わったとき、新規計算IDでない場合はstagedCalculationをクリア
      if (
        state.activeCalculationId !== id &&
        id &&
        !id.startsWith('new-calculation-')
      ) {
        return { activeCalculationId: id, stagedCalculation: null };
      }
      return { activeCalculationId: id };
    }),

  // 新規計算の一時状態
  stagedCalculation: null,
  setStagedCalculation: (calc: CalculationInstance | null) =>
    set({ stagedCalculation: calc }),
  clearStagedCalculation: () => set({ stagedCalculation: null }),

  // 統合されたアクション実装
  selectCalculation: (calculationId: string | null) => {
    set({ activeCalculationId: calculationId });
    if (calculationId) {
      // 既存計算を選択した場合はstagedをクリア
      set({ stagedCalculation: null });
    }
  },

  createNewCalculation: () => {
    const defaultParams: QuantumCalculationRequest = {
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
    };

    const newId = `new-calculation-${Date.now()}`;
    const newCalculation: CalculationInstance = {
      id: newId,
      name: '',
      status: 'pending',
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      parameters: defaultParams,
      results: undefined,
    };

    set({
      stagedCalculation: newCalculation,
      activeCalculationId: newId,
    });

    // UIStoreのアクションを呼び出し
    useUIStore.getState().setCurrentPage('calculation-settings');
    useUIStore.getState().closeSidebar();
  },

  createNewFromExisting: (
    originalCalc: CalculationInstance,
    newParams: QuantumCalculationRequest
  ) => {
    const newId = `new-calculation-${Date.now()}`;
    const newCalculation: CalculationInstance = {
      ...originalCalc,
      id: newId,
      status: 'pending',
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      parameters: newParams,
      results: undefined,
      workingDirectory: undefined,
      errorMessage: undefined,
    };

    set({
      stagedCalculation: newCalculation,
      activeCalculationId: newId,
    });
  },

  updateStagedCalculation: (updatedCalculation: CalculationInstance) => {
    const { stagedCalculation } = get();
    if (stagedCalculation && updatedCalculation.id === stagedCalculation.id) {
      set({ stagedCalculation: updatedCalculation });
    }
  },

  clearStaged: () => {
    set({ stagedCalculation: null });
  },

  isStagedCalculation: (id: string | null) => {
    return id ? id.startsWith('new-calculation-') : false;
  },
}));
