import { create } from 'zustand';
import { CalculationInstance } from '../types/api-types';

interface CalculationState {
  // アクティブ計算ID管理
  activeCalculationId: string | null;
  setActiveCalculationId: (id: string | null) => void;

  // 新規計算作成時の一時状態管理
  stagedCalculation: CalculationInstance | null;
  setStagedCalculation: (calc: CalculationInstance | null) => void;
  clearStagedCalculation: () => void;
}

export const useCalculationStore = create<CalculationState>(set => ({
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
}));
