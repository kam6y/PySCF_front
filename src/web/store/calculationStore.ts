// src/web/store/calculationStore.ts

import { create } from 'zustand';
import { CalculationInstance } from '../types/api-types';

interface CalculationState {
  activeCalculationId: string | null;
  stagedCalculation: CalculationInstance | null;
  setActiveCalculationId: (id: string | null) => void;
  setStagedCalculation: (calc: CalculationInstance | null) => void;
  clearStagedCalculation: () => void;
}

export const useCalculationStore = create<CalculationState>(set => ({
  activeCalculationId: null,
  stagedCalculation: null,

  setActiveCalculationId: id =>
    set(state => {
      // IDが切り替わったら編集中の一時データはクリア（新規計算IDは除く）
      if (state.activeCalculationId !== id) {
        // 新規計算IDの場合はstagedCalculationをクリアしない
        if (id && id.startsWith('new-calculation-')) {
          return { activeCalculationId: id };
        }
        return { activeCalculationId: id, stagedCalculation: null };
      }
      return { activeCalculationId: id };
    }),

  setStagedCalculation: calc => set({ stagedCalculation: calc }),

  clearStagedCalculation: () => set({ stagedCalculation: null }),
}));
