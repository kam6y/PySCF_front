import { useUIStore } from '../store/uiStore';
import { useCalculationStore } from '../store/calculationStore';
import { useShallow } from 'zustand/react/shallow';

/**
 * UIStoreとCalculationStoreを統合したアプリケーション状態へのアクセス
 *
 * このフックは、UIとCalculation両方の状態とアクションを
 * 単一のインターフェースで提供し、App.tsxの複雑さを軽減
 */
export const useAppState = () => {
  const uiState = useUIStore(useShallow(state => state));
  const calculationState = useCalculationStore(useShallow(state => state));

  return {
    // UI状態
    ui: {
      // Sidebar
      isSidebarOpen: uiState.isSidebarOpen,
      isDropdownOpen: uiState.isDropdownOpen,
      sidebarView: uiState.sidebarView,
      toggleSidebar: uiState.toggleSidebar,
      closeSidebar: uiState.closeSidebar,
      toggleDropdown: uiState.toggleDropdown,
      closeDropdown: uiState.closeDropdown,
      setSidebarView: uiState.setSidebarView,

      // Page navigation
      currentPage: uiState.currentPage,
      currentPageTitle: uiState.getCurrentPageTitle(),
      setCurrentPage: uiState.setCurrentPage,

      // AI Agent
      isAIAgentEnabled: uiState.isAIAgentEnabled,
      setAIAgentEnabled: uiState.setAIAgentEnabled,

      // Search
      searchQuery: uiState.searchQuery,
      setSearchQuery: uiState.setSearchQuery,

      // Modals
      isUserMenuOpen: uiState.isUserMenuOpen,
      isSettingsOpen: uiState.isSettingsOpen,
      toggleUserMenu: uiState.toggleUserMenu,
      openSettings: uiState.openSettings,
      closeSettings: uiState.closeSettings,
    },

    // Calculation状態
    calculation: {
      // Basic state
      activeCalculationId: calculationState.activeCalculationId,
      stagedCalculation: calculationState.stagedCalculation,

      // Actions
      selectCalculation: calculationState.selectCalculation,
      createNewCalculation: calculationState.createNewCalculation,
      createNewFromExisting: calculationState.createNewFromExisting,
      updateStagedCalculation: calculationState.updateStagedCalculation,
      clearStagedCalculation: calculationState.clearStagedCalculation,
    },

    // 統合されたアクション（複数ストアにまたがる操作）
    actions: {
      // 新規計算作成（UI状態も更新、Agent画面をオフ）
      handleCreateNew: () => {
        calculationState.createNewCalculation();
        uiState.setCurrentPage('calculation-settings');
        uiState.closeSidebar();
        uiState.setAIAgentEnabled(false);
      },

      // 計算選択（UI状態も更新、Agent画面をオフ）
      handleCalculationSelect: (calculationId: string) => {
        calculationState.selectCalculation(calculationId);
        uiState.closeSidebar();
        uiState.setAIAgentEnabled(false);
      },

      // ドロップダウントグル（Agent画面をオフ）
      handleDropdownToggle: () => {
        uiState.toggleDropdown();
        uiState.setAIAgentEnabled(false);
      },
    },
  };
};
