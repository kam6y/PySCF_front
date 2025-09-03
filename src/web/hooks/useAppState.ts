import { useUIStore } from '../store/uiStore';
import { useCalculationStore } from '../store/calculationStore';

/**
 * UIStoreとCalculationStoreを統合したアプリケーション状態へのアクセス
 * 
 * このフックは、UIとCalculation両方の状態とアクションを
 * 単一のインターフェースで提供し、App.tsxの複雑さを軽減
 */
export const useAppState = () => {
  const uiState = useUIStore();
  const calculationState = useCalculationStore();

  return {
    // UI状態
    ui: {
      // Sidebar
      isSidebarOpen: uiState.isSidebarOpen,
      isDropdownOpen: uiState.isDropdownOpen,
      toggleSidebar: uiState.toggleSidebar,
      closeSidebar: uiState.closeSidebar,
      toggleDropdown: uiState.toggleDropdown,
      closeDropdown: uiState.closeDropdown,
      
      // Page navigation
      currentPage: uiState.currentPage,
      currentPageTitle: uiState.getCurrentPageTitle(),
      setCurrentPage: uiState.setCurrentPage,
      
      // Search
      searchQuery: uiState.searchQuery,
      setSearchQuery: uiState.setSearchQuery,
      
      // Modals
      isUserMenuOpen: uiState.isUserMenuOpen,
      isSettingsOpen: uiState.isSettingsOpen,
      toggleUserMenu: uiState.toggleUserMenu,
      closeUserMenu: uiState.closeUserMenu,
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
      clearStaged: calculationState.clearStaged,
      isStagedCalculation: calculationState.isStagedCalculation,
    },

    // 統合されたアクション（複数ストアにまたがる操作）
    actions: {
      // 新規計算作成（UI状態も更新）
      handleCreateNew: calculationState.createNewCalculation,
      
      // 計算選択（UI状態も更新）
      handleCalculationSelect: (calculationId: string) => {
        calculationState.selectCalculation(calculationId);
        calculationState.clearStaged();
        uiState.closeSidebar();
      },
    },
  };
};