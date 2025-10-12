// src/web/App.tsx

import React, { useMemo } from 'react';
import './App.css';
import styles from './App.module.css';
import { Header } from './components/Header';
import { Sidebar } from './components/Sidebar';
import { ToastContainer } from './components/ToastContainer';
import { CalculationSettingsPage } from './pages/CalculationSettingsPage';
import { CalculationResultsPage } from './pages/CalculationResultsPage';
import { SettingsPage } from './pages/SettingsPage';
import { AgentPage } from './pages/AgentPage';
import { useAppState } from './hooks/useAppState';
import { useCalculationData } from './hooks/useCalculationData';
import { useCalculationActions } from './hooks/useCalculationActions';
import { useUnifiedWebSocket } from './hooks/useUnifiedWebSocket';

export const App = () => {
  // 統合された状態管理
  const appState = useAppState();
  const calculationData = useCalculationData();
  const calculationActions = useCalculationActions();

  // 統合WebSocketによるリアルタイム更新（グローバル + アクティブ計算監視）
  useUnifiedWebSocket({
    activeCalculationId: calculationData.activeCalculation?.id || null,
  });

  // 検索機能によるフィルタリング
  const filteredCalculations = useMemo(() => {
    const searchQuery = appState.ui.searchQuery;
    if (!searchQuery.trim()) {
      return calculationData.sidebarCalculations;
    }

    const query = searchQuery.toLowerCase();
    return calculationData.sidebarCalculations.filter(
      calc =>
        calc.name.toLowerCase().includes(query) ||
        calc.status.toLowerCase().includes(query)
    );
  }, [calculationData.sidebarCalculations, appState.ui.searchQuery]);

  // イベントハンドラー（統合されたアクションを使用）
  const handleSearchChange = (query: string) => {
    appState.ui.setSearchQuery(query);
  };

  const handleCreateNew = () => {
    appState.calculation.createNewCalculation();
  };

  const handleCalculationSelect = (calculationId: string) => {
    appState.actions.handleCalculationSelect(calculationId);
  };

  const renderCurrentPage = () => {
    // AI Agentが有効な場合はAgentPageを表示
    if (appState.ui.isAIAgentEnabled) {
      return <AgentPage />;
    }

    switch (appState.ui.currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage
            activeCalculation={calculationData.activeCalculation || undefined}
            onCalculationUpdate={
              calculationData.isStagedCalculation
                ? appState.calculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={
              appState.calculation.createNewFromExisting
            }
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage
            activeCalculation={calculationData.activeCalculation || undefined}
            isLoadingDetails={calculationData.detailsLoading}
            detailsError={null}
            onCalculationUpdate={calculationActions.handleCalculationUpdate}
          />
        );
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={calculationData.activeCalculation || undefined}
            onCalculationUpdate={
              calculationData.isStagedCalculation
                ? appState.calculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={
              appState.calculation.createNewFromExisting
            }
          />
        );
    }
  };

  return (
    <div className={styles.appContainer}>
      <button
        className={`${styles.independentSidebarToggle} ${appState.ui.isSidebarOpen ? styles.sidebarOpen : ''}`}
        onClick={appState.ui.toggleSidebar}
        aria-label="Toggle sidebar"
      >
        <svg
          width="16"
          height="16"
          viewBox="0 0 16 16"
          fill="none"
          xmlns="http://www.w3.org/2000/svg"
        >
          {appState.ui.isSidebarOpen ? (
            <path
              d="M10 4L6 8L10 12"
              stroke="currentColor"
              strokeWidth="3"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          ) : (
            <path
              d="M6 4L10 8L6 12"
              stroke="currentColor"
              strokeWidth="3"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          )}
        </svg>
      </button>

      <Header
        onDropdownToggle={appState.ui.toggleDropdown}
        onPlusClick={handleCreateNew}
        isDropdownOpen={appState.ui.isDropdownOpen}
        isSidebarOpen={appState.ui.isSidebarOpen}
        currentPageTitle={appState.ui.currentPageTitle}
        currentPage={appState.ui.currentPage}
        onDropdownOptionSelect={appState.ui.setCurrentPage}
        onDropdownClose={appState.ui.closeDropdown}
        isAIAgentEnabled={appState.ui.isAIAgentEnabled}
        onAIAgentToggle={appState.ui.setAIAgentEnabled}
      />

      <Sidebar
        isOpen={appState.ui.isSidebarOpen}
        onClose={appState.ui.closeSidebar}
        calculations={filteredCalculations}
        activeCalculationId={calculationData.activeCalculation?.id || null}
        calculationsLoading={calculationData.calculationsLoading}
        calculationsError={
          calculationData.calculationsError
            ? calculationData.calculationsError.message
            : null
        }
        onCalculationSelect={handleCalculationSelect}
        onCalculationDelete={async (calculationId: string) => {
          await calculationActions.handleCalculationDelete(calculationId);
          // 削除された計算がアクティブだった場合はクリア
          if (calculationData.activeCalculationId === calculationId) {
            appState.calculation.selectCalculation(null);
            appState.calculation.clearStaged();
            appState.ui.setCurrentPage('calculation-settings');
          }
        }}
        onCreateNew={handleCreateNew}
        searchQuery={appState.ui.searchQuery}
        onSearchChange={handleSearchChange}
        onUserMenuToggle={appState.ui.toggleUserMenu}
        isUserMenuOpen={appState.ui.isUserMenuOpen}
        onSettingsOpen={appState.ui.openSettings}
      />

      <main
        className={`${styles.appContent} ${appState.ui.isSidebarOpen ? styles.sidebarOpen : ''}`}
      >
        {renderCurrentPage()}
      </main>

      {/* Settings Overlay */}
      {appState.ui.isSettingsOpen && (
        <div className={styles.settingsOverlay}>
          <div
            className={styles.settingsModalBackdrop}
            onClick={appState.ui.closeSettings}
          />
          <div className={styles.settingsModal}>
            <div className={styles.settingsModalHeader}>
              <button
                onClick={appState.ui.closeSettings}
                className={styles.settingsCloseButton}
                aria-label="Close settings"
              >
                <svg
                  width="20"
                  height="20"
                  viewBox="0 0 20 20"
                  fill="none"
                  xmlns="http://www.w3.org/2000/svg"
                >
                  <path
                    d="M15 5L5 15M5 5L15 15"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </svg>
              </button>
            </div>
            <SettingsPage />
          </div>
        </div>
      )}

      <ToastContainer />
    </div>
  );
};
