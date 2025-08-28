// src/web/App.tsx

import React, { useState, useMemo } from 'react';
import './App.css';
import styles from './App.module.css';
import { Header } from './components/Header';
import { Sidebar } from './components/Sidebar';
import { ToastContainer } from './components/ToastContainer';
import { CalculationSettingsPage } from './pages/CalculationSettingsPage';
import { CalculationResultsPage } from './pages/CalculationResultsPage';
import { DrawMoleculePage } from './pages/DrawMoleculePage';
import {
  useSidebarState,
  usePageNavigation,
  useActiveCalculation,
  useCalculationActions,
  useStagedCalculation,
  useCalculationWebSocket,
} from './hooks';
import { CalculationInstance } from './types/api-types';

export const App = () => {
  // 新しい状態管理
  const [searchQuery, setSearchQuery] = useState('');
  const [isUserMenuOpen, setIsUserMenuOpen] = useState(false);

  // 統一された状態管理フック
  const sidebarState = useSidebarState();
  const pageNavigation = usePageNavigation();
  const calculationActions = useCalculationActions();
  const stagedCalculation = useStagedCalculation();

  // 統一されたアクティブ計算状態（複雑な導出ロジックが内部で処理される）
  const {
    activeCalculation,
    activeCalculationId,
    isStagedCalculation,
    selectCalculation,
    isLoading: calculationLoading,
    detailsLoading,
    calculationsLoading,
    calculationsError,
    sidebarCalculations,
  } = useActiveCalculation();

  // WebSocketによるリアルタイム更新（簡素化されたインターフェース）
  useCalculationWebSocket(
    activeCalculation?.id || null,
    activeCalculation?.status
  );

  // 検索機能によるフィルタリング
  const filteredCalculations = useMemo(() => {
    if (!searchQuery.trim()) {
      return sidebarCalculations;
    }

    const query = searchQuery.toLowerCase();
    return sidebarCalculations.filter(
      calc =>
        calc.name.toLowerCase().includes(query) ||
        calc.status.toLowerCase().includes(query)
    );
  }, [sidebarCalculations, searchQuery]);

  // イベントハンドラー
  const handleSearchChange = (query: string) => {
    setSearchQuery(query);
  };

  const handleCreateNew = () => {
    stagedCalculation.createNewCalculation(
      pageNavigation.handleDropdownOptionSelect,
      sidebarState.handleSidebarClose
    );
  };

  const handleUserMenuToggle = () => {
    setIsUserMenuOpen(prev => !prev);
  };

  const renderCurrentPage = () => {
    switch (pageNavigation.currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation || undefined}
            onCalculationUpdate={
              isStagedCalculation
                ? stagedCalculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={
              stagedCalculation.createNewFromExisting
            }
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage
            activeCalculation={activeCalculation || undefined}
            isLoadingDetails={detailsLoading}
            detailsError={null}
            onCalculationUpdate={calculationActions.handleCalculationUpdate}
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation || undefined}
            onCalculationUpdate={
              isStagedCalculation
                ? stagedCalculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={
              stagedCalculation.createNewFromExisting
            }
          />
        );
    }
  };

  return (
    <div className={styles.appContainer}>
      <button
        className={`${styles.independentSidebarToggle} ${sidebarState.isSidebarOpen ? styles.sidebarOpen : ''}`}
        onClick={sidebarState.handleSidebarToggle}
        aria-label="Toggle sidebar"
      >
        <svg
          width="16"
          height="16"
          viewBox="0 0 16 16"
          fill="none"
          xmlns="http://www.w3.org/2000/svg"
        >
          {sidebarState.isSidebarOpen ? (
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
        onDropdownToggle={sidebarState.handleDropdownToggle}
        onPlusClick={() =>
          stagedCalculation.createNewCalculation(
            pageNavigation.handleDropdownOptionSelect,
            sidebarState.handleSidebarClose
          )
        }
        isDropdownOpen={sidebarState.isDropdownOpen}
        isSidebarOpen={sidebarState.isSidebarOpen}
        currentPageTitle={pageNavigation.currentPageTitle}
        currentPage={pageNavigation.currentPage}
        onDropdownOptionSelect={pageNavigation.handleDropdownOptionSelect}
        onDropdownClose={sidebarState.handleDropdownClose}
      />

      <Sidebar
        isOpen={sidebarState.isSidebarOpen}
        onClose={sidebarState.handleSidebarClose}
        calculations={filteredCalculations}
        activeCalculationId={activeCalculation?.id || null}
        calculationsLoading={calculationsLoading}
        calculationsError={calculationsError ? calculationsError.message : null}
        onCalculationSelect={(calculationId: string) => {
          selectCalculation(calculationId);
          stagedCalculation.clearStaged();
          sidebarState.handleSidebarClose();
        }}
        onCalculationDelete={async (calculationId: string) => {
          await calculationActions.handleCalculationDelete(calculationId);
          // 削除された計算がアクティブだった場合はクリア
          if (activeCalculationId === calculationId) {
            selectCalculation(null);
            stagedCalculation.clearStaged();
            pageNavigation.handleDropdownOptionSelect('calculation-settings');
          }
        }}
        onCreateNew={handleCreateNew}
        searchQuery={searchQuery}
        onSearchChange={handleSearchChange}
        onUserMenuToggle={handleUserMenuToggle}
        isUserMenuOpen={isUserMenuOpen}
      />

      <main
        className={`${styles.appContent} ${sidebarState.isSidebarOpen ? styles.sidebarOpen : ''}`}
      >
        {renderCurrentPage()}
      </main>

      <ToastContainer />
    </div>
  );
};
