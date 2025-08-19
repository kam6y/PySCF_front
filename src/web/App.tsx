// src/web/App.tsx

import { useQueryClient } from '@tanstack/react-query';
import './App.css';
import { Header } from './components/Header';
import { Sidebar } from './components/Sidebar';
import { CalculationSettingsPage } from './pages/CalculationSettingsPage';
import { CalculationResultsPage } from './pages/CalculationResultsPage';
import { DrawMoleculePage } from './pages/DrawMoleculePage';
import { useCalculationSubscription } from './hooks/useCalculationSubscription';
import {
  useSidebarState,
  usePageNavigation,
  useCalculationOperations,
  useActiveCalculation,
} from './hooks';
import { CalculationInstance } from './types/api-types';

export const App = () => {
  const queryClient = useQueryClient();

  // 責務分散されたカスタムフック
  const sidebarState = useSidebarState();
  const pageNavigation = usePageNavigation();
  const calculationOperations = useCalculationOperations(
    pageNavigation.handleDropdownOptionSelect
  );
  const activeCalculationData = useActiveCalculation(
    pageNavigation.handleDropdownOptionSelect,
    sidebarState.handleSidebarClose
  );

  // WebSocketによるリアルタイム更新
  useCalculationSubscription({
    calculationId: activeCalculationData.activeCalculation?.id || null,
    status: activeCalculationData.activeCalculation?.status,
    onUpdate: (updatedCalculation: CalculationInstance) => {
      // 個別計算詳細のキャッシュを更新
      queryClient.setQueryData(['calculation', updatedCalculation.id], {
        calculation: updatedCalculation,
      });

      // 計算リストのキャッシュも直接更新して即座に反映
      queryClient.setQueryData(['calculations'], (oldData: any) => {
        if (!oldData?.calculations) return oldData;

        const updatedCalculations = oldData.calculations.map((calc: any) =>
          calc.id === updatedCalculation.id
            ? { ...calc, status: updatedCalculation.status }
            : calc
        );

        return {
          ...oldData,
          calculations: updatedCalculations,
        };
      });

      // 念のため計算リストも無効化（バックグラウンドでの最新データ取得用）
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    },
    onError: (error: string) => {
      console.error('WebSocket error:', error);
    },
  });

  const renderCurrentPage = () => {
    switch (pageNavigation.currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculationData.activeCalculation}
            onCalculationUpdate={
              activeCalculationData.handleActiveCalculationUpdate
            }
            onStartCalculation={calculationOperations.handleStartCalculation}
            onCalculationRename={calculationOperations.handleCalculationRename}
            createNewCalculationFromExisting={
              activeCalculationData.handleCreateNewFromExisting
            }
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage
            activeCalculation={activeCalculationData.activeCalculation}
            isLoadingDetails={false}
            detailsError={null}
            onCalculationUpdate={
              activeCalculationData.handleActiveCalculationUpdate
            }
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculationData.activeCalculation}
            onCalculationUpdate={
              activeCalculationData.handleActiveCalculationUpdate
            }
            onStartCalculation={calculationOperations.handleStartCalculation}
            onCalculationRename={calculationOperations.handleCalculationRename}
            createNewCalculationFromExisting={
              activeCalculationData.handleCreateNewFromExisting
            }
          />
        );
    }
  };

  return (
    <div className="app-container">
      <button
        className={`independent-sidebar-toggle ${sidebarState.isSidebarOpen ? 'sidebar-open' : ''}`}
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
        onPlusClick={activeCalculationData.handleNewCalculation}
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
        calculations={activeCalculationData.sidebarCalculations}
        activeCalculationId={
          activeCalculationData.activeCalculation?.id || null
        }
        calculationsLoading={activeCalculationData.calculationsLoading}
        calculationsError={
          activeCalculationData.calculationsError
            ? activeCalculationData.calculationsError.message
            : null
        }
        onCalculationSelect={activeCalculationData.handleCalculationSelect}
        onCalculationDelete={calculationOperations.handleCalculationDelete}
      />

      <main
        className={`app-content ${sidebarState.isSidebarOpen ? 'sidebar-open' : ''}`}
      >
        {renderCurrentPage()}
      </main>
    </div>
  );
};
