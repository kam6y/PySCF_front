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
  useActiveCalculationId,
  useCalculationData,
  useCalculationActions,
  useStagedCalculation,
} from './hooks';
import { CalculationInstance } from './types/api-types';

export const App = () => {
  const queryClient = useQueryClient();

  // 新しいリファクタリングされたフック（責務分離）
  const sidebarState = useSidebarState();
  const pageNavigation = usePageNavigation();
  const { activeCalculationId, selectCalculation } = useActiveCalculationId();
  const calculationData = useCalculationData(activeCalculationId);
  const calculationActions = useCalculationActions();
  const stagedCalculation = useStagedCalculation();

  // アクティブ計算の決定：staged計算を優先、なければサーバーから詳細データ、最後にリストから基本データ
  const activeCalculation = 
    stagedCalculation.stagedCalculation ||
    calculationData.activeCalculationDetail ||
    calculationData.activeCalculationBasic;

  // WebSocketによるリアルタイム更新
  useCalculationSubscription({
    calculationId: activeCalculation?.id || null,
    status: activeCalculation?.status,
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
            activeCalculation={activeCalculation}
            onCalculationUpdate={
              stagedCalculation.isStagedCalculation(activeCalculation?.id || null)
                ? stagedCalculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={stagedCalculation.createNewFromExisting}
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage
            activeCalculation={activeCalculation}
            isLoadingDetails={calculationData.detailsLoading}
            detailsError={null}
            onCalculationUpdate={calculationActions.handleCalculationUpdate}
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation}
            onCalculationUpdate={
              stagedCalculation.isStagedCalculation(activeCalculation?.id || null)
                ? stagedCalculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={stagedCalculation.createNewFromExisting}
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
        onPlusClick={() => stagedCalculation.createNewCalculation(
          pageNavigation.handleDropdownOptionSelect,
          sidebarState.handleSidebarClose
        )}
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
        calculations={calculationData.sidebarCalculations}
        activeCalculationId={activeCalculation?.id || null}
        calculationsLoading={calculationData.calculationsLoading}
        calculationsError={
          calculationData.calculationsError
            ? calculationData.calculationsError.message
            : null
        }
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
      />

      <main
        className={`app-content ${sidebarState.isSidebarOpen ? 'sidebar-open' : ''}`}
      >
        {renderCurrentPage()}
      </main>
    </div>
  );
};
