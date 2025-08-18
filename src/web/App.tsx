// src/web/App.tsx

import { useState, useEffect } from 'react';
import { useQueryClient } from '@tanstack/react-query';
import './App.css';
import { Header } from './components/Header';
import { Sidebar } from './components/Sidebar';
import { DropdownOption } from './components/DropdownMenu';
import { CalculationSettingsPage } from './pages/CalculationSettingsPage';
import { CalculationResultsPage } from './pages/CalculationResultsPage';
import { DrawMoleculePage } from './pages/DrawMoleculePage';
import { useCalculationSubscription } from './hooks/useCalculationSubscription';
import { useCalculationStore } from './store/calculationStore';
import {
  useGetCalculations,
  useGetCalculationDetails,
  useDeleteCalculation,
  useUpdateCalculationName,
  useStartCalculation,
} from './hooks/useCalculationQueries';
import {
  CalculationInstance,
  QuantumCalculationRequest,
  ApiError,
} from './types/api-types';

export const App = () => {
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const [currentPage, setCurrentPage] = useState<DropdownOption>(
    'calculation-settings'
  );

  const queryClient = useQueryClient();

  // Zustandストアからの状態
  const {
    activeCalculationId,
    stagedCalculation,
    setActiveCalculationId,
    setStagedCalculation,
    clearStagedCalculation,
  } = useCalculationStore();

  // TanStack Queryフック
  const {
    data: calculationsData,
    isLoading: calculationsLoading,
    error: calculationsError,
  } = useGetCalculations();
  const { data: activeCalculationDetails } =
    useGetCalculationDetails(activeCalculationId);
  const deleteCalculationMutation = useDeleteCalculation();
  const updateCalculationNameMutation = useUpdateCalculationName();
  const startCalculationMutation = useStartCalculation();

  // 計算リストのデータ変換
  const calculations =
    calculationsData?.calculations.map(c => ({
      id: c.id,
      name: c.name,
      status: c.status,
      createdAt: new Date(c.date).toISOString(),
      updatedAt: new Date(c.date).toISOString(),
      parameters: {} as any,
      results: undefined,
    })) || [];

  // アクティブな計算（詳細データがあればそれを使用、なければリストから取得）
  const activeCalculation =
    stagedCalculation ||
    activeCalculationDetails?.calculation ||
    calculations.find(c => c.id === activeCalculationId);

  // 最初の計算を自動選択（新規計算作成中は除く）
  useEffect(() => {
    if (!activeCalculationId && calculations.length > 0 && !stagedCalculation) {
      setActiveCalculationId(calculations[0].id);
    }
  }, [
    calculations,
    activeCalculationId,
    setActiveCalculationId,
    stagedCalculation,
  ]);

  // WebSocketによるリアルタイム更新
  useCalculationSubscription({
    calculationId: activeCalculationId,
    status: activeCalculation?.status,
    onUpdate: (updatedCalculation: CalculationInstance) => {
      // Queryキャッシュを直接更新
      queryClient.setQueryData(['calculation', updatedCalculation.id], {
        calculation: updatedCalculation,
      });
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    },
    onError: (error: string) => {
      console.error('WebSocket error:', error);
    },
  });

  const getPageTitle = (page: DropdownOption): string => {
    const titles: Record<DropdownOption, string> = {
      'calculation-settings': 'Calculation Settings',
      'calculation-results': 'Calculation Results',
      'draw-molecule': 'Draw Molecule',
    };
    return titles[page] || 'PySCF_front';
  };

  /**
   * エラーを適切に処理し、ユーザーフレンドリーなメッセージを表示するヘルパー関数
   */
  const handleApiError = (error: unknown, defaultMessage: string) => {
    console.error(defaultMessage, error);

    if (error instanceof ApiError) {
      // ApiErrorの場合、ユーザーフレンドリーなメッセージを表示
      alert(error.getUserMessage());
    } else if (error instanceof Error) {
      // 通常のErrorの場合、エラーメッセージを表示
      alert(`${defaultMessage}: ${error.message}`);
    } else {
      // 不明なエラーの場合、デフォルトメッセージを表示
      alert(`${defaultMessage}: 不明なエラーが発生しました。`);
    }
  };

  const handleSidebarToggle = () => setIsSidebarOpen(!isSidebarOpen);
  const handleSidebarClose = () => setIsSidebarOpen(false);
  const handleDropdownToggle = () => setIsDropdownOpen(!isDropdownOpen);
  const handleDropdownClose = () => setIsDropdownOpen(false);
  const handleDropdownOptionSelect = (option: DropdownOption) =>
    setCurrentPage(option);

  const handleCalculationSelect = (calculationId: string) => {
    setActiveCalculationId(calculationId);
    clearStagedCalculation();
    handleSidebarClose();
  };

  const handleNewCalculation = () => {
    const defaultParams: QuantumCalculationRequest = {
      calculation_method: 'DFT',
      basis_function: '6-31G(d)',
      exchange_correlation: 'B3LYP',
      charges: 0,
      spin_multiplicity: 1,
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

    setStagedCalculation(newCalculation);
    setActiveCalculationId(newId);
    setCurrentPage('calculation-settings');
    handleSidebarClose();
  };

  const handleStartCalculation = async (params: QuantumCalculationRequest) => {
    try {
      const response = await startCalculationMutation.mutateAsync(params);
      const runningCalculation = response.calculation;

      clearStagedCalculation();
      setActiveCalculationId(runningCalculation.id);

      return runningCalculation;
    } catch (error) {
      handleApiError(error, '計算の開始に失敗しました');
      throw error;
    }
  };

  const handleCalculationRename = async (
    calculationId: string,
    newName: string
  ) => {
    try {
      await updateCalculationNameMutation.mutateAsync({
        id: calculationId,
        newName,
      });
    } catch (error) {
      handleApiError(error, '計算名の変更に失敗しました');
    }
  };

  const handleCalculationDelete = async (calculationId: string) => {
    try {
      await deleteCalculationMutation.mutateAsync(calculationId);
      if (activeCalculationId === calculationId) {
        setActiveCalculationId(null);
        clearStagedCalculation();
        setCurrentPage('calculation-settings');
      }
    } catch (error) {
      handleApiError(error, '計算の削除に失敗しました');
    }
  };

  const handleActiveCalculationUpdate = (
    updatedCalculation: CalculationInstance
  ) => {
    if (stagedCalculation && updatedCalculation.id === stagedCalculation.id) {
      setStagedCalculation(updatedCalculation);
    } else {
      // Queryキャッシュを更新
      queryClient.setQueryData(['calculation', updatedCalculation.id], {
        calculation: updatedCalculation,
      });
      queryClient.invalidateQueries({ queryKey: ['calculations'] });
    }
  };

  const handleCreateNewFromExisting = (
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

    setStagedCalculation(newCalculation);
    setActiveCalculationId(newId);
  };

  const renderCurrentPage = () => {
    switch (currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation}
            onCalculationUpdate={handleActiveCalculationUpdate}
            onStartCalculation={handleStartCalculation}
            onCalculationRename={handleCalculationRename}
            createNewCalculationFromExisting={handleCreateNewFromExisting}
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage
            activeCalculation={activeCalculation}
            isLoadingDetails={false}
            detailsError={null}
            onCalculationUpdate={handleActiveCalculationUpdate}
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation}
            onCalculationUpdate={handleActiveCalculationUpdate}
            onStartCalculation={handleStartCalculation}
            onCalculationRename={handleCalculationRename}
            createNewCalculationFromExisting={handleCreateNewFromExisting}
          />
        );
    }
  };

  const sidebarCalculations = calculations.filter(
    calc =>
      calc.status !== 'pending' ||
      (calc.parameters &&
        calc.parameters.xyz &&
        calc.parameters.xyz.trim() !== '')
  );

  return (
    <div className="app-container">
      <button
        className={`independent-sidebar-toggle ${isSidebarOpen ? 'sidebar-open' : ''}`}
        onClick={handleSidebarToggle}
        aria-label="Toggle sidebar"
      >
        <svg
          width="16"
          height="16"
          viewBox="0 0 16 16"
          fill="none"
          xmlns="http://www.w3.org/2000/svg"
        >
          {isSidebarOpen ? (
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
        onDropdownToggle={handleDropdownToggle}
        onPlusClick={handleNewCalculation}
        isDropdownOpen={isDropdownOpen}
        isSidebarOpen={isSidebarOpen}
        currentPageTitle={getPageTitle(currentPage)}
        currentPage={currentPage}
        onDropdownOptionSelect={handleDropdownOptionSelect}
        onDropdownClose={handleDropdownClose}
      />

      <Sidebar
        isOpen={isSidebarOpen}
        onClose={handleSidebarClose}
        calculations={sidebarCalculations}
        activeCalculationId={activeCalculationId}
        calculationsLoading={calculationsLoading}
        calculationsError={calculationsError ? calculationsError.message : null}
        onCalculationSelect={handleCalculationSelect}
        onCalculationDelete={handleCalculationDelete}
      />

      <main className={`app-content ${isSidebarOpen ? 'sidebar-open' : ''}`}>
        {renderCurrentPage()}
      </main>
    </div>
  );
};
