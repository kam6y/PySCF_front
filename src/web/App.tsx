// src/web/App.tsx

import { useState, useEffect } from "react";
import "./App.css";
import { Header } from "./components/Header";
import { Sidebar } from "./components/Sidebar";
import { DropdownOption } from "./components/DropdownMenu";
import { CalculationSettingsPage } from "./pages/CalculationSettingsPage";
import { CalculationResultsPage } from "./pages/CalculationResultsPage";
import { DrawMoleculePage } from "./pages/DrawMoleculePage";
import { useCalculations, useActiveCalculation, useCalculationSubscription } from "./hooks";
import { CalculationInstance, QuantumCalculationRequest } from "./types/api-types";
import { startCalculation } from "./apiClient"; // apiClientからstartCalculationを直接インポート

export const App = () => {
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const [currentPage, setCurrentPage] = useState<DropdownOption>('calculation-settings');

  const {
    calculations,
    isLoading: calculationsLoading,
    error: calculationsError,
    refreshCalculations,
    updateCalculationName,
    deleteCalculation,
    createCalculation,
    updateCalculation,
    createNewCalculationFromExisting,
  } = useCalculations();

  const {
    activeCalculationId,
    activeCalculation,
    isLoadingDetails,
    detailsError,
    setActiveCalculationById,
    loadCalculationDetails,
    clearActiveCalculation,
    updateActiveCalculationInCache,
  } = useActiveCalculation(calculations);

  useEffect(() => {
    if (!activeCalculationId && calculations.length > 0) {
      setActiveCalculationById(calculations[0].id);
    }
  }, [calculations, activeCalculationId, setActiveCalculationById]);

  const getPageTitle = (page: DropdownOption): string => {
    const titles: Record<DropdownOption, string> = {
      'calculation-settings': 'Calculation Settings',
      'calculation-results': 'Calculation Results',
      'draw-molecule': 'Draw Molecule',
    };
    return titles[page] || 'PySCF_front';
  };

  const handleSidebarToggle = () => setIsSidebarOpen(!isSidebarOpen);
  const handleSidebarClose = () => setIsSidebarOpen(false);
  const handleDropdownToggle = () => setIsDropdownOpen(!isDropdownOpen);
  const handleDropdownClose = () => setIsDropdownOpen(false);
  const handleDropdownOptionSelect = (option: DropdownOption) => setCurrentPage(option);

  const handleCalculationSelect = (calculationId: string) => {
    setActiveCalculationById(calculationId);
    loadCalculationDetails(calculationId, true);
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
      name: '', // no default name (must be set by user)
      tddft_nstates: 10,
      tddft_method: 'TDDFT',
      tddft_analyze_nto: false
    };
    const newCalculation = createCalculation('', defaultParams); // no default name (must be set by user)
    setActiveCalculationById(newCalculation.id);
    setCurrentPage('calculation-settings');
    handleSidebarClose();
  };
  
  const handleStartCalculation = async (params: QuantumCalculationRequest) => {
    try {
      // 1. APIを呼び出して計算を開始し、新しい永続的な計算インスタンスを取得
      const response = await startCalculation(params);
      const runningCalculation = response.calculation;

      // 2. calculationsリストの状態を更新 (一時IDを永続IDに置き換え)
      updateCalculation(runningCalculation);

      // 3. アクティブな計算IDを新しい永続IDに直接設定
      setActiveCalculationById(runningCalculation.id);

      return runningCalculation; // ポーリング開始のために返す
    } catch (error) {
      console.error("Failed to start calculation:", error);
      alert(`Error starting calculation: ${error instanceof Error ? error.message : 'Unknown error'}`);
      throw error;
    }
  };


  const handleCalculationRename = async (calculationId: string, newName: string) => {
    try {
      await updateCalculationName(calculationId, newName);
      // Refresh calculations to get updated display name
      await refreshCalculations();
    } catch (error) {
      console.error('Failed to rename calculation:', error);
      alert(`Error: Could not rename calculation. ${error instanceof Error ? error.message : ''}`);
    }
  };

  const handleCalculationDelete = async (calculationId: string) => {
    try {
      await deleteCalculation(calculationId);
      if (activeCalculationId === calculationId) {
        clearActiveCalculation();
        setCurrentPage('calculation-settings');
      }
    } catch (error) {
      console.error('Failed to delete calculation:', error);
      alert(`Error: Could not delete calculation. ${error instanceof Error ? error.message : ''}`);
    }
  };
  
  const handleActiveCalculationUpdate = (updatedCalculation: CalculationInstance) => {
    updateCalculation(updatedCalculation);
    updateActiveCalculationInCache(updatedCalculation);
  };

  useCalculationSubscription({
    calculationId: activeCalculationId,
    status: activeCalculation?.status,
    onUpdate: handleActiveCalculationUpdate,
    onError: (error: string) => {
      console.error('WebSocket error:', error);
    }
  });

  const handleCreateNewFromExisting = (originalCalc: CalculationInstance, newParams: QuantumCalculationRequest) => {
    const newInstance = createNewCalculationFromExisting(originalCalc, newParams);
    setActiveCalculationById(newInstance.id);
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
            isLoadingDetails={isLoadingDetails}
            detailsError={detailsError}
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
    (calc) => calc.status !== 'pending' || (calc.parameters && calc.parameters.xyz && calc.parameters.xyz.trim() !== '')
  );

  return (
    <div className="app-container">
      <button
        className={`independent-sidebar-toggle ${isSidebarOpen ? 'sidebar-open' : ''}`}
        onClick={handleSidebarToggle}
        aria-label="Toggle sidebar"
      >
        <svg width="16" height="16" viewBox="0 0 16 16" fill="none" xmlns="http://www.w3.org/2000/svg">
          {isSidebarOpen ? (
            <path d="M10 4L6 8L10 12" stroke="currentColor" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" />
          ) : (
            <path d="M6 4L10 8L6 12" stroke="currentColor" strokeWidth="3" strokeLinecap="round" strokeLinejoin="round" />
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
        calculationsError={calculationsError}
        onCalculationSelect={handleCalculationSelect}
        onCalculationDelete={handleCalculationDelete}
      />

      <main className={`app-content ${isSidebarOpen ? 'sidebar-open' : ''}`}>
        {renderCurrentPage()}
      </main>
    </div>
  );
};
