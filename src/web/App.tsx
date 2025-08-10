import { useState, useEffect } from "react";
import "./App.css";
import { Header } from "./components/Header";
import { Sidebar } from "./components/Sidebar";
import { DropdownOption } from "./components/DropdownMenu";
import { CalculationSettingsPage } from "./pages/CalculationSettingsPage";
import { CalculationResultsPage } from "./pages/CalculationResultsPage";
import { DrawMoleculePage } from "./pages/DrawMoleculePage";
import { useCalculations, useActiveCalculation } from "./hooks";
import { CalculationInstance, CalculationParameters } from "./types/calculation";

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
    return titles[page] || 'PySCF Front';
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
    const defaultParams: CalculationParameters = {
      calculation_method: 'DFT',
      basis_function: '6-31G(d)',
      exchange_correlation: 'B3LYP',
      charges: 0,
      spin_multiplicity: 1,
      solvent_method: 'none',
      solvent: '-',
      xyz: '',
      molecule_name: '' // no default name
    };
    const newCalculation = createCalculation('', defaultParams);
    setActiveCalculationById(newCalculation.id);
    setCurrentPage('calculation-settings');
    handleSidebarClose();
  };
  
  const handleCalculationSuccess = async (completedCalculation: CalculationInstance) => {
    updateCalculation(completedCalculation);
    await refreshCalculations(); 
    setActiveCalculationById(completedCalculation.id);
    setCurrentPage('calculation-results');
  };

  const handleCalculationRename = async (calculationId: string, newName: string) => {
    try {
      const response = await updateCalculationName(calculationId, newName);
      // If the active calculation was renamed, update its ID
      if (activeCalculationId === response.old_id) {
        setActiveCalculationById(response.new_id);
      }
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
  };

  const handleCreateNewFromExisting = (originalCalc: CalculationInstance, newParams: CalculationParameters) => {
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
            onCalculationSuccess={handleCalculationSuccess}
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
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation}
            onCalculationUpdate={handleActiveCalculationUpdate}
            onCalculationSuccess={handleCalculationSuccess}
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