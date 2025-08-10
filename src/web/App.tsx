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
  // --- State Management ---
  const [isSidebarOpen, setIsSidebarOpen] = useState(true);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const [currentPage, setCurrentPage] = useState<DropdownOption>('calculation-settings');

  const {
    calculations,
    isLoading: calculationsLoading,
    error: calculationsError,
    refreshCalculations,
    updateCalculationName,
    deleteCalculation,
  } = useCalculations();

  const {
    activeCalculationId,
    activeCalculation,
    isLoadingDetails,
    detailsError,
    setActiveCalculation,
    setActiveCalculationById,
    loadCalculationDetails,
    clearActiveCalculation,
  } = useActiveCalculation(calculations);

  // --- Effect Hooks ---
  useEffect(() => {
    if (!activeCalculationId && calculations.length > 0) {
      const firstCalculation = calculations[0];
      setActiveCalculationById(firstCalculation.id);
    }
  }, [calculations, activeCalculationId, setActiveCalculationById]);

  // --- Helper Functions ---
  const getPageTitle = (page: DropdownOption): string => {
    const titles: Record<DropdownOption, string> = {
      'calculation-settings': 'Calculation Settings',
      'calculation-results': 'Calculation Results',
      'draw-molecule': 'Draw Molecule',
    };
    return titles[page] || 'PySCF Front';
  };

  // --- Event Handlers ---
  const handleSidebarToggle = () => setIsSidebarOpen(!isSidebarOpen);
  const handleSidebarClose = () => setIsSidebarOpen(false);
  const handleDropdownToggle = () => setIsDropdownOpen(!isDropdownOpen);
  const handleDropdownClose = () => setIsDropdownOpen(false);
  const handleDropdownOptionSelect = (option: DropdownOption) => setCurrentPage(option);

  /**
   * Handles selecting a calculation from the sidebar.
   * This function now ONLY sets the active calculation and loads its details.
   * It no longer changes the current page.
   */
  const handleCalculationSelect = (calculationId: string) => {
    setActiveCalculationById(calculationId);
    loadCalculationDetails(calculationId, true); // Force refresh details
    handleSidebarClose();
  };

  const handleNewCalculation = () => {
    const now = new Date().toISOString();
    const tempId = `new-calculation-${Date.now()}`;
    const newCalculation: CalculationInstance = {
      id: tempId,
      name: 'New Calculation',
      status: 'pending',
      createdAt: now,
      updatedAt: now,
      parameters: {
        calculation_method: 'DFT',
        basis_function: '6-31G(d)',
        exchange_correlation: 'B3LYP',
        charges: 0,
        spin_multiplicity: 1,
        solvent_method: 'none',
        solvent: '-',
        xyz: '',
        molecule_name: 'New Calculation'
      },
    };
    setActiveCalculation(newCalculation);
    setCurrentPage('calculation-settings');
    handleSidebarClose();
  };
  
  const handleCalculationSuccess = async (completedCalculation: CalculationInstance) => {
    await refreshCalculations(); 
    setActiveCalculation(completedCalculation);
    setCurrentPage('calculation-results');
  };

  const handleCalculationRename = async (calculationId: string, newName: string) => {
    try {
      await updateCalculationName(calculationId, newName);
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

  // --- Page Rendering Logic ---
  const renderCurrentPage = () => {
    switch (currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage
            activeCalculation={activeCalculation}
            onCalculationUpdate={setActiveCalculation}
            onCalculationSuccess={handleCalculationSuccess}
            refreshCalculations={refreshCalculations}
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
            onCalculationUpdate={setActiveCalculation}
            onCalculationSuccess={handleCalculationSuccess}
            refreshCalculations={refreshCalculations}
          />
        );
    }
  };

  // --- Main Component Render ---
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
        calculations={calculations}
        activeCalculationId={activeCalculationId}
        calculationsLoading={calculationsLoading}
        calculationsError={calculationsError}
        onCalculationSelect={handleCalculationSelect}
        onCalculationRename={handleCalculationRename}
        onCalculationDelete={handleCalculationDelete}
      />

      <main className={`app-content ${isSidebarOpen ? 'sidebar-open' : ''}`}>
        {renderCurrentPage()}
      </main>
    </div>
  );
};