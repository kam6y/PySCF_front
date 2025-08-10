import { useState } from "react";
import "./App.css";
import { Header } from "./components/Header";
import { Sidebar } from "./components/Sidebar";
import { DropdownOption } from "./components/DropdownMenu";
import { CalculationSettingsPage } from "./pages/CalculationSettingsPage";
import { CalculationResultsPage } from "./pages/CalculationResultsPage";
import { DrawMoleculePage } from "./pages/DrawMoleculePage";
import { useCalculations, useActiveCalculation } from "./hooks";
import { CalculationParameters } from "./types/calculation";

export const App = () => {
  // UI state
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const [currentPage, setCurrentPage] = useState<DropdownOption>('calculation-settings');

  // Calculation management state
  const {
    calculations,
    isLoading: calculationsLoading,
    error: calculationsError,
    refreshCalculations,
    getCalculationById,
    updateCalculationName,
    deleteCalculation,
    createCalculation
  } = useCalculations();

  const {
    activeCalculationId,
    activeCalculation,
    setActiveCalculation,
    setActiveCalculationById,
    loadCalculationDetails,
    clearActiveCalculation
  } = useActiveCalculation(calculations);

  // Get page title based on current page
  const getPageTitle = (page: DropdownOption): string => {
    switch (page) {
      case 'calculation-settings':
        return 'Calculation settings';
      case 'calculation-results':
        return 'Calculation results';
      case 'draw-molecule':
        return 'Draw molecule';
      default:
        return 'Calculation settings';
    }
  };

  // Event handlers for UI components
  const handleSidebarToggle = () => {
    setIsSidebarOpen(!isSidebarOpen);
  };

  const handleSidebarClose = () => {
    setIsSidebarOpen(false);
  };

  const handleDropdownToggle = () => {
    setIsDropdownOpen(!isDropdownOpen);
  };

  const handleDropdownClose = () => {
    setIsDropdownOpen(false);
  };

  const handleDropdownOptionSelect = (option: DropdownOption) => {
    setCurrentPage(option);
  };

  // Handle calculation management actions
  const handleCalculationSelect = (calculationId: string) => {
    setActiveCalculationById(calculationId);
    setIsSidebarOpen(false); // Close sidebar when selecting calculation
  };

  const handleCalculationRename = async (calculationId: string, newName: string) => {
    try {
      await updateCalculationName(calculationId, newName);
    } catch (error) {
      console.error('Failed to rename calculation:', error);
    }
  };

  const handleCalculationDelete = async (calculationId: string) => {
    try {
      await deleteCalculation(calculationId);
      // If deleted calculation was active, clear it
      if (activeCalculationId === calculationId) {
        clearActiveCalculation();
      }
    } catch (error) {
      console.error('Failed to delete calculation:', error);
    }
  };

  const handlePlusClick = () => {
    // Create a new calculation and make it active
    const defaultParams: CalculationParameters = {
      calculation_method: 'DFT',
      basis_function: '6-31G(d)',
      exchange_correlation: 'B3LYP',
      charges: 0,
      spin_multiplicity: 1,
      solvent_method: 'none',
      solvent: '-',
      xyz: '',
      molecule_name: 'New Calculation'
    };

    const newCalculation = createCalculation('New Calculation', defaultParams);
    setActiveCalculation(newCalculation);
    setCurrentPage('calculation-settings'); // Navigate to settings page
  };

  // Render current page component
  const renderCurrentPage = () => {
    switch (currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage 
            activeCalculation={activeCalculation}
            onCalculationUpdate={(updatedCalculation) => {
              setActiveCalculation(updatedCalculation);
            }}
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage 
            activeCalculation={activeCalculation}
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage 
            activeCalculation={activeCalculation}
            onCalculationUpdate={(updatedCalculation) => {
              setActiveCalculation(updatedCalculation);
            }}
          />
        );
    }
  };

  return (
    <div className="app-container">
      {/* Independent Sidebar Toggle */}
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
            // Arrow pointing left (close sidebar)
            <path
              d="M10 4L6 8L10 12"
              stroke="currentColor"
              strokeWidth="3"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          ) : (
            // Arrow pointing right (open sidebar)
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

      {/* Header */}
      <Header
        onDropdownToggle={handleDropdownToggle}
        onPlusClick={handlePlusClick}
        isDropdownOpen={isDropdownOpen}
        isSidebarOpen={isSidebarOpen}
        currentPageTitle={getPageTitle(currentPage)}
        currentPage={currentPage}
        onDropdownOptionSelect={handleDropdownOptionSelect}
        onDropdownClose={handleDropdownClose}
      />


      {/* Sidebar */}
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

      {/* Main Content */}
      <div className={`app-content ${isSidebarOpen ? 'sidebar-open' : ''}`}>
        {renderCurrentPage()}
      </div>

    </div>
  );
};