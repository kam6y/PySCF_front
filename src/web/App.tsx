import { useState } from "react";
import "./App.css";
import { Header } from "./components/Header";
import { Sidebar } from "./components/Sidebar";
import { DropdownOption } from "./components/DropdownMenu";
import { CalculationSettingsPage } from "./pages/CalculationSettingsPage";
import { CalculationResultsPage } from "./pages/CalculationResultsPage";
import { DrawMoleculePage } from "./pages/DrawMoleculePage";

export const App = () => {
  // UI state
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const [currentPage, setCurrentPage] = useState<DropdownOption>('calculation-settings');

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

  const handlePlusClick = () => {
    // TODO: Implement add new calculation functionality
    console.log('Add new calculation');
  };

  // Render current page component
  const renderCurrentPage = () => {
    switch (currentPage) {
      case 'calculation-settings':
        return <CalculationSettingsPage />;
      case 'calculation-results':
        return <CalculationResultsPage />;
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return <CalculationSettingsPage />;
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
      />

      {/* Main Content */}
      <div className={`app-content ${isSidebarOpen ? 'sidebar-open' : ''}`}>
        {renderCurrentPage()}
      </div>

      {/* Footer */}
      <footer className="app-footer">
        <p>Powered by 3Dmol.js • Built with React and TypeScript</p>
      </footer>
    </div>
  );
};