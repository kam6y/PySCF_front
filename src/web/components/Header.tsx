import React, { useEffect, useState } from "react";

interface HeaderProps {
  onSidebarToggle: () => void;
  onDropdownToggle: () => void;
  onPlusClick: () => void;
  isDropdownOpen: boolean;
  isSidebarOpen: boolean;
}

export const Header: React.FC<HeaderProps> = ({
  onSidebarToggle,
  onDropdownToggle,
  onPlusClick,
  isDropdownOpen,
  isSidebarOpen,
}) => {
  const [isMaximized, setIsMaximized] = useState(false);

  useEffect(() => {
    // Check if electronAPI is available (we're in Electron)
    if (typeof window !== 'undefined' && window.electronAPI) {
      // Get initial window state
      window.electronAPI.isWindowMaximized().then(setIsMaximized);

      // Listen for window state changes
      window.electronAPI.onWindowStateChange(setIsMaximized);

      // Cleanup
      return () => {
        window.electronAPI.removeAllListeners('window-state-changed');
      };
    }
  }, []);

  const handleCloseWindow = () => {};
  return (
    <header className="app-header-new">
      {/* Left Section */}
      <div className="header-left">
        <button 
          className="sidebar-toggle"
          onClick={onSidebarToggle}
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
                strokeWidth="1.5"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            ) : (
              // Arrow pointing right (open sidebar)
              <path
                d="M6 4L10 8L6 12"
                stroke="currentColor"
                strokeWidth="1.5"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            )}
          </svg>
        </button>
        <h1 className="app-title">PySCF_app</h1>
      </div>

      {/* Right Section */}
      <div className="header-right">
        <button 
          className={`dropdown-button ${isDropdownOpen ? 'active' : ''}`}
          onClick={onDropdownToggle}
        >
          <svg
            width="16"
            height="16"
            viewBox="0 0 16 16"
            fill="none"
            xmlns="http://www.w3.org/2000/svg"
          >
            <path
              d="M8 1.5C8 1.22386 7.77614 1 7.5 1H6.5C6.22386 1 6 1.22386 6 1.5V2.5C6 2.77614 6.22386 3 6.5 3H7.5C7.77614 3 8 2.77614 8 2.5V1.5Z"
              fill="currentColor"
            />
            <path
              d="M10.5 6C10.7761 6 11 6.22386 11 6.5V7.5C11 7.77614 10.7761 8 10.5 8H9.5C9.22386 8 9 7.77614 9 7.5V6.5C9 6.22386 9.22386 6 9.5 6H10.5Z"
              fill="currentColor"
            />
            <path
              d="M6 9.5C6 9.22386 6.22386 9 6.5 9H7.5C7.77614 9 8 9.22386 8 9.5V10.5C8 10.77614 7.77614 11 7.5 11H6.5C6.22386 11 6 10.77614 6 10.5V9.5Z"
              fill="currentColor"
            />
            <path
              d="M14.5 8C14.7761 8 15 8.22386 15 8.5V9.5C15 9.77614 14.7761 10 14.5 10H13.5C13.22386 10 13 9.77614 13 9.5V8.5C13 8.22386 13.22386 8 13.5 8H14.5Z"
              fill="currentColor"
            />
            <path
              d="M2 6.5C2 6.22386 2.22386 6 2.5 6H3.5C3.77614 6 4 6.22386 4 6.5V7.5C4 7.77614 3.77614 8 3.5 8H2.5C2.22386 8 2 7.77614 2 7.5V6.5Z"
              fill="currentColor"
            />
            <path
              d="M10.5 2C10.7761 2 11 2.22386 11 2.5V3.5C11 3.77614 10.7761 4 10.5 4H9.5C9.22386 4 9 3.77614 9 3.5V2.5C9 2.22386 9.22386 2 9.5 2H10.5Z"
              fill="currentColor"
            />
            <path
              d="M14 4.5C14 4.22386 13.7761 4 13.5 4H12.5C12.22386 4 12 4.22386 12 4.5V5.5C12 5.77614 12.22386 6 12.5 6H13.5C13.7761 6 14 5.77614 14 5.5V4.5Z"
              fill="currentColor"
            />
            <path
              d="M2.5 10C2.77614 10 3 10.2239 3 10.5V11.5C3 11.7761 2.77614 12 2.5 12H1.5C1.22386 12 1 11.7761 1 11.5V10.5C1 10.2239 1.22386 10 1.5 10H2.5Z"
              fill="currentColor"
            />
          </svg>
          <span>Calculation settings</span>
          <svg
            width="12"
            height="12"
            viewBox="0 0 12 12"
            fill="none"
            xmlns="http://www.w3.org/2000/svg"
            className={`dropdown-arrow ${isDropdownOpen ? 'open' : ''}`}
          >
            <path
              d="M3 4.5L6 7.5L9 4.5"
              stroke="currentColor"
              strokeWidth="1.5"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </button>

        <button 
          className="plus-button"
          onClick={onPlusClick}
          aria-label="Add new calculation"
        >
          <svg
            width="16"
            height="16"
            viewBox="0 0 16 16"
            fill="none"
            xmlns="http://www.w3.org/2000/svg"
          >
            <path
              d="M8 3V13M3 8H13"
              stroke="currentColor"
              strokeWidth="2"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </button>
      </div>
    </header>
  );
};