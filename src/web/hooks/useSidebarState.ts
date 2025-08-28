import { useState } from 'react';

export interface SidebarState {
  isSidebarOpen: boolean;
  isDropdownOpen: boolean;
  handleSidebarToggle: () => void;
  handleSidebarClose: () => void;
  handleDropdownToggle: () => void;
  handleDropdownClose: () => void;
}

export const useSidebarState = (): SidebarState => {
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);

  const handleSidebarToggle = () => setIsSidebarOpen(!isSidebarOpen);
  const handleSidebarClose = () => setIsSidebarOpen(false);
  const handleDropdownToggle = () => setIsDropdownOpen(!isDropdownOpen);
  const handleDropdownClose = () => setIsDropdownOpen(false);

  return {
    isSidebarOpen,
    isDropdownOpen,
    handleSidebarToggle,
    handleSidebarClose,
    handleDropdownToggle,
    handleDropdownClose,
  };
};
