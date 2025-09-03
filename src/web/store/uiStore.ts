import { create } from 'zustand';
import { DropdownOption } from '../components/DropdownMenu';

interface UIState {
  // Sidebar状態
  isSidebarOpen: boolean;
  isDropdownOpen: boolean;
  
  // Page navigation状態  
  currentPage: DropdownOption;
  
  // Search状態
  searchQuery: string;
  
  // Modal状態
  isUserMenuOpen: boolean;
  isSettingsOpen: boolean;
  
  // Sidebar actions
  toggleSidebar: () => void;
  closeSidebar: () => void;
  toggleDropdown: () => void;
  closeDropdown: () => void;
  
  // Page navigation actions
  setCurrentPage: (page: DropdownOption) => void;
  getCurrentPageTitle: () => string;
  
  // Search actions
  setSearchQuery: (query: string) => void;
  
  // Modal actions
  toggleUserMenu: () => void;
  closeUserMenu: () => void;
  openSettings: () => void;
  closeSettings: () => void;
}

const getPageTitle = (page: DropdownOption): string => {
  const titles: Record<DropdownOption, string> = {
    'calculation-settings': 'Calculation Settings',
    'calculation-results': 'Calculation Results',
    'draw-molecule': 'Draw Molecule',
  };
  return titles[page] || 'PySCF_front';
};

export const useUIStore = create<UIState>((set, get) => ({
  // Initial state
  isSidebarOpen: false,
  isDropdownOpen: false,
  currentPage: 'calculation-settings',
  searchQuery: '',
  isUserMenuOpen: false,
  isSettingsOpen: false,

  // Sidebar actions
  toggleSidebar: () => set((state) => ({ isSidebarOpen: !state.isSidebarOpen })),
  closeSidebar: () => set({ isSidebarOpen: false }),
  toggleDropdown: () => set((state) => ({ isDropdownOpen: !state.isDropdownOpen })),
  closeDropdown: () => set({ isDropdownOpen: false }),

  // Page navigation actions
  setCurrentPage: (page: DropdownOption) => set({ currentPage: page }),
  getCurrentPageTitle: () => getPageTitle(get().currentPage),

  // Search actions  
  setSearchQuery: (query: string) => set({ searchQuery: query }),

  // Modal actions
  toggleUserMenu: () => set((state) => ({ isUserMenuOpen: !state.isUserMenuOpen })),
  closeUserMenu: () => set({ isUserMenuOpen: false }),
  openSettings: () => set({ 
    isSettingsOpen: true, 
    isUserMenuOpen: false // Close user menu when opening settings
  }),
  closeSettings: () => set({ isSettingsOpen: false }),
}));