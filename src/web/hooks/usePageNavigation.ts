import { useState } from 'react';
import { DropdownOption } from '../components/DropdownMenu';

export interface PageNavigation {
  currentPage: DropdownOption;
  currentPageTitle: string;
  handleDropdownOptionSelect: (option: DropdownOption) => void;
}

export const usePageNavigation = (): PageNavigation => {
  const [currentPage, setCurrentPage] = useState<DropdownOption>(
    'calculation-settings'
  );

  const getPageTitle = (page: DropdownOption): string => {
    const titles: Record<DropdownOption, string> = {
      'calculation-settings': 'Calculation Settings',
      'calculation-results': 'Calculation Results',
      'draw-molecule': 'Draw Molecule',
    };
    return titles[page] || 'PySCF_front';
  };

  const handleDropdownOptionSelect = (option: DropdownOption) =>
    setCurrentPage(option);

  return {
    currentPage,
    currentPageTitle: getPageTitle(currentPage),
    handleDropdownOptionSelect,
  };
};
