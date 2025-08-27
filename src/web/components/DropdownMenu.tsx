import React from 'react';
import styles from './DropdownMenu.module.css';

export type DropdownOption =
  | 'calculation-settings'
  | 'calculation-results'
  | 'draw-molecule';

interface DropdownMenuProps {
  isOpen: boolean;
  selectedOption: DropdownOption;
  onOptionSelect: (option: DropdownOption) => void;
  onClose: () => void;
}

export const DropdownMenu: React.FC<DropdownMenuProps> = ({
  isOpen,
  selectedOption,
  onOptionSelect,
  onClose,
}) => {
  if (!isOpen) return null;

  const handleOptionClick = (option: DropdownOption) => {
    onOptionSelect(option);
    onClose();
  };

  return (
    <>
      {/* Backdrop to close dropdown when clicking outside */}
      <div className={styles.dropdownBackdrop} onClick={onClose} />

      <div className={styles.dropdownMenu}>
        <div
          className={`${styles.dropdownItem} ${selectedOption === 'calculation-settings' ? styles.selected : ''}`}
          onClick={() => handleOptionClick('calculation-settings')}
        >
          <svg
            width="16"
            height="16"
            viewBox="0 0 16 16"
            fill="none"
            xmlns="http://www.w3.org/2000/svg"
            className={styles.dropdownIcon}
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
        </div>

        <div
          className={`${styles.dropdownItem} ${selectedOption === 'calculation-results' ? styles.selected : ''}`}
          onClick={() => handleOptionClick('calculation-results')}
        >
          <svg
            width="16"
            height="16"
            viewBox="0 0 16 16"
            fill="none"
            xmlns="http://www.w3.org/2000/svg"
            className={styles.dropdownIcon}
          >
            <path
              d="M3 2C2.44772 2 2 2.44772 2 3V13C2 13.5523 2.44772 14 3 14H13C13.5523 14 14 13.5523 14 13V5.41421C14 5.01639 13.842 4.63486 13.5607 4.35355L10.6464 1.43934C10.3651 1.15803 9.98361 1 9.58579 1H3Z"
              stroke="currentColor"
              strokeWidth="1.5"
              fill="none"
            />
            <path
              d="M10 1V4C10 4.55228 10.4477 5 11 5H14"
              stroke="currentColor"
              strokeWidth="1.5"
              fill="none"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
          <span>Calculation results</span>
        </div>

        <div
          className={`${styles.dropdownItem} ${selectedOption === 'draw-molecule' ? styles.selected : ''}`}
          onClick={() => handleOptionClick('draw-molecule')}
        >
          <svg
            width="16"
            height="16"
            viewBox="0 0 16 16"
            fill="none"
            xmlns="http://www.w3.org/2000/svg"
            className={styles.dropdownIcon}
          >
            <path
              d="M12.854 1.854a.5.5 0 0 0-.708-.708L10.5 2.793 8.354.646a.5.5 0 1 0-.708.708L9.293 3 1.146 11.146a.5.5 0 0 0-.146.354V14a.5.5 0 0 0 .5.5h2.5a.5.5 0 0 0 .354-.146L12.5 6.207l1.647 1.647a.5.5 0 0 0 .708-.708L12.707 5l1.647-1.647a.5.5 0 0 0 0-.708L12.854 1.854zM11.207 4L4 11.207V13h1.793L13 5.793 11.207 4z"
              fill="currentColor"
            />
          </svg>
          <span>Draw molecule</span>
        </div>
      </div>
    </>
  );
};
