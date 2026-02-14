import React from 'react';
import styles from '../Sidebar.module.css';

interface SidebarFooterProps {
  onUserMenuToggle: () => void;
  isUserMenuOpen: boolean;
  onSettingsOpen: () => void;
}

export const SidebarFooter: React.FC<SidebarFooterProps> = ({
  onUserMenuToggle,
  isUserMenuOpen,
  onSettingsOpen,
}) => {
  return (
    <div className={styles.sidebarBottomSection}>
      <div
        className={`${styles.userInfoSection} ${isUserMenuOpen ? styles.userMenuOpen : ''}`}
        onClick={onUserMenuToggle}
      >
        {isUserMenuOpen && (
          <div className={styles.userMenu}>
            <button className={styles.userMenuItem} onClick={onSettingsOpen}>
              Settings
            </button>
            <button
              className={styles.userMenuItem}
              onClick={() => {
                window.electronAPI.showAboutDialog();
              }}
            >
              About
            </button>
          </div>
        )}
        <div className={styles.userInfoMain}>
          <span className={styles.menuLabel}>menu</span>
          <span
            className={`${styles.userMenuToggle} ${isUserMenuOpen ? styles.rotated : ''}`}
          >
            <svg
              width="24"
              height="24"
              viewBox="0 0 16 16"
              fill="none"
              xmlns="http://www.w3.org/2000/svg"
              style={{ transformOrigin: 'center center' }}
            >
              <path
                d="M4 6L8 10L12 6"
                stroke="currentColor"
                strokeWidth="1.5"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </span>
        </div>
      </div>
    </div>
  );
};
