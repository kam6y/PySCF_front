import React from 'react';
import { CalculationSummary, ChatSessionSummary } from '../types/api-types';
import { SidebarView } from '../store/uiStore';
import { ChatHistoryList } from './ChatHistoryList';
import { CalculationsList } from './sidebar/CalculationsList';
import { DeleteModals } from './sidebar/DeleteModals';
import { SidebarFooter } from './sidebar/SidebarFooter';
import { useDeleteModals } from './sidebar/useDeleteModals';
import styles from './Sidebar.module.css';

interface SidebarProps {
  isOpen: boolean;
  onClose: () => void;
  calculations: CalculationSummary[];
  activeCalculationId: string | null;
  calculationsLoading: boolean;
  calculationsError: string | null;
  onCalculationSelect: (calculationId: string) => void;
  onCalculationDelete: (calculationId: string) => Promise<void>;
  onCreateNew: () => void;
  searchQuery: string;
  onSearchChange: (query: string) => void;
  onUserMenuToggle: () => void;
  isUserMenuOpen: boolean;
  onSettingsOpen: () => void;
  sidebarView: SidebarView;
  onSidebarViewChange: (view: SidebarView) => void;
  onChatSessionSelect: (sessionId: string) => void;
  filteredChatSessions: ChatSessionSummary[];
}

export const Sidebar: React.FC<SidebarProps> = ({
  isOpen,
  onClose,
  calculations,
  activeCalculationId,
  calculationsLoading,
  calculationsError,
  onCalculationSelect,
  onCalculationDelete,
  onCreateNew,
  searchQuery,
  onSearchChange,
  onUserMenuToggle,
  isUserMenuOpen,
  onSettingsOpen,
  sidebarView,
  onSidebarViewChange,
  onChatSessionSelect,
  filteredChatSessions,
}) => {
  const deleteModals = useDeleteModals({ onCalculationDelete });

  return (
    <>
      {isOpen && (
        <div
          className={styles.sidebarBackdrop}
          onClick={onClose}
          aria-label="Close sidebar"
        />
      )}

      <aside className={`${styles.sidebar} ${isOpen ? styles.open : ''}`}>
        <div className={styles.sidebarContent}>
          <div className={styles.sidebarTopSection}>
            <button className={styles.createNewButton} onClick={onCreateNew}>
              + New Calculation
            </button>
            <div className={styles.searchContainer}>
              <svg
                className={styles.searchIcon}
                width="16"
                height="16"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              >
                <circle cx="11" cy="11" r="8"></circle>
                <path d="m21 21-4.35-4.35"></path>
              </svg>
              <input
                type="text"
                placeholder="Value"
                value={searchQuery}
                onChange={e => onSearchChange(e.target.value)}
                className={styles.searchInput}
              />
              <button
                className={styles.searchClearButton}
                onClick={e => {
                  e.preventDefault();
                  onSearchChange('');
                }}
                style={{ display: searchQuery ? 'block' : 'none' }}
              >
                ×
              </button>
            </div>
          </div>

          <div className={styles.sidebarMainContent}>
            <div className={styles.sidebarHeader}>
              <div className={styles.sidebarViewTabs}>
                <button
                  className={`${styles.sidebarViewTab} ${styles.calculationsTab} ${
                    sidebarView === 'calculations' ? styles.active : ''
                  }`}
                  onClick={() => onSidebarViewChange('calculations')}
                >
                  Calculations
                </button>
                <button
                  className={`${styles.sidebarViewTab} ${styles.chatsTab} ${
                    sidebarView === 'chats' ? styles.active : ''
                  }`}
                  onClick={() => onSidebarViewChange('chats')}
                >
                  Chats
                </button>
              </div>
            </div>

            <div className={styles.sidebarCalculations}>
              {sidebarView === 'calculations' ? (
                <CalculationsList
                  calculations={calculations}
                  activeCalculationId={activeCalculationId}
                  calculationsLoading={calculationsLoading}
                  calculationsError={calculationsError}
                  onCalculationSelect={onCalculationSelect}
                  onRequestDelete={deleteModals.handleRequestDelete}
                  onBulkDeleteError={deleteModals.handleBulkDeleteError}
                />
              ) : (
                <ChatHistoryList
                  onSessionSelect={onChatSessionSelect}
                  onRequestDelete={deleteModals.handleRequestChatDelete}
                  filteredSessions={filteredChatSessions}
                  searchQuery={searchQuery}
                />
              )}
            </div>
          </div>

          <SidebarFooter
            onUserMenuToggle={onUserMenuToggle}
            isUserMenuOpen={isUserMenuOpen}
            onSettingsOpen={onSettingsOpen}
          />
        </div>
      </aside>

      <DeleteModals {...deleteModals} />
    </>
  );
};
