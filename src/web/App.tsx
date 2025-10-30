// src/web/App.tsx

import React, { useMemo, useState, useEffect } from 'react';
import './App.css';
import styles from './App.module.css';
import { Header } from './components/Header';
import { Sidebar } from './components/Sidebar';
import { ToastContainer } from './components/ToastContainer';
import { InitialSetupDialog } from './components/InitialSetupDialog';
import { CalculationSettingsPage } from './pages/CalculationSettingsPage';
import { CalculationResultsPage } from './pages/CalculationResultsPage';
import { SettingsPage } from './pages/SettingsPage';
import { AgentPage } from './pages/AgentPage';
import { DrawMoleculePage } from './pages/DrawMoleculePage';
import { useAppState } from './hooks/useAppState';
import { useCalculationData } from './hooks/useCalculationData';
import { useCalculationActions } from './hooks/useCalculationActions';
import { useUnifiedWebSocket } from './hooks/useUnifiedWebSocket';
import { useChatHistoryStore } from './store/chatHistoryStore';
import { useAgentStore } from './store/agentStore';
import {
  useGetChatSessionDetail,
  useGetChatSessions,
} from './hooks/useChatHistoryQueries';
import { useAppSettings } from './hooks/useAppSettings';

export const App = () => {
  // 統合された状態管理
  const appState = useAppState();
  const calculationData = useCalculationData();
  const calculationActions = useCalculationActions();

  // アプリ設定の取得
  const { settings, updateSettings } = useAppSettings();

  // 初回セットアップダイアログの表示状態
  const [showInitialSetup, setShowInitialSetup] = useState(false);
  const [defaultDirectory, setDefaultDirectory] = useState('');

  // プラットフォーム情報の取得
  const [platform, setPlatform] = useState<string>('');
  const [isFullScreen, setIsFullScreen] = useState<boolean>(false);

  useEffect(() => {
    const fetchPlatform = async () => {
      try {
        const platformName = await window.electronAPI.getPlatform();
        setPlatform(platformName);
      } catch (error) {
        console.error('Failed to get platform:', error);
      }
    };
    fetchPlatform();
  }, []);

  // 全画面状態の取得と監視
  useEffect(() => {
    const fetchFullScreen = async () => {
      try {
        const fullScreen = await window.electronAPI.isFullScreen();
        setIsFullScreen(fullScreen);
      } catch (error) {
        console.error('Failed to get fullscreen state:', error);
      }
    };
    fetchFullScreen();

    // 全画面状態の変更を監視
    const cleanup = window.electronAPI.onFullScreenChange(fullScreen => {
      setIsFullScreen(fullScreen);
    });

    return cleanup;
  }, []);

  // 初回起動チェック
  useEffect(() => {
    if (settings) {
      // LocalStorageでセットアップ完了フラグを確認
      const setupCompleted = localStorage.getItem('pyscf_setup_completed');

      if (!setupCompleted) {
        // 初回起動の場合、デフォルトディレクトリを設定してダイアログを表示
        setDefaultDirectory(settings.calculations_directory);
        setShowInitialSetup(true);
        console.log('[App] First run detected, showing initial setup dialog');
      } else {
        setShowInitialSetup(false);
        console.log(
          '[App] Setup already completed, skipping initial setup dialog'
        );
      }
    }
  }, [settings]);

  // 初回セットアップ完了時の処理
  const handleSetupComplete = async (calculationsDirectory: string) => {
    if (settings) {
      await updateSettings({
        ...settings,
        calculations_directory: calculationsDirectory,
      });
    }

    // セットアップ完了フラグをLocalStorageに保存
    localStorage.setItem('pyscf_setup_completed', 'true');
    console.log('[App] Initial setup completed, saved flag to localStorage');

    setShowInitialSetup(false);

    // 初回セットアップ完了後、自動的に新規計算を作成
    console.log('[App] Auto-creating new calculation after setup completion');
    appState.calculation.createNewCalculation();
  };

  // チャットセッションのデータ取得
  const { data: chatSessionsData } = useGetChatSessions();

  // アプリ起動時に自動的に新規計算を作成
  useEffect(() => {
    // 初回セットアップダイアログが閉じており、settingsがロードされ、
    // まだ新規計算が作成されていない場合のみ実行
    if (
      !showInitialSetup &&
      settings &&
      !appState.calculation.stagedCalculation
    ) {
      console.log('[App] Auto-creating new calculation on app startup');
      appState.calculation.createNewCalculation();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [showInitialSetup, settings]);

  // AIチャット画面を開いた時にサイドバーのタブを"chats"に自動切り替え
  // AI Agent ページがOFFの場合は"calculations"に戻す
  useEffect(() => {
    if (appState.ui.isAIAgentEnabled) {
      appState.ui.setSidebarView('chats');
    } else {
      appState.ui.setSidebarView('calculations');
    }
  }, [appState.ui.isAIAgentEnabled, appState.ui.setSidebarView]);

  // 統合WebSocketによるリアルタイム更新（グローバル + アクティブ計算監視）
  useUnifiedWebSocket({
    activeCalculationId: calculationData.activeCalculation?.id || null,
  });

  // 検索機能によるフィルタリング（計算インスタンス）
  const filteredCalculations = useMemo(() => {
    const searchQuery = appState.ui.searchQuery;
    if (!searchQuery.trim()) {
      return calculationData.sidebarCalculations;
    }

    const query = searchQuery.toLowerCase();
    return calculationData.sidebarCalculations.filter(
      calc =>
        calc.name.toLowerCase().includes(query) ||
        calc.status.toLowerCase().includes(query)
    );
  }, [calculationData.sidebarCalculations, appState.ui.searchQuery]);

  // 検索機能によるフィルタリング（チャットセッション）
  const filteredChatSessions = useMemo(() => {
    const searchQuery = appState.ui.searchQuery;
    const sessions = chatSessionsData?.sessions || [];

    if (!searchQuery.trim()) {
      return sessions;
    }

    const query = searchQuery.toLowerCase();
    return sessions.filter(session =>
      session.name.toLowerCase().includes(query)
    );
  }, [chatSessionsData?.sessions, appState.ui.searchQuery]);

  // イベントハンドラー（統合されたアクションを使用）
  const handleSearchChange = (query: string) => {
    appState.ui.setSearchQuery(query);
  };

  const handleCalculationSelect = (calculationId: string) => {
    appState.actions.handleCalculationSelect(calculationId);
  };

  // チャット履歴のハンドラー
  const setActiveSessionId = useChatHistoryStore(
    state => state.setActiveSessionId
  );
  const setHistory = useAgentStore(state => state.setHistory);

  const handleChatSessionSelect = async (sessionId: string) => {
    try {
      // セッションIDを設定
      setActiveSessionId(sessionId);

      // AI Agentモードに切り替え
      appState.ui.setAIAgentEnabled(true);

      // セッション詳細を取得してagentStoreのhistoryを更新
      // （実際の読み込みはAgentPageで行う）
    } catch (error) {
      console.error('Failed to load chat session:', error);
    }
  };

  const renderCurrentPage = () => {
    // AI Agentが有効な場合はAgentPageを表示
    if (appState.ui.isAIAgentEnabled) {
      return <AgentPage />;
    }

    switch (appState.ui.currentPage) {
      case 'calculation-settings':
        return (
          <CalculationSettingsPage
            activeCalculation={calculationData.activeCalculation || undefined}
            onCalculationUpdate={
              calculationData.isStagedCalculation
                ? appState.calculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={
              appState.calculation.createNewFromExisting
            }
          />
        );
      case 'calculation-results':
        return (
          <CalculationResultsPage
            activeCalculation={calculationData.activeCalculation || undefined}
            isLoadingDetails={calculationData.detailsLoading}
            detailsError={null}
            onCalculationUpdate={calculationActions.handleCalculationUpdate}
          />
        );
      case 'draw-molecule':
        return <DrawMoleculePage />;
      default:
        return (
          <CalculationSettingsPage
            activeCalculation={calculationData.activeCalculation || undefined}
            onCalculationUpdate={
              calculationData.isStagedCalculation
                ? appState.calculation.updateStagedCalculation
                : calculationActions.handleCalculationUpdate
            }
            onStartCalculation={calculationActions.handleStartCalculation}
            onCalculationRename={calculationActions.handleCalculationRename}
            createNewCalculationFromExisting={
              appState.calculation.createNewFromExisting
            }
          />
        );
    }
  };

  return (
    <div className={styles.appContainer}>
      <button
        className={`${styles.independentSidebarToggle} ${appState.ui.isSidebarOpen ? styles.sidebarOpen : ''} ${
          platform === 'linux'
            ? styles.platformLinux
            : platform === 'darwin' && isFullScreen
              ? styles.platformMacFullScreen
              : styles.platformMac
        }`}
        onClick={appState.ui.toggleSidebar}
        aria-label="Toggle sidebar"
      >
        <svg
          width="16"
          height="16"
          viewBox="0 0 16 16"
          fill="none"
          xmlns="http://www.w3.org/2000/svg"
        >
          {appState.ui.isSidebarOpen ? (
            <path
              d="M10 4L6 8L10 12"
              stroke="currentColor"
              strokeWidth="3"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          ) : (
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

      <Header
        onDropdownToggle={appState.actions.handleDropdownToggle}
        onPlusClick={appState.actions.handleCreateNew}
        isDropdownOpen={appState.ui.isDropdownOpen}
        isSidebarOpen={appState.ui.isSidebarOpen}
        currentPageTitle={appState.ui.currentPageTitle}
        currentPage={appState.ui.currentPage}
        onDropdownOptionSelect={appState.ui.setCurrentPage}
        onDropdownClose={appState.ui.closeDropdown}
        isAIAgentEnabled={appState.ui.isAIAgentEnabled}
        onAIAgentToggle={appState.ui.setAIAgentEnabled}
      />

      <Sidebar
        isOpen={appState.ui.isSidebarOpen}
        onClose={appState.ui.closeSidebar}
        calculations={filteredCalculations}
        activeCalculationId={calculationData.activeCalculation?.id || null}
        calculationsLoading={calculationData.calculationsLoading}
        calculationsError={
          calculationData.calculationsError
            ? calculationData.calculationsError.message
            : null
        }
        onCalculationSelect={handleCalculationSelect}
        onCalculationDelete={async (calculationId: string) => {
          await calculationActions.handleCalculationDelete(calculationId);
          // 削除された計算がアクティブだった場合はクリア
          if (calculationData.activeCalculationId === calculationId) {
            appState.calculation.selectCalculation(null);
            appState.calculation.clearStaged();
            appState.ui.setCurrentPage('calculation-settings');
          }
        }}
        onCreateNew={appState.actions.handleCreateNew}
        searchQuery={appState.ui.searchQuery}
        onSearchChange={handleSearchChange}
        onUserMenuToggle={appState.ui.toggleUserMenu}
        isUserMenuOpen={appState.ui.isUserMenuOpen}
        onSettingsOpen={appState.ui.openSettings}
        sidebarView={appState.ui.sidebarView}
        onSidebarViewChange={appState.ui.setSidebarView}
        onChatSessionSelect={handleChatSessionSelect}
        filteredChatSessions={filteredChatSessions}
      />

      <main
        className={`${styles.appContent} ${appState.ui.isSidebarOpen ? styles.sidebarOpen : ''}`}
      >
        {renderCurrentPage()}
      </main>

      {/* Settings Overlay */}
      {appState.ui.isSettingsOpen && (
        <div className={styles.settingsOverlay}>
          <div
            className={styles.settingsModalBackdrop}
            onClick={appState.ui.closeSettings}
          />
          <div className={styles.settingsModal}>
            <div className={styles.settingsModalHeader}>
              <button
                onClick={appState.ui.closeSettings}
                className={styles.settingsCloseButton}
                aria-label="Close settings"
              >
                <svg
                  width="20"
                  height="20"
                  viewBox="0 0 20 20"
                  fill="none"
                  xmlns="http://www.w3.org/2000/svg"
                >
                  <path
                    d="M15 5L5 15M5 5L15 15"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </svg>
              </button>
            </div>
            <SettingsPage />
          </div>
        </div>
      )}

      {/* Initial Setup Dialog */}
      {showInitialSetup && settings && (
        <InitialSetupDialog
          onComplete={handleSetupComplete}
          defaultDirectory={settings.calculations_directory || defaultDirectory}
        />
      )}

      <ToastContainer />
    </div>
  );
};
