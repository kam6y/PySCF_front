export interface ElectronAPI {
  // URLパラメータから取得したFlaskポート番号（preloadで設定）
  flaskPort: number | null;
  // IPC経由で受信した認証トークンを取得（トークンが届くまで待機）
  getAuthToken: () => Promise<string | null>;
  openExternalUrl: (
    url: string
  ) => Promise<{ success: boolean; error?: string }>;
  showAboutDialog: () => Promise<void>;
  selectFolder: () => Promise<{
    canceled: boolean;
    filePath: string | null;
    error?: string;
  }>;
  getPlatform: () => Promise<string>;
  isFullScreen: () => Promise<boolean>;
  onFullScreenChange: (callback: (isFullScreen: boolean) => void) => () => void;
}

declare global {
  interface Window {
    electronAPI: ElectronAPI;
    flaskPort?: number;
  }
}
