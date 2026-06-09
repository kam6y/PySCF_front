export interface ElectronAPI {
  // Backend port number from URL parameters (set in preload)
  backendPort: number | null;
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
  }
}
