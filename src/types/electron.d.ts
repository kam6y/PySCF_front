export interface ElectronAPI {
  closeWindow: () => Promise<void>;
  minimizeWindow: () => Promise<void>;
  maximizeWindow: () => Promise<void>;
  isWindowMaximized: () => Promise<boolean>;
  onWindowStateChange: (callback: (isMaximized: boolean) => void) => void;
  removeAllListeners: (channel: string) => void;
  onSetFlaskPort: (callback: (port: number) => void) => () => void;
  getFlaskPort: () => Promise<number | null>;
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
