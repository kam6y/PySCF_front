export interface ElectronAPI {
  closeWindow: () => Promise<void>;
  minimizeWindow: () => Promise<void>;
  maximizeWindow: () => Promise<void>;
  isWindowMaximized: () => Promise<boolean>;
  onWindowStateChange: (callback: (isMaximized: boolean) => void) => void;
  removeAllListeners: (channel: string) => void;
}

declare global {
  interface Window {
    electronAPI: ElectronAPI;
  }
}