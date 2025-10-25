// src/preload.ts

import { contextBridge, ipcRenderer } from 'electron';

// Expose protected methods that allow the renderer process to use
// the ipcRenderer without exposing the entire object
contextBridge.exposeInMainWorld('electronAPI', {
  // We can also expose variables, not just functions
  onSetFlaskPort: (callback: (port: number) => void) => {
    const handler = (_event: any, port: number) => callback(port);
    ipcRenderer.on('set-flask-port', handler);
    // Return a cleanup function
    return () => {
      ipcRenderer.removeListener('set-flask-port', handler);
    };
  },
  getFlaskPort: () => ipcRenderer.invoke('get-flask-port'),
  openExternalUrl: (url: string) =>
    ipcRenderer.invoke('open-external-url', url),
  showAboutDialog: () => ipcRenderer.invoke('show-about-dialog'),
  selectFolder: () =>
    ipcRenderer.invoke('dialog:select-folder') as Promise<{
      canceled: boolean;
      filePath: string | null;
      error?: string;
    }>,
  getPlatform: () => ipcRenderer.invoke('get-platform'),
  isFullScreen: () => ipcRenderer.invoke('get-fullscreen'),
  onFullScreenChange: (callback: (isFullScreen: boolean) => void) => {
    const handler = (_event: any, isFullScreen: boolean) =>
      callback(isFullScreen);
    ipcRenderer.on('fullscreen-changed', handler);
    // Return a cleanup function
    return () => {
      ipcRenderer.removeListener('fullscreen-changed', handler);
    };
  },
});
