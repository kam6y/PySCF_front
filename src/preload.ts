import { contextBridge, ipcRenderer } from "electron";

console.log("preloaded!");

// Expose window control APIs to the renderer process
contextBridge.exposeInMainWorld("electronAPI", {
  // Window controls
  closeWindow: () => ipcRenderer.invoke("window-close"),
  minimizeWindow: () => ipcRenderer.invoke("window-minimize"),
  maximizeWindow: () => ipcRenderer.invoke("window-maximize"),
  
  // Check if window is maximized
  isWindowMaximized: () => ipcRenderer.invoke("window-is-maximized"),
  
  // Listen to window state changes
  onWindowStateChange: (callback: (isMaximized: boolean) => void) => {
    ipcRenderer.on("window-state-changed", (event, isMaximized) => callback(isMaximized));
  },
  
  // Remove listeners
  removeAllListeners: (channel: string) => {
    ipcRenderer.removeAllListeners(channel);
  }
});