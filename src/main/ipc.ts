import { ipcMain, shell, dialog, BrowserWindow } from 'electron';
import { showAboutDialog } from './menu';

/**
 * IPCハンドラーを登録する関数
 */
export const registerIpcHandlers = (
  getMainWindow: () => BrowserWindow | null
): void => {
  // IPC handler for getting platform information
  ipcMain.handle('get-platform', () => process.platform);

  // IPC handler for getting fullscreen state
  ipcMain.handle('get-fullscreen', () => {
    const mainWindow = getMainWindow();
    if (mainWindow && !mainWindow.isDestroyed()) {
      return mainWindow.isFullScreen();
    }
    return false;
  });

  // IPC handler for opening external URLs in default browser
  ipcMain.handle('open-external-url', async (_event, url: string) => {
    try {
      // Security: Only allow http and https protocols
      const parsedUrl = new URL(url);
      if (parsedUrl.protocol !== 'http:' && parsedUrl.protocol !== 'https:') {
        console.warn(`Blocked attempt to open non-http(s) URL: ${url}`);
        return {
          success: false,
          error: 'Only HTTP and HTTPS URLs are allowed',
        };
      }

      console.log(`Opening external URL in default browser: ${url}`);
      await shell.openExternal(url);
      return { success: true };
    } catch (error) {
      console.error(`Failed to open external URL: ${url}`, error);
      return { success: false, error: String(error) };
    }
  });

  // IPC handler for showing About dialog
  ipcMain.handle('show-about-dialog', () => {
    showAboutDialog();
  });

  // IPC handler for selecting a folder
  ipcMain.handle('dialog:select-folder', async () => {
    try {
      const mainWindow = getMainWindow();
      const result = await dialog.showOpenDialog(
        mainWindow && !mainWindow.isDestroyed()
          ? mainWindow
          : (undefined as any),
        {
          properties: ['openDirectory', 'createDirectory'],
          title: 'Select Calculations Directory',
          buttonLabel: 'Select',
        }
      );

      if (result.canceled || result.filePaths.length === 0) {
        return { canceled: true, filePath: null };
      }

      return { canceled: false, filePath: result.filePaths[0] };
    } catch (error) {
      console.error('Failed to open folder selection dialog:', error);
      return { canceled: true, filePath: null, error: String(error) };
    }
  });
};
