import { ipcMain, shell, dialog, BrowserWindow } from 'electron';
import { showAboutDialog } from './menu';
import { assertAllowedIpcSender, validateExternalUrl } from './ipc-security';

/**
 * Register IPC handlers for renderer ↔ main communication.
 *
 * Every privileged handler (anything that accesses OS APIs, opens dialogs, or
 * exposes platform state) validates the sender via {@link assertAllowedIpcSender}
 * before proceeding. This prevents rogue renderers or compromised webviews from
 * invoking capabilities they should not have access to.
 */
export const registerIpcHandlers = (
  getMainWindow: () => BrowserWindow | null
): void => {
  // IPC handler for getting platform information
  ipcMain.handle('get-platform', (event) => {
    assertAllowedIpcSender(event, getMainWindow);
    return process.platform;
  });

  // IPC handler for getting fullscreen state
  ipcMain.handle('get-fullscreen', (event) => {
    assertAllowedIpcSender(event, getMainWindow);
    const mainWindow = getMainWindow();
    if (mainWindow && !mainWindow.isDestroyed()) {
      return mainWindow.isFullScreen();
    }
    return false;
  });

  // IPC handler for opening external URLs in default browser
  ipcMain.handle('open-external-url', async (event, url: unknown) => {
    assertAllowedIpcSender(event, getMainWindow);

    const validation = validateExternalUrl(url);
    if (!validation.valid) {
      if (typeof url === 'string') {
        console.warn(`Blocked attempt to open invalid URL: ${url}`);
      }
      return { success: false, error: validation.error };
    }

    try {
      await shell.openExternal(validation.url);
      return { success: true };
    } catch (error) {
      console.error(`Failed to open external URL: ${validation.url}`, error);
      return { success: false, error: 'Failed to open URL' };
    }
  });

  // IPC handler for showing About dialog
  ipcMain.handle('show-about-dialog', (event) => {
    assertAllowedIpcSender(event, getMainWindow);
    showAboutDialog();
  });

  // IPC handler for selecting a folder
  ipcMain.handle('dialog:select-folder', async (event) => {
    assertAllowedIpcSender(event, getMainWindow);

    try {
      const mainWindow = getMainWindow();
      const dialogOptions: Electron.OpenDialogOptions = {
        properties: ['openDirectory', 'createDirectory'],
        title: 'Select Calculations Directory',
        buttonLabel: 'Select',
      };
      const result =
        mainWindow && !mainWindow.isDestroyed()
          ? await dialog.showOpenDialog(mainWindow, dialogOptions)
          : await dialog.showOpenDialog(dialogOptions);

      if (result.canceled || result.filePaths.length === 0) {
        return { canceled: true, filePath: null };
      }

      return { canceled: false, filePath: result.filePaths[0] };
    } catch (error) {
      console.error('Failed to open folder selection dialog:', error);
      return { canceled: true, filePath: null, error: 'Dialog failed' };
    }
  });
};
