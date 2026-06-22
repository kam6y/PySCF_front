import { ipcMain, shell, dialog, BrowserWindow } from 'electron';
import { showAboutDialog } from './menu';
import { assertAllowedIpcSender, validateExternalUrl } from './ipc-security';
import { sanitizeForLog } from './app-protocol';

/**
 * Dependencies for the confirm-and-open-external flow, extracted for testability (J6).
 */
export type ConfirmAndOpenExternalDeps = {
  showDialog: (
    options: Electron.MessageBoxOptions,
    parentWindow?: BrowserWindow | null
  ) => Promise<{ response: number }>;
  openExternal: (url: string) => Promise<void>;
  parentWindow?: BrowserWindow | null;
};

/**
 * Result type for the open-external-url flow.
 */
export type OpenExternalResult = { success: true } | { success: false; error: string };

/**
 * Show confirmation dialog and open external URL if user confirms.
 *
 * Extracted from the IPC handler so the three paths (cancel, confirm, error)
 * can be unit-tested with stubs (J6).
 *
 * @param url - The validated URL string to open.
 * @param deps - Injectable dependencies for dialog and shell.openExternal.
 */
export const confirmAndOpenExternal = async (
  url: string,
  deps: ConfirmAndOpenExternalDeps
): Promise<OpenExternalResult> => {
  const dialogOptions: Electron.MessageBoxOptions = {
    type: 'question',
    buttons: ['Open in Browser', 'Cancel'],
    defaultId: 1,
    cancelId: 1,
    title: 'Open External Link',
    message: 'Do you want to open this link in your default browser?',
    detail: url,
  };

  // H4: Wrap dialog call in try/catch to handle TOCTOU race — parentWindow
  // can be destroyed between the isDestroyed() check and the await.
  // Fail-closed: treat any dialog error as cancel.
  let response: number;
  try {
    const result = await deps.showDialog(dialogOptions, deps.parentWindow);
    response = result.response;
  } catch {
    // Window destroyed during dialog or other dialog failure — treat as cancel
    return { success: false, error: 'Dialog failed' };
  }

  if (response !== 0) {
    return { success: false, error: 'User cancelled' };
  }

  try {
    await deps.openExternal(url);
    return { success: true };
  } catch (error) {
    // J9: Log sanitized URL and error message only (not raw error object)
    const errMsg = error instanceof Error ? error.message : 'Unknown error';
    console.error(
      `Failed to open external URL: ${sanitizeForLog(url)}: ${sanitizeForLog(errMsg)}`
    );
    return { success: false, error: 'Failed to open URL' };
  }
};

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
  ipcMain.handle('get-platform', (event) => {
    assertAllowedIpcSender(event, getMainWindow);
    return process.platform;
  });

  ipcMain.handle('get-fullscreen', (event) => {
    assertAllowedIpcSender(event, getMainWindow);
    const mainWindow = getMainWindow();
    if (mainWindow && !mainWindow.isDestroyed()) {
      return mainWindow.isFullScreen();
    }
    return false;
  });

  // IPC handler for opening external URLs in default browser.
  // Shows an interstitial confirmation dialog so the user sees the full
  // destination before anything is opened (M-001 mitigation).
  ipcMain.handle('open-external-url', async (event, url: unknown) => {
    assertAllowedIpcSender(event, getMainWindow);

    const validation = validateExternalUrl(url);
    if (!validation.valid) {
      if (typeof url === 'string') {
        // J9: Use sanitizeForLog for consistent truncation + control-char stripping
        console.warn(`Blocked attempt to open invalid URL: ${sanitizeForLog(url)}`);
      }
      return { success: false, error: validation.error };
    }

    const parentWindow = getMainWindow();
    return confirmAndOpenExternal(validation.url, {
      showDialog: (options, parent) => {
        return parent && !parent.isDestroyed()
          ? dialog.showMessageBox(parent, options)
          : dialog.showMessageBox(options);
      },
      openExternal: (u) => shell.openExternal(u),
      parentWindow,
    });
  });

  ipcMain.handle('show-about-dialog', (event) => {
    assertAllowedIpcSender(event, getMainWindow);
    showAboutDialog();
  });

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
      // J9: Log error message only, not raw error object
      const errMsg = error instanceof Error ? error.message : 'Unknown error';
      console.error(`Failed to open folder selection dialog: ${errMsg}`);
      return { canceled: true, filePath: null, error: 'Dialog failed' };
    }
  });
};
