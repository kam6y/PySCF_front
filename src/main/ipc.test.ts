const assert: typeof import('node:assert/strict') = require('node:assert/strict');

/**
 * Tests for the extracted confirmAndOpenExternal function (J6).
 *
 * Covers the three paths:
 * 1. User clicks Cancel (or dialog returns non-0) -> MUST NOT open
 * 2. User clicks Open (response === 0) -> opens the URL
 * 3. Dialog throws (window destroyed race) -> fail-closed cancel
 */

type MessageBoxOptions = {
  type?: string;
  buttons?: string[];
  defaultId?: number;
  cancelId?: number;
  title?: string;
  message?: string;
  detail?: string;
};

type ConfirmAndOpenExternalDeps = {
  showDialog: (
    options: MessageBoxOptions,
    parentWindow?: unknown
  ) => Promise<{ response: number }>;
  openExternal: (url: string) => Promise<void>;
  parentWindow?: unknown;
};

type OpenExternalResult = { success: true } | { success: false; error: string };

type IpcModule = {
  confirmAndOpenExternal: (
    url: string,
    deps: ConfirmAndOpenExternalDeps
  ) => Promise<OpenExternalResult>;
};

const loadIpc = (): IpcModule => {
  return require('./ipc') as IpcModule;
};

// ============================================================
// J6: confirmAndOpenExternal tests
// ============================================================

const testConfirmAndOpen_userCancels_doesNotOpen = async (): Promise<void> => {
  const { confirmAndOpenExternal } = loadIpc();

  let openCalled = false;
  const result = await confirmAndOpenExternal('https://example.com/', {
    showDialog: async () => ({ response: 1 }), // Cancel button (index 1)
    openExternal: async () => { openCalled = true; },
  });

  assert.equal(result.success, false, 'Must return failure on cancel');
  assert.equal(openCalled, false, 'openExternal must NOT be called on cancel');
  if (!result.success) {
    assert.equal(result.error, 'User cancelled');
  }
};

const testConfirmAndOpen_userConfirms_opensUrl = async (): Promise<void> => {
  const { confirmAndOpenExternal } = loadIpc();

  let openedUrl = '';
  const result = await confirmAndOpenExternal('https://example.com/', {
    showDialog: async () => ({ response: 0 }), // Open button (index 0)
    openExternal: async (url) => { openedUrl = url; },
  });

  assert.equal(result.success, true, 'Must return success on confirm');
  assert.equal(openedUrl, 'https://example.com/', 'Must open the correct URL');
};

const testConfirmAndOpen_dialogThrows_failsClosed = async (): Promise<void> => {
  const { confirmAndOpenExternal } = loadIpc();

  let openCalled = false;
  const result = await confirmAndOpenExternal('https://example.com/', {
    showDialog: async () => { throw new Error('Window destroyed'); },
    openExternal: async () => { openCalled = true; },
  });

  assert.equal(result.success, false, 'Must return failure on dialog throw');
  assert.equal(openCalled, false, 'openExternal must NOT be called on dialog failure');
  if (!result.success) {
    assert.equal(result.error, 'Dialog failed');
  }
};

const testConfirmAndOpen_openExternalThrows_returnsFailure = async (): Promise<void> => {
  const { confirmAndOpenExternal } = loadIpc();

  const result = await confirmAndOpenExternal('https://example.com/', {
    showDialog: async () => ({ response: 0 }), // User confirms
    openExternal: async () => { throw new Error('Shell error'); },
  });

  assert.equal(result.success, false, 'Must return failure when openExternal throws');
  if (!result.success) {
    assert.equal(result.error, 'Failed to open URL');
  }
};

const testConfirmAndOpen_dialogOptions_correct = async (): Promise<void> => {
  const { confirmAndOpenExternal } = loadIpc();

  let capturedOptions: MessageBoxOptions | null = null;
  await confirmAndOpenExternal('https://example.com/path', {
    showDialog: async (options) => {
      capturedOptions = options;
      return { response: 1 }; // Cancel
    },
    openExternal: async () => {},
  });

  assert.ok(capturedOptions !== null, 'Dialog must be called');
  const opts = capturedOptions as MessageBoxOptions;
  assert.equal(opts.defaultId, 1, 'defaultId must be Cancel (1)');
  assert.equal(opts.cancelId, 1, 'cancelId must be Cancel (1)');
  assert.equal(opts.detail, 'https://example.com/path', 'detail must show the URL');
};

// ============================================================
// Runner
// ============================================================

const run = async (): Promise<void> => {
  await testConfirmAndOpen_userCancels_doesNotOpen();
  await testConfirmAndOpen_userConfirms_opensUrl();
  await testConfirmAndOpen_dialogThrows_failsClosed();
  await testConfirmAndOpen_openExternalThrows_returnsFailure();
  await testConfirmAndOpen_dialogOptions_correct();

  console.log('ipc tests passed (5 tests)');
};

run().catch((err) => {
  console.error(err);
  process.exit(1);
});

export {};
