const assert: typeof import('node:assert/strict') = require('node:assert/strict');
const { EventEmitter } = require('node:events') as typeof import('node:events');

const Module = require('node:module') as {
  _load: (...args: unknown[]) => unknown;
};
const originalLoad = Module._load;

const electronState = {
  quitCount: 0,
  errorBoxes: [] as Array<{ title: string; message: string }>,
};

Module._load = (...args: unknown[]): unknown => {
  const [request] = args;

  if (request === 'electron') {
    return {
      app: {
        isPackaged: false,
        quit: () => {
          electronState.quitCount += 1;
        },
      },
      dialog: {
        showErrorBox: (title: string, message: string) => {
          electronState.errorBoxes.push({ title, message });
        },
      },
    };
  }

  return originalLoad(...args);
};

type MockProcess = InstanceType<typeof EventEmitter> & {
  stdout: InstanceType<typeof EventEmitter>;
  stderr: InstanceType<typeof EventEmitter>;
  exitCode: number | null;
  killed: boolean;
  killCalls: string[];
  kill: (signal: string) => void;
};

const createMockProcess = (): MockProcess => {
  const proc = new EventEmitter() as MockProcess;
  proc.stdout = new EventEmitter();
  proc.stderr = new EventEmitter();
  proc.exitCode = null;
  proc.killed = false;
  proc.killCalls = [];
  proc.kill = (signal: string): void => {
    proc.killCalls.push(signal);
  };
  return proc;
};

const createContext = () => ({
  pythonExecutablePath: '/mock/conda/bin/python',
  pythonPath: '/mock/src/python',
  serverPort: 5050,
});

const resetElectronState = (): void => {
  electronState.quitCount = 0;
  electronState.errorBoxes = [];
};

const loadBackendProcess = (): typeof import('./backend-process') => {
  return require('./backend-process') as typeof import('./backend-process');
};

const asChildProcess = (
  proc: MockProcess
): import('child_process').ChildProcess => {
  return proc as unknown as import('child_process').ChildProcess;
};

const testStopMarksExpectedCloseWithoutDialog = (): void => {
  resetElectronState();
  const { BackendProcessController } = loadBackendProcess();
  const controller = new BackendProcessController();
  const proc = createMockProcess();

  controller.attach(asChildProcess(proc), createContext(), error => {
    throw error;
  });

  assert.equal(controller.isRunning(), true);

  controller.stop();
  proc.emit('close', 0, null);

  assert.deepEqual(proc.killCalls, ['SIGTERM']);
  assert.equal(controller.isRunning(), false);
  assert.equal(electronState.quitCount, 0);
  assert.deepEqual(electronState.errorBoxes, []);
};

const testUnexpectedCloseShowsDialogAndQuits = (): void => {
  resetElectronState();
  const { BackendProcessController } = loadBackendProcess();
  const controller = new BackendProcessController();
  const proc = createMockProcess();

  controller.attach(asChildProcess(proc), createContext(), error => {
    throw error;
  });
  proc.emit('close', 1, null);

  assert.equal(electronState.quitCount, 1);
  assert.equal(electronState.errorBoxes.length, 1);
  assert.equal(electronState.errorBoxes[0].title, 'Backend Process Error');
  assert.match(electronState.errorBoxes[0].message, /unexpectedly stopped/);
};

const run = (): void => {
  try {
    testStopMarksExpectedCloseWithoutDialog();
    testUnexpectedCloseShowsDialogAndQuits();
    console.log('backend process tests passed');
  } finally {
    Module._load = originalLoad;
  }
};

run();

export {};
