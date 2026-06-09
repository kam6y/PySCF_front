const assert: typeof import('node:assert/strict') = require('node:assert/strict');

const MAX_WAIT_TICKS = 20;

type StreamChatWithAgent = typeof import('./agent').streamChatWithAgent;
type AbortMode = 'onerror' | 'reject-only';

interface MockFetchEventSourceOptions {
  signal: AbortSignal;
  onerror?: (error: unknown) => void;
}

const Module = require('node:module') as {
  _load: (...args: unknown[]) => unknown;
};
const originalLoad = Module._load;
const pendingFetches: Promise<void>[] = [];
let abortMode: AbortMode = 'onerror';

const createAbortError = (): Error => {
  const error = new Error('The operation was aborted.');
  error.name = 'AbortError';
  return error;
};

const rejectWithAbort = (
  options: MockFetchEventSourceOptions,
  reject: (reason?: unknown) => void
): void => {
  const abortError = createAbortError();

  if (abortMode === 'onerror') {
    try {
      options.onerror?.(abortError);
    } catch (error) {
      reject(error);
      return;
    }
  }

  reject(abortError);
};

const mockFetchEventSource = async (
  _url: string,
  options: MockFetchEventSourceOptions
): Promise<void> => {
  const pendingFetch = new Promise<void>((_resolve, reject) => {
    const abort = (): void => rejectWithAbort(options, reject);

    if (options.signal.aborted) {
      queueMicrotask(abort);
      return;
    }

    options.signal.addEventListener('abort', abort, { once: true });
  });

  pendingFetches.push(pendingFetch);
  return pendingFetch;
};

Module._load = (...args: unknown[]): unknown => {
  const [request] = args;

  if (request === '@microsoft/fetch-event-source') {
    return { fetchEventSource: mockFetchEventSource };
  }

  return originalLoad(...args);
};

(globalThis as any).window = {
  electronAPI: {
    backendPort: 5000,
  },
};

const { streamChatWithAgent } = require('./agent') as {
  streamChatWithAgent: StreamChatWithAgent;
};

const waitForFetchStart = async (): Promise<void> => {
  for (let tick = 0; tick < MAX_WAIT_TICKS; tick += 1) {
    if (pendingFetches.length > 0) {
      return;
    }

    await new Promise(resolve => setTimeout(resolve, 0));
  }

  throw new Error('fetchEventSource was not called.');
};

const waitForStreamToSettle = async (): Promise<void> => {
  await Promise.allSettled(pendingFetches);
  await new Promise(resolve => setTimeout(resolve, 0));
};

const resetMockState = (mode: AbortMode): void => {
  pendingFetches.length = 0;
  abortMode = mode;
};

const runCancelledStream = async (
  mode: AbortMode
): Promise<{ closeCount: number; errors: Error[] }> => {
  resetMockState(mode);

  const errors: Error[] = [];
  let closeCount = 0;
  const cancel = streamChatWithAgent('hello', [], null, {
    onMessage: () => undefined,
    onClose: () => {
      closeCount += 1;
    },
    onError: error => {
      errors.push(error);
    },
  });

  await waitForFetchStart();
  cancel();
  await waitForStreamToSettle();

  return { closeCount, errors };
};

const testCancellationIgnoresFetchEventSourceOnError =
  async (): Promise<void> => {
    const { closeCount, errors } = await runCancelledStream('onerror');

    assert.equal(closeCount, 0);
    assert.deepEqual(errors, []);
  };

const testCancellationIgnoresRejectedAbortError = async (): Promise<void> => {
  const { closeCount, errors } = await runCancelledStream('reject-only');

  assert.equal(closeCount, 0);
  assert.deepEqual(errors, []);
};

const run = async (): Promise<void> => {
  await testCancellationIgnoresFetchEventSourceOnError();
  await testCancellationIgnoresRejectedAbortError();
  console.log('agent SSE cancellation tests passed');
};

run()
  .catch(error => {
    console.error(error);
    process.exitCode = 1;
  })
  .finally(() => {
    Module._load = originalLoad;
  });

export {};
