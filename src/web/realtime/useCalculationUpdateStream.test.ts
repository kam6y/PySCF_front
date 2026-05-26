const assert: typeof import('node:assert/strict') = require('node:assert/strict');

type ParseCalculationUpdateStreamEvent =
  typeof import('./useCalculationUpdateStream').parseCalculationUpdateStreamEvent;

(globalThis as any).window = {
  electronAPI: {
    backendPort: 5000,
    getAuthToken: async () => null,
  },
};

const { parseCalculationUpdateStreamEvent } =
  require('./useCalculationUpdateStream') as {
    parseCalculationUpdateStreamEvent: ParseCalculationUpdateStreamEvent;
  };

const testCalculationUpdateReturnsTypeAndCalculation = (): void => {
  const calculation = {
    id: 'calc-1',
    name: 'Water',
    status: 'running',
    updatedAt: '2026-05-26T01:23:45Z',
  };

  const parsed = parseCalculationUpdateStreamEvent(
    JSON.stringify({
      type: 'calculation_update',
      payload: { calculation },
    })
  );

  assert.equal(parsed.type, 'calculation_update');
  assert.deepEqual(parsed.calculation, calculation);
};

const testHeartbeatReturnsDefaultWithoutCalculation = (): void => {
  const parsed = parseCalculationUpdateStreamEvent(
    JSON.stringify({ type: 'heartbeat' })
  );

  assert.equal(parsed.type, 'heartbeat');
  assert.equal(parsed.calculation, undefined);
};

const testErrorReturnsErrorMessage = (): void => {
  const parsed = parseCalculationUpdateStreamEvent(
    JSON.stringify({
      type: 'error',
      payload: { message: 'stream failed' },
    })
  );

  assert.equal(parsed.type, 'error');
  assert.equal(parsed.errorMessage, 'stream failed');
};

const run = (): void => {
  testCalculationUpdateReturnsTypeAndCalculation();
  testHeartbeatReturnsDefaultWithoutCalculation();
  testErrorReturnsErrorMessage();
  console.log('calculation update SSE parser tests passed');
};

run();

export {};
