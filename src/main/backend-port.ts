import { findAvailablePort } from './port-manager';

const DEFAULT_PORT_RANGE_END = 5100;
const EXTENDED_PORT_RANGE_SIZE = 100;

export const resolveBackendPort = async (
  defaultPort: number
): Promise<number> => {
  try {
    console.log(
      `Auto-detecting available port in range ${defaultPort}-${DEFAULT_PORT_RANGE_END}...`
    );
    const port = await findAvailablePort(defaultPort, DEFAULT_PORT_RANGE_END);
    console.log(`✓ Found available port: ${port}`);
    return port;
  } catch (error) {
    console.log(`⚠️  Auto-detection failed: ${error}`);
  }

  try {
    console.log(
      `Searching in extended range ${DEFAULT_PORT_RANGE_END + 1}-${DEFAULT_PORT_RANGE_END + EXTENDED_PORT_RANGE_SIZE}...`
    );
    const port = await findAvailablePort(
      DEFAULT_PORT_RANGE_END + 1,
      DEFAULT_PORT_RANGE_END + EXTENDED_PORT_RANGE_SIZE
    );
    console.log(`✓ Found port in extended range: ${port}`);
    return port;
  } catch (extendedError) {
    throw new Error(
      `CRITICAL: No available ports found in any range. This may indicate:\n1. Too many services running on localhost\n2. Firewall blocking port access\n3. System resource limitations\n\nTried ranges: ${defaultPort}-${DEFAULT_PORT_RANGE_END}, ${DEFAULT_PORT_RANGE_END + 1}-${DEFAULT_PORT_RANGE_END + EXTENDED_PORT_RANGE_SIZE}`
    );
  }
};
