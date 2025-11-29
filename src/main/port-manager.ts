import net from 'net';

/**
 * ポート検出を行う統一関数
 */
export const findAvailablePort = async (
  startPort: number,
  endPort: number
): Promise<number> => {
  const attemptedPorts: number[] = [];

  for (let port = startPort; port <= endPort; port++) {
    attemptedPorts.push(port);
    try {
      await new Promise((resolve, reject) => {
        const server = net.createServer();
        server.listen(port, '127.0.0.1', () => {
          server.close(resolve);
        });
        server.on('error', reject);
      });
      return port;
    } catch (error) {
      // Port is in use, try next
      if (endPort - startPort <= 5) {
        // 少ない範囲の場合は詳細ログを出力
        console.log(`  Port ${port} is in use, trying next...`);
      }
      continue;
    }
  }

  const rangeSize = endPort - startPort + 1;
  const errorDetails =
    rangeSize <= 10
      ? ` (tried: ${attemptedPorts.join(', ')})`
      : ` (checked ${rangeSize} ports)`;

  throw new Error(
    `No available port found in range ${startPort}-${endPort}${errorDetails}`
  );
};
