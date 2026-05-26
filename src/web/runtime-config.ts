type BackendRuntimeConfig = {
  backendPort: number | null;
  backendBaseUrl: string;
};

const isValidBackendPort = (
  port: number | null | undefined
): port is number => {
  return (
    typeof port === 'number' &&
    Number.isInteger(port) &&
    port > 0 &&
    port < 65536
  );
};

const getBackendPort = (): number | null => {
  const backendPort = window.electronAPI?.backendPort ?? null;
  return isValidBackendPort(backendPort) ? backendPort : null;
};

export const getBackendRuntimeConfig = (): BackendRuntimeConfig => {
  const backendPort = getBackendPort();
  return {
    backendPort,
    backendBaseUrl: backendPort ? `http://127.0.0.1:${backendPort}` : '',
  };
};
