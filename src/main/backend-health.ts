import http from 'node:http';
import { updateSplashStatus } from './splash-window-manager';

export type BackendHealthDiagnosticOptions = {
  isPackaged: boolean;
  port: number;
  url: string;
  retries: number;
};

export type BackendHealthCheckOptions = {
  port: number;
  authToken: string;
  retries?: number;
  delay?: number;
  isPackaged: boolean;
};

export const buildBackendHealthDiagnosticMessage = ({
  isPackaged,
  port,
  url,
  retries,
}: BackendHealthDiagnosticOptions): string => {
  if (isPackaged) {
    return `Python backend failed to start after ${retries} attempts.\n\nDiagnostic information:\n- Port: ${port}\n- Health endpoint: ${url}\n\nThis may indicate:\n1. Bundled Python environment is corrupted\n2. Port ${port} is blocked by firewall\n3. Python dependencies are missing\n\nPlease report this issue with the console output.`;
  }

  return `Python backend failed to start after ${retries} attempts.\n\nDiagnostic information:\n- Port: ${port}\n- Health endpoint: ${url}\n- Environment: Development mode\n\nTroubleshooting steps:\n1. Check if conda environment 'pyscf-env' is activated\n2. Verify all dependencies are installed: conda env create -f .github/environment.yml\n3. Test the FastAPI backend manually: cd src/python && python app.py\n4. Check if port ${port} is available\n\nFor more details, see CLAUDE.md`;
};

export const checkBackendHealth = ({
  port,
  authToken,
  retries = 20,
  delay = 500,
  isPackaged,
}: BackendHealthCheckOptions): Promise<void> => {
  return new Promise((resolve, reject) => {
    let attempts = 0;
    const url = `http://127.0.0.1:${port}/health`;
    const options = {
      headers: {
        'X-Auth-Token': authToken,
      },
    };

    const rejectWithDiagnostic = (): void => {
      const diagnosticMessage = buildBackendHealthDiagnosticMessage({
        isPackaged,
        port,
        url,
        retries,
      });
      reject(new Error(diagnosticMessage));
    };

    const handleFailedAttempt = (message: string): void => {
      attempts += 1;
      console.log(message);
      updateSplashStatus('health-check', 'Waiting for server...', attempts);

      if (attempts >= retries) {
        clearInterval(interval);
        rejectWithDiagnostic();
      }
    };

    const interval = setInterval(() => {
      http
        .get(url, options, res => {
          if (res.statusCode === 200) {
            clearInterval(interval);
            console.log('Python server is healthy.');
            resolve();
            return;
          }

          res.resume();
          handleFailedAttempt(
            `Health check attempt ${attempts + 1}/${retries} failed with status: ${res.statusCode}`
          );
        })
        .on('error', _err => {
          handleFailedAttempt(
            `Health check attempt ${attempts + 1}/${retries} failed for ${url}`
          );
        });
    }, delay);
  });
};
