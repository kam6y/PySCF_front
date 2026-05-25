import type { ChildProcess } from 'child_process';
import { app, dialog } from 'electron';

export type BackendProcessContext = {
  pythonExecutablePath: string;
  pythonPath: string;
  serverPort: number;
};

export const attachBackendOutputHandlers = (
  proc: ChildProcess,
  serverPort: number
): void => {
  proc.stdout?.on('data', data => {
    const output = data.toString().trim();
    console.log(`[PYTHON STDOUT] ${output}`);
    if (output.includes('Starting gunicorn')) {
      console.log('✓ Gunicorn is starting up...');
    }
    if (output.includes('Listening at:')) {
      console.log('✓ Server is listening for connections');
    }
    if (output.includes('Booting worker')) {
      console.log('✓ Gunicorn worker is starting...');
    }
    if (output.includes('Application object must be callable')) {
      console.log('✗ CRITICAL: Python/FastAPI backend object error detected');
    }
    if (
      output.includes('ModuleNotFoundError') ||
      output.includes('ImportError')
    ) {
      console.log(`✗ CRITICAL: Python import error detected - ${output}`);
    }
  });

  proc.stderr?.on('data', data => {
    const errorOutput = data.toString().trim();
    console.log(`[PYTHON STDERR] ${errorOutput}`);
    if (errorOutput.includes('ModuleNotFoundError')) {
      console.log(`✗ CRITICAL: Missing Python module - ${errorOutput}`);
    }
    if (errorOutput.includes('ImportError')) {
      console.log(`✗ CRITICAL: Python import error - ${errorOutput}`);
    }
    if (errorOutput.includes('gunicorn')) {
      console.log(`⚠️  Gunicorn-related error - ${errorOutput}`);
    }
    if (errorOutput.toLowerCase().includes('fastapi')) {
      console.log(`⚠️  FastAPI-related error - ${errorOutput}`);
    }
    if (errorOutput.includes('Address already in use')) {
      console.log(`✗ CRITICAL: Port ${serverPort} is already in use`);
    }
    if (
      errorOutput.includes('[CRITICAL]') ||
      errorOutput.includes('CRITICAL')
    ) {
      console.log(`✗ CRITICAL ERROR FROM PYTHON: ${errorOutput}`);
    }
  });
};

export class BackendProcessController {
  private proc: ChildProcess | null = null;
  private isQuitting = false;
  private forceKillTimer: NodeJS.Timeout | null = null;

  public isRunning(): boolean {
    return (
      this.proc !== null && !this.proc.killed && this.proc.exitCode === null
    );
  }

  public attach(
    proc: ChildProcess,
    context: BackendProcessContext,
    reject: (err: Error) => void
  ): void {
    this.proc = proc;
    this.isQuitting = false;
    attachBackendOutputHandlers(proc, context.serverPort);
    this.attachLifecycleHandlers(proc, context, reject);
  }

  public stop(): void {
    const processToStop = this.proc;
    if (processToStop && processToStop.exitCode === null) {
      console.log('Stopping Python/FastAPI backend...');
      this.isQuitting = true;
      processToStop.kill('SIGTERM');

      if (this.forceKillTimer !== null) {
        clearTimeout(this.forceKillTimer);
      }
      this.forceKillTimer = setTimeout(() => {
        if (processToStop.exitCode === null) {
          console.log('Force killing Python server...');
          processToStop.kill('SIGKILL');
        }
        this.forceKillTimer = null;
      }, 5000);
    }
  }

  private attachLifecycleHandlers(
    proc: ChildProcess,
    context: BackendProcessContext,
    reject: (err: Error) => void
  ): void {
    const { pythonExecutablePath, pythonPath, serverPort } = context;

    proc.on('error', error => {
      const processError = error as NodeJS.ErrnoException;
      console.error(`✗ CRITICAL: Failed to start Python server process`);
      console.error(`Error details: ${processError.message}`);
      console.error(`Error code: ${processError.code || 'N/A'}`);
      console.error(`Error errno: ${processError.errno || 'N/A'}`);
      console.error(`Error syscall: ${processError.syscall || 'N/A'}`);
      console.error(`Python executable path: ${pythonExecutablePath}`);
      console.error(`Working directory: ${pythonPath}`);
      console.error(`Environment variables:`, {
        CONDA_DEFAULT_ENV: process.env.CONDA_DEFAULT_ENV,
        PATH: process.env.PATH?.split(':')
          .filter(p => p.includes('conda'))
          .slice(0, 3),
        PYTHONPATH: process.env.PYTHONPATH || 'Not set',
      });
      reject(error);
    });

    proc.on('close', (code, signal) => {
      console.log(`=== Python Server Process Terminated ===`);
      console.log(`Exit code: ${code}`);
      console.log(`Signal: ${signal || 'None'}`);
      console.log(`Was quitting: ${this.isQuitting}`);
      console.log(`Server port: ${serverPort}`);
      console.log(`Python executable: ${pythonExecutablePath}`);
      console.log(`Working directory: ${pythonPath}`);

      if (this.forceKillTimer !== null) {
        clearTimeout(this.forceKillTimer);
        this.forceKillTimer = null;
      }
      if (this.proc === proc) {
        this.proc = null;
      }

      if (!this.isQuitting) {
        const errorMessage = `The Python backend process has unexpectedly stopped (exit code: ${code})${signal ? `, signal: ${signal}` : ''}.\n\nDebugging Information:\n• Python executable: ${pythonExecutablePath}\n• Working directory: ${pythonPath}\n• Server port: ${serverPort}\n• Packaged mode: ${app.isPackaged}\n\nPlease check the console output for detailed error messages and restart the application.`;
        console.log(`✗ CRITICAL: Showing error dialog to user`);
        dialog.showErrorBox('Backend Process Error', errorMessage);
        app.quit();
      }
    });
  }
}
