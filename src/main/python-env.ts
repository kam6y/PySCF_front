import path from 'node:path';
import fs from 'fs';
import { app } from 'electron';
import { execFile } from 'child_process';
import { promisify } from 'util';

const execFilePromise = promisify(execFile);

/**
 * 開発時conda環境のPythonパスを簡素化して検出する
 * @returns conda環境のPythonパス、見つからなければnull
 */
const detectCondaEnvironmentPath = async (): Promise<string | null> => {
  const envName = 'pyscf-env';

  // 1. 環境変数での指定をチェック
  const envPath = process.env.CONDA_ENV_PATH;
  if (envPath) {
    const pythonPath = path.join(envPath, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Using conda environment from CONDA_ENV_PATH: ${pythonPath}`);
      return pythonPath;
    }
  }

  // 2. conda info --base で環境パスを取得
  try {
    const { stdout } = await execFilePromise('conda', ['info', '--base'], {
      timeout: 3000,
      encoding: 'utf8',
    });
    const basePath = stdout.trim();
    const pythonPath = path.join(basePath, 'envs', envName, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Found conda environment via command: ${pythonPath}`);
      return pythonPath;
    }
  } catch (error) {
    console.log(
      'conda command unavailable in PATH, trying fallback locations...'
    );
  }

  // 3. 一般的なインストールパスを探索（フォールバック）
  const homeDir = app.getPath('home');
  const commonLocations = [
    // macOS / Linux
    path.join(homeDir, 'miniforge3'),
    path.join(homeDir, 'miniconda3'),
    path.join(homeDir, 'anaconda3'),
    path.join(homeDir, 'opt', 'miniforge3'),
    path.join(homeDir, 'opt', 'miniconda3'),
    path.join(homeDir, 'opt', 'anaconda3'),
    '/opt/miniconda3',
    '/opt/anaconda3',
    // Windows (typically handled by different path structure, but good to have)
    path.join(homeDir, 'Miniconda3'),
    path.join(homeDir, 'Anaconda3'),
    'C:\\ProgramData\\miniconda3',
    'C:\\ProgramData\\Anaconda3',
    path.join(homeDir, 'AppData', 'Local', 'Continuum', 'miniconda3'),
    path.join(homeDir, 'AppData', 'Local', 'Continuum', 'anaconda3'),
  ];

  for (const basePath of commonLocations) {
    const pythonPath = path.join(basePath, 'envs', envName, 'bin', 'python');
    if (fs.existsSync(pythonPath)) {
      console.log(`Found conda environment at fallback path: ${pythonPath}`);
      return pythonPath;
    }
  }

  console.log(
    `conda environment '${envName}' not found in PATH or common locations`
  );
  return null;
};

/**
 * Python環境のパスを簡素化して検出する統合関数
 * 1. パッケージ時: 同梱conda環境のみ
 * 2. 開発時: pyscf-env環境のみ
 * @returns Python環境のパス、見つからなければnull
 */
export const detectPythonEnvironmentPath = async (): Promise<string | null> => {
  console.log('=== Python Environment Detection (Simplified) ===');
  console.log(`Running in packaged mode: ${app.isPackaged}`);

  // 1. パッケージ時：同梱conda環境のみ
  if (app.isPackaged) {
    const bundledCondaPath = path.join(
      process.resourcesPath,
      'conda_env',
      'bin',
      'python'
    );
    const bundledGunicornPath = path.join(
      process.resourcesPath,
      'conda_env',
      'bin',
      'gunicorn'
    );

    console.log(`Checking bundled conda environment: ${bundledCondaPath}`);

    if (fs.existsSync(bundledCondaPath) && fs.existsSync(bundledGunicornPath)) {
      console.log(`✓ Using bundled conda environment: ${bundledCondaPath}`);
      return bundledCondaPath;
    } else {
      console.log(`✗ Bundled conda environment incomplete or missing`);
      console.log(`  - Python exists: ${fs.existsSync(bundledCondaPath)}`);
      console.log(`  - Gunicorn exists: ${fs.existsSync(bundledGunicornPath)}`);
      return null;
    }
  }

  // 2. 開発時：pyscf-env環境のみ
  console.log('Detecting development conda environment...');
  const condaPath = await detectCondaEnvironmentPath();
  if (condaPath) {
    console.log(`✓ Using conda environment: ${condaPath}`);
    return condaPath;
  }

  console.log('✗ No Python environment found');
  return null;
};

/**
 * Create a clean environment for the Python process using a whitelist approach.
 * This prevents user's local Python environment (pyenv, conda, venv) from interfering.
 *
 * @param condaBinDir - The bin directory of the conda environment to use
 * @param serverPort - The port number for the Python backend server
 * @param authToken - The authentication token for the Python backend server
 * @returns A clean environment object
 */
export const createCleanEnvironment = (
  condaBinDir: string,
  serverPort: number,
  authToken: string
): Record<string, string> => {
  // Whitelist of environment variables to pass through
  const ALLOWED_ENV_VARS = [
    // System basics
    'HOME',
    'USER',
    'TMPDIR',
    'SHELL',
    'TERM',
    'LANG',
    'LC_ALL',
    // Display / GUI
    'DISPLAY',
    'XAUTHORITY',
    // MacOS specific
    '__CF_USER_TEXT_ENCODING',
    // Dynamic Linker
    'LD_LIBRARY_PATH',
    'DYLD_LIBRARY_PATH',
    // Scientific Computing
    'OMP_NUM_THREADS',
    'MKL_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    'RDBASE',
    // Build / System
    'PKG_CONFIG_PATH',
  ];

  const cleanEnv: Record<string, string> = {};

  // 1. Pass through whitelisted variables
  for (const key of ALLOWED_ENV_VARS) {
    if (process.env[key] !== undefined) {
      cleanEnv[key] = process.env[key]!;
    }
  }

  // 2. Construct a clean PATH
  // We want to keep system paths but exclude any user-land Python paths
  const originalPath = process.env.PATH || '';
  const pathEntries = originalPath.split(':').filter(p => {
    // Exclude common Python environment paths
    // We want to avoid using the user's local python environments
    // but we should be careful not to exclude system paths that might be needed
    if (p.includes('/.pyenv/versions/')) return false;
    if (p.includes('/anaconda') && p.includes('/bin')) return false;
    if (p.includes('/miniconda') && p.includes('/bin')) return false;
    if (p.includes('virtualenvs')) return false;
    return true;
  });

  // Prepend our conda bin directory to ensure it takes precedence
  cleanEnv.PATH = `${condaBinDir}:${pathEntries.join(':')}`;

  // 3. Set application-specific variables
  cleanEnv.CONDA_DEFAULT_ENV = 'pyscf-env';
  cleanEnv.PYSCF_SERVER_PORT = String(serverPort);
  cleanEnv.PYSCF_PARENT_PID = String(process.pid);
  cleanEnv.PYSCF_AUTH_TOKEN = authToken;

  // Explicitly unset potentially conflicting Python variables
  // (They are already not in the whitelist, but this documents intent)
  // PYTHONPATH, PYTHONHOME, VIRTUAL_ENV are NOT copied.

  // 4. Set environment-specific variables
  if (app.isPackaged) {
    cleanEnv.PYSCF_RESOURCES_PATH = process.resourcesPath;
    cleanEnv.PYSCF_ENV = 'production';
    console.log(`Setting PYSCF_RESOURCES_PATH=${process.resourcesPath}`);
  } else {
    cleanEnv.PYSCF_ENV = 'development';
    console.log('Setting PYSCF_ENV=development');
  }

  // 5. Platform-specific fixes
  if (process.platform === 'darwin' && process.arch === 'arm64') {
    cleanEnv.OPENBLAS_CORETYPE = 'ARMV8';
    if (!cleanEnv.LC_ALL) cleanEnv.LC_ALL = 'C';
  }

  return cleanEnv;
};
